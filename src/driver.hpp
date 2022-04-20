#pragma once

#include <cstddef>
#include <cstring>
#include <span>
#include <string_view>
#include <vector>

#include <filesystem>
#include <fstream>
#include <sol/sol.hpp>

#include <legion.h>
#include <new>

#include <nlopt.hpp>
#include <shoccs.hpp>

static_assert(sizeof(std::byte) == sizeof(char));

namespace fs = std::filesystem;

namespace detail
{

template <typename C>
class driver_
{
    C buf;
    std::optional<sol::state> lua_;

    char* script_string()
    {
        int dims = opt_dims();
        return (char*)buf.data() + (dims + 1) * sizeof(double);
    }

    int opt_dims() const { return *(const int*)buf.data(); }

    double* dims_data() { return (double*)buf.data() + 1; }

    sol::state& lua_state()
    {
        if (!lua_) {
            lua_ = sol::state{};
            lua_->open_libraries(
                sol::lib::base, sol::lib::string, sol::lib::package, sol::lib::math);
            lua_->safe_script(script_string());
        }
        return *lua_;
    }

public:
    driver_(const std::string& file) : lua_{sol::state{}}
    {
        // need to grab dims from file for allocation
        auto& lua = *lua_;
        lua.open_libraries(
            sol::lib::base, sol::lib::string, sol::lib::package, sol::lib::math);
        lua.script_file(file);

        int dims = lua["NLopt"]["dims"];

        auto sz = fs::file_size(file);
        auto dim_sz = (dims + 1) * sizeof(double);
        buf.resize(dim_sz + sz + 1);
        *(int*)buf.data() = dims;

        char* file_buf = (char*)(buf.data()) + dim_sz;
        std::ifstream is{file};
        is.read(file_buf, sz);
        file_buf[sz] = '\0';

        std::cout << "buf_size = " << buf.size() << '\n';
    }

    driver_(std::byte* buf, size_t buf_size) : buf{buf, buf_size} {}

    static driver_<std::span<std::byte>> from_task(const Legion::Task* task)
    {
        if (task->local_arglen > 0)
            return {(std::byte*)task->local_args, task->local_arglen};
        else
            return {(std::byte*)task->args, task->arglen};
    }

    nlopt::opt opt(Legion::Logger& log)
    {
        auto& lua = lua_state();

        nlopt::opt opt{};
        sol::table t = lua["NLopt"];

        std::string algorithm = t["algorithm"];
        int dims = t["dims"];

        if (algorithm == "LN_COBYLA") {
            opt = nlopt::opt(nlopt::LN_COBYLA, dims);
        } else {
            log.fatal("unknown nlopt algorithm");
        }

        sol::optional<double> xtol_rel = t["xtol_rel"];
        sol::optional<double> xtol_abs = t["xtol_abs"];
        sol::optional<int> maxeval = t["maxeval"];
        sol::optional<double> initial_step = t["initial_step"];

        if (xtol_rel) opt.set_xtol_rel(*xtol_rel);
        if (xtol_abs) opt.set_xtol_abs(*xtol_abs);
        if (maxeval) opt.set_maxeval(*maxeval);
        if (initial_step) opt.set_initial_step(*initial_step);

        return opt;
    }

    void set_data(std::span<const double> x)
    {
        memcpy(dims_data(), x.data(), x.size_bytes());
    }

    std::vector<double> guess() { return std::vector<double>(opt_dims()); }
    std::vector<double> params()
    {
        int dims = opt_dims();
        auto x = std::vector<double>(dims);
        memcpy(x.data(), dims_data(), dims * sizeof(double));
        return x;
    }

    double run(int i)
    {
        auto& lua = lua_state();
        sol::table t = lua["Simulations"][1];
        sol::function set_values = t["set_values"];
        sol::function result = t["result"];
        auto x = params();

        set_values(t, i + 1, x);
        auto r = ccs::simulation_run(t["simulations"][i + 1]);
        assert(r);
        return result(t, *r);
    }

    double result(std::span<const double> res) {
        auto& lua = lua_state();
        sol::table t = lua["Simulations"][1];
        return t["aggregate"](t, res);
    }

    int simulation_size()
    {
        auto& lua = lua_state();
        // really need to abstract this
        sol::table t = lua["Simulations"][1]["simulations"];
        assert(t.valid());
        return t.size();
    }

    operator Legion::TaskArgument() const { return {buf.data(), buf.size()}; }
};

} // namespace detail

using driver = detail::driver_<std::vector<std::byte>>;
using driver_span = detail::driver_<std::span<std::byte>>;
