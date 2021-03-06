#pragma once

#include <cstddef>
#include <cstring>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

#include <filesystem>
#include <fstream>
#include <functional>
#include <random>

#include <sol/sol.hpp>

#include <legion.h>
#include <new>

#include <nlopt.hpp>
#include <shoccs.hpp>

static_assert(sizeof(std::byte) == sizeof(char));
static_assert(2 * sizeof(int) <= sizeof(double));

namespace fs = std::filesystem;

namespace detail
{

//
// need to make the member function explicit rather than a template parameter
// since the nlopt member functions form an overload set
template <typename T>
inline void set_opt(const sol::optional<T>& v,
                    void (nlopt::opt::*fn)(T),
                    nlopt::opt& o,
                    std::optional<nlopt::opt>& local_o)
{
    if (v) {
        std::invoke(fn, o, *v);
        if (local_o) std::invoke(fn, *local_o, *v);
    }
}

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
    int task_index_() { return *((int*)buf.data() + 1); }

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
        if (task->local_arglen > task->arglen)
            return {(std::byte*)task->local_args, task->local_arglen};
        else
            return {(std::byte*)task->args, task->arglen};
    }

    static driver_<std::span<std::byte>> from_task(int& i, const Legion::Task* task)
    {
        i = *((const int*)task->args);
        return {(std::byte*)task->local_args, task->local_arglen};
    }

    nlopt::opt opt(Legion::Logger& log)
    {
        auto& lua = lua_state();

        nlopt::opt opt{};
        std::optional<nlopt::opt> local_opt{};

        sol::table t = lua["NLopt"];

        std::string algorithm = t["algorithm"];
        int dims = t["dims"];

        if (algorithm == "LN_COBYLA") {
            opt = nlopt::opt(nlopt::LN_COBYLA, dims);
        } else if (algorithm == "LN_SBPLX") {
            local_opt = nlopt::opt(nlopt::LN_SBPLX, dims);
            opt = nlopt::opt(nlopt::AUGLAG, dims);
        } else {
            log.fatal("unknown nlopt algorithm");
        }

        sol::optional<double> xtol_rel = t["xtol_rel"];
        sol::optional<double> xtol_abs = t["xtol_abs"];
        sol::optional<double> ftol_rel = t["ftol_rel"];
        sol::optional<double> ftol_abs = t["ftol_abs"];
        sol::optional<int> maxeval = t["maxeval"];
        sol::optional<double> initial_step = t["initial_step"];

        // What are the implications of maxeval and initial_step for the local_opt?
        // Will this launch a bunch of very expensive local optimizations?
        set_opt(xtol_rel, &nlopt::opt::set_xtol_rel, opt, local_opt);
        set_opt(xtol_abs, &nlopt::opt::set_xtol_abs, opt, local_opt);
        set_opt(ftol_rel, &nlopt::opt::set_ftol_rel, opt, local_opt);
        set_opt(ftol_abs, &nlopt::opt::set_ftol_abs, opt, local_opt);
        set_opt(initial_step, &nlopt::opt::set_initial_step, opt, local_opt);
        set_opt(maxeval, &nlopt::opt::set_maxeval, opt, local_opt);

        if (local_opt) opt.set_local_optimizer(*local_opt);

        return opt;
    }

    void set_data(std::span<const double> x)
    {
        memcpy(dims_data(), x.data(), x.size_bytes());
    }

    int& task_index() { return *((int*)buf.data() + 1); }

    std::vector<double> guess()
    {
        std::default_random_engine u{};
        std::random_device rd{};
        u.seed(rd());
        std::uniform_real_distribution<> d{-1.0, 1.0};
        auto v = std::vector<double>(opt_dims());
        for (auto&& x : v) x = d(u);
        return v;
    }
    std::vector<double> params()
    {
        int dims = opt_dims();
        auto x = std::vector<double>(dims);
        memcpy(x.data(), dims_data(), dims * sizeof(double));
        return x;
    }

    double run(int idx, int i)
    {
        auto& lua = lua_state();
        sol::table t = lua["Simulations"][idx + 1];
        sol::function set_values = t["set_values"];
        sol::function result = t["result"];
        auto x = params();

        set_values(t, i + 1, x);
        auto r = ccs::simulation_run(t["simulations"][i + 1]);
        assert(r);
        return result(t, *r);
    }

    double constraint()
    {
        auto& lua = lua_state();
        sol::table t = lua["Constraints"][1];

        sol::function set_values = t["set_values"];
        sol::function result = t["result"];
        sol::table sims = t["simulations"];

        auto x = params();
        std::vector<double> r(sims.size());

        for (int i = 0; i < sims.size(); i++) {
            set_values(t, i + 1, x);
            r[i] = result(t, *ccs::simulation_run(sims[i + 1]));
        }
        return t["aggregate"](t, r);
    }

    double result(std::span<const double> res)
    {
        auto& lua = lua_state();
        sol::table t = lua["Simulations"];
        return t["aggregate"](t, res);
    }

    double result(int i, std::span<const double> res)
    {
        auto& lua = lua_state();
        sol::table t = lua["Simulations"][i + 1];
        return t["aggregate"](t, res);
    }

    int simulation_size()
    {
        auto& lua = lua_state();
        // really need to abstract this
        sol::table t = lua["Simulations"];
        assert(t.valid());
        return t.size();
    }

    int simulation_size(int i)
    {
        auto& lua = lua_state();
        // really need to abstract this
        sol::table t = lua["Simulations"][i + 1]["simulations"];
        assert(t.valid());
        return t.size();
    }

    bool accept(double v)
    {
        auto& lua = lua_state();
        sol::table t = lua["Simulations"];
        return t["accept"](t, v);
    }

    double time_limit()
    {
        auto& lua = lua_state();
        sol::optional<double> t = lua["wallclock_hours"];
        // convert hours to seconds for comparison with Realm::Clock
        return t ? *t * 3600 : std::numeric_limits<double>::max();
    }

    operator Legion::TaskArgument() const { return {buf.data(), buf.size()}; }
};

} // namespace detail

using driver = detail::driver_<std::vector<std::byte>>;
using driver_span = detail::driver_<std::span<std::byte>>;
