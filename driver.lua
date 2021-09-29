local inspect = require "inspect"

local eigenvalue_base = {
    {
        logging = false,
        mesh = {
            index_extents = {21},
            domain_bounds = {1}
        },
        shapes = {
            {
                type = "yz_rect",
                psi = 0.001,
                normal = 1,
                boundary_condition = "dirichlet"
            },
            {
                type = "yz_rect",
                psi = 0.9,
                normal = -1,
                boundary_condition = "floating"
            }
        },
        scheme = {
            order = 1,
            type = "E2-poly",
            floating_alpha = {13 / 100, 7 / 50, 3 / 20, 4 / 25, 17 / 100, 9 / 50},
            dirichlet_alpha = {3 / 25, 13 / 100, 7 / 50}
        },
        system = {
            type = "eigenvalues"
        }
    }
}

local function simple_copy(obj)
    if type(obj) ~= "table" then
        return obj
    end
    local res = {}
    for k, v in pairs(obj) do
        res[simple_copy(k)] = simple_copy(v)
    end
    return res
end

local function nested_set(tbl, v, x, ...)
    if ... == nil then
        tbl[x] = v
    else
        return nested_set(tbl[x], v, ...)
    end
end

local function add_entries(tbl, vs, ...)
    local out = {}
    for _, t in ipairs(tbl) do
        for _, v in ipairs(vs) do
            local i = #out + 1
            out[i] = simple_copy(t)
            nested_set(out[i], v, ...)
        end
    end
    return out
end

local psi_values = {1e-12, 1e-6, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 - 1e-8}

-- only add a bunch of psi constraints on the left boundary since we're only
Constraints = {
    {
        simulations = add_entries(
            add_entries(eigenvalue_base, {21, 31, 41, 51, 61, 71, 81, 91, 101, 111}, "mesh", "index_extents", 1),
            psi_values,
            "shapes",
            1,
            "psi"
        ),
        set_values = function(i, lst)
            assert(#lst >= 9)
            local s = Constraints[1].simulations[i].scheme
            for j = 1, 6 do
                s.floating_alpha[j] = lst[j]
            end
            for j = 1, 3 do
                s.dirichlet_alpha[j] = lst[j + 6]
            end
        end,
        result = function(lst)
            return lst[1]
        end,
        aggregate = function(lst)
            local result = -math.huge
            for _, v in ipairs(lst) do
                result = math.max(result, v)
            end
            return result
        end,
        debug = function(i)
            print("tbl: ", i, inspect(Constraints[1].simulations[i]))
        end
    }
}

Simulation = {
    logging = false,
    mesh = {
        index_extents = {51},
        domain_bounds = {2}
    },
    domain_boundaries = {
        xmin = "dirichlet"
    },
    shapes = {
        {
            type = "yz_rect",
            psi = 0.99,
            normal = 1,
            boundary_condition = "dirichlet"
        }
    },
    scheme = {
        order = 1,
        type = "E2-poly",
        alpha = {}
    },
    system = {
        type = "scalar wave",
        center = {-1},
        radius = 0,
        max_error = 2.0
    },
    integrator = {
        type = "rk4"
    },
    step_controller = {
        max_time = 100,
        --max_step = 2,
        cfl = {
            hyperbolic = 0.8,
            parabolic = 0.2
        }
    }
}
