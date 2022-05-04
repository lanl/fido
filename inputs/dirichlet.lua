-- parameters for driving nlopt
NLopt = {
   algorithm = "LN_COBYLA",
   dims = 3,
   xtol_rel = 1e-5,
   xtol_abs = 1e-8,
   maxeval = 50,
   initial_step = 0.1,
}

local DIMS = NLopt.dims

-- use the successful floating alphas from the uniform test
local floating_success = { -1.0234993183468384, 0.6236227747564413, -0.6492048678049246,
                           1.2671119217697635, 0.8246634386357885, 0.7017989818620513 }

local uniform = 1 - 1e-6 -- psi for uniform walls


-- base tables for eigenvalue constraints and simulation base objective
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
            psi = uniform,
            normal = -1,
            boundary_condition = "floating"
         }
      },
      scheme = {
         order = 1,
         type = "E2-poly",
         floating_alpha = floating_success,
         dirichlet_alpha = {3 / 25, 13 / 100, 7 / 50}
      },
      system = {
         type = "eigenvalues"
      }
   }
}

local wave_base = {
   {
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
         --alpha = {}
         floating_alpha = floating_success,
         dirichlet_alpha = {}
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
            hyperbolic = 0.8
         }
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
local eigenvalue_constraint = {
   simulations = add_entries(
      add_entries(eigenvalue_base, {31, 41, 51, 61, 71, 81, 91, 101, 111}, "mesh", "index_extents", 1),
      psi_values,
      "shapes",
      1,
      "psi"
   ),
   set_values = function(self, i, lst)
      assert(#lst >= DIMS)
      local s = self.simulations[i].scheme
      for j = 1, 3 do
         s.dirichlet_alpha[j] = lst[j]
      end
   end,
   result = function(self, lst)
      return lst[2]
   end,
   aggregate = function(self, lst)
      local result = -math.huge
      for _, v in ipairs(lst) do
         result = math.max(result, v)
      end
      return result
   end
}

local wave_simulation = {
   simulations = add_entries(
      add_entries(add_entries(wave_base, {51, 71}, "mesh", "index_extents", 1), psi_values, "shapes", 1, "psi"),
      {1.0, 0.5},
      "step_controller",
      "cfl",
      "hyperbolic"
   ),
   --simulations = add_entries(wave_base, {51, 71, 91}, "mesh", "index_extents", 1),
   set_values = function(self, i, lst)
      assert(#lst >= DIMS)
      local s = self.simulations[i].scheme
      for j = 1, 3 do
         s.dirichlet_alpha[j] = lst[j]
      end
   end,
   result = function(self, lst)
      -- lst = {time, error, ...}
      local max_time = wave_base[1].step_controller.max_time
      if lst[1] < max_time then
         return 20.0 * lst[1] / max_time
      else
         return 20.0 - math.log(lst[2])
      end
   end,
   aggregate = function(self, lst)
      local result = 0.0
      for _, v in ipairs(lst) do
         result = result + v * v
      end
      return math.sqrt(result) / #lst
   end
}
Simulations = {
   wave_simulation,

   aggregate = function(self, lst)
      return lst[1]
   end,

   accept = function(self, v)
      return v > 21;
   end
}

Constraints = {
   eigenvalue_constraint
}

-- extra parameters
wallclock_hours = 5 / 60
