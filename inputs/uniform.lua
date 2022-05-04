local floating_dims = 6
local dirichlet_dims = 3
local DIMS = floating_dims + dirichlet_dims
local uniform = 1 - 1e-8

-- parameters for driving nlopt

NLopt = {
   algorithm = "LN_COBYLA",
   dims = DIMS,
   xtol_rel = 1e-5,
   xtol_abs = 1e-8,
   maxeval = 100,
   initial_step = 0.1,
}



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
            psi = uniform,
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
         floating_alpha = {},
         dirichlet_alpha = {}
      },
      system = {
         type = "eigenvalues"
      }
   }
}

local wave_base_1d = {
   {
      logging = false,
      mesh = {
         index_extents = {51},
         domain_bounds = {2}
      },
      domain_boundaries = {
         xmin = "dirichlet"
      },
      scheme = {
         order = 1,
         type = "E2-poly",
         floating_alpha = {},
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
         max_time = 50,
         cfl = {
            hyperbolic = 0.8
         }
      }
   }
}


local wave_base_2d = {
   {
      logging = false,
      mesh = {
         index_extents = {51, 51},
         domain_bounds = {2, 2}
      },
      domain_boundaries = {
         xmin = "dirichlet",
         ymin = "dirichlet"
      },
      scheme = {
         order = 1,
         type = "E2-poly",
         floating_alpha = {},
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
         max_time = 50,
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


-- will be reusing this function in our tables
local function set_values (self, i, lst)
   assert(#lst >= DIMS)
   local s = self.simulations[i].scheme
   for j = 1, floating_dims do
      s.floating_alpha[j] = lst[j]
   end
   for j = 1, dirichlet_dims do
      s.dirichlet_alpha[j] = lst[j + floating_dims]
   end
end

-- only add a bunch of psi constraints on the left boundary since we're only
local eigenvalue_constraint = {
   simulations = add_entries(
      eigenvalue_base,
      {31, 41, 51, 61, 71, 81, 91, 101, 111},
      "mesh",
      "index_extents",
      1),
   set_values = set_values,
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

local wave_1d = {
   simulations = add_entries(
      add_entries(wave_base_1d, {51, 61, 71, 81, 91, 101}, "mesh", "index_extents", 1),
      {0.01, 0.1, 0.8},
      "step_controller",
      "cfl",
      "hyperbolic"
   ),
   --simulations = add_entries(wave_base, {51, 71, 91}, "mesh", "index_extents", 1),
   set_values = set_values,
   result = function(self, lst)
      -- lst = {time, error, ...}
      local max_time = wave_base_1d[1].step_controller.max_time
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
      return math.sqrt(result / #lst)
   end
}

local wave_2d = {
   simulations = add_entries(
      add_entries(wave_base_2d, { {91, 91}, {71, 71}, {51, 51}}, "mesh", "index_extents"),
      {0.8, 0.05},
      "step_controller",
      "cfl",
      "hyperbolic"
   ),
   set_values = set_values,
   result = function(self, lst)
      -- lst = {time, error, ...}
      local max_time = wave_base_2d[1].step_controller.max_time
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
      return math.sqrt(result / #lst)
   end
}

Simulations = {
   wave_1d,
   wave_2d,

   aggregate = function(self, lst)
      local result = 0.0
      for _, v in ipairs(lst) do
         result = result + v * v
      end
      return math.sqrt(result / #lst)
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
