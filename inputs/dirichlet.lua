local floating_dims = 0
local dirichlet_dims = 3
local DIMS = floating_dims + dirichlet_dims
local uniform = 1 - 1e-8
local max_time_default = 500

-- parameters for driving nlopt
NLopt = {
   algorithm = "LN_COBYLA",
   dims = 3,
   xtol_rel = 1e-5,
   xtol_abs = 1e-8,
   maxeval = 2,
   initial_step = 0.1,
   ftol_abs = 1e-2,
   ftol_rel = 1e-4
}


-- use the successful floating alphas from the uniform test
local floating_success = {
   0.7039278390946743,
   0.5390086175376538,
      -0.647109821986589,
   0.2051508287133347,
   0.6062051039572746,
   0.8148425279273044
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
         floating_alpha = floating_success,
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
      shapes = {
         {
            type = "yz_rect",
            psi = uniform,
            normal = 1,
            boundary_condition = "dirichlet"
         }
      },
      scheme = {
         order = 1,
         type = "E2-poly",
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
         max_time = max_time_default,
         --max_step = 2,
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
      shapes = {
         {
            type = "sphere",
            center = {1.053, 0.901},
            radius = math.pi / 10,
            boundary_condition = "dirichlet"
         }
      },
      scheme = {
         order = 1,
         type = "E2-poly",
         floating_alpha = floating_success,
         dirichlet_alpha = {}
      },
      system = {
         type = "scalar wave",
         max_error = 2.0
      },
      integrator = {
         type = "rk4"
      },
      step_controller = {
         max_time = max_time_default,
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

--
-- will be reusing these functions in our tables
--
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

local function l2 (self, lst)
   local result = 0.0
   for _, v in ipairs(lst) do
      result = result + v * v
   end
   return math.sqrt(result / #lst)
end

local function linf (self, lst)
   local result = -math.huge
   for _, v in ipairs(lst) do
      result = math.max(result, v)
   end
   return result
end



local psi_values = {1e-12, 1e-6, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, uniform}

local function add_psi (t)
   return add_entries(t, psi_values, "shapes", 1, "psi")
end

local function add_centers (t, lst)
   return add_entries(t, lst, "shapes", 1, "center")
end

local function add_mesh (t, lst)
   return add_entries(t, lst, "mesh", "index_extents")
end

local function add_cfl (t, lst)
   return add_entries(t, lst, "step_controller", "cfl", "hyperbolic")
end

function listify (lst)
   local v = {}
   for j = 1, #lst do
      v[j] = {lst[j]}
   end
   return v
end


-- only add a bunch of psi constraints on the left boundary since we're only
local eigenvalue_constraint = {
   simulations = add_psi(
      add_mesh(eigenvalue_base, listify({31, 41, 51, 61, 71, 81, 91, 101, 111}))
   ),

   set_values = set_values,

   result = function(self, lst)
      return lst[2]
   end,

   aggregate = linf
}

--
-- Set up 1d wave equation constraints with varying psi on left
--
local wave_1d = {
   simulations = add_cfl(
      add_psi(
         add_mesh(wave_base_1d, listify({51, 71}))),
      {0.1, 0.8}
   ),

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

   aggregate = l2
}

--
-- Set up 2d wave equation constraints with varying center positions
--

local wave_2d = {
   simulations = add_cfl(
      add_centers(
         add_mesh(wave_base_2d, { {51, 51}, {61, 61} }),
         {
            {1 + 1/14, 1 - 1/11},
            {1 - 1/17, 1 + 11/173},
            {1 - 1/19, 1 - 1/63},
            {1 + 7/111, 1 + 1/33}
         }
      ),
      {0.1,  0.8}
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
   aggregate = l2
}


Simulations = {
   wave_1d,
   wave_2d,

   aggregate = l2,

   accept = function(self, v)
      return v > 21
   end
}

Constraints = {
   eigenvalue_constraint
}

-- extra parameters
wallclock_hours = 5 / 60
