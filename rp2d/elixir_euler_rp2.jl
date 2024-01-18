using Downloads: download
using OrdinaryDiffEq
using Trixi
include("../common/lw_cfl.jl")

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

"""
    initial_condition_double_mach_reflection(x, t, equations::CompressibleEulerEquations2D)

Compressible Euler setup for a double Mach reflection problem.
Involves strong shock interactions as well as steady / unsteady flow structures.
Also exercises special boundary conditions along the bottom of the domain that is a mixture of
Dirichlet and slip wall.
See Section IV c on the paper below for details.

- Paul Woodward and Phillip Colella (1984)
  The Numerical Simulation of Two-Dimensional Fluid Flows with Strong Shocks.
  [DOI: 10.1016/0021-9991(84)90142-6](https://doi.org/10.1016/0021-9991(84)90142-6)
"""
@inline function initial_condition_rp(x_, t, equations::CompressibleEulerEquations2D)
  x, y = x_[1], x_[2]

  if x >= 0.5 && y >= 0.5
     ρ  = 0.5313
     v1 = 0.0
     v2 = 0.0
     p  = 0.4
  elseif x < 0.5 && y >= 0.5
     ρ  = 1.0
     v1 = 0.7276
     v2 = 0.0
     p  = 1.0
  elseif x < 0.5 && y < 0.5
     ρ  = 0.8
     v1 = 0.0
     v2 = 0.0
     p  = 1.0
  elseif x >= 0.5 && y < 0.5
     ρ  = 1.0
     v1 = 0.0
     v2 = 0.7276
     p  = 1.0
  end
  rho = ρ

  prim = SVector(rho, v1, v2, p)
  return prim2cons(prim, equations)
end

initial_condition = initial_condition_rp


boundary_condition_inflow = BoundaryConditionDirichlet(initial_condition)

# Supersonic outflow boundary condition. Solution is taken entirely from the internal state.
# See `examples/p4est_2d_dgsem/elixir_euler_forward_step_amr.jl` for complete documentation.
@inline function boundary_condition_outflow(u_inner, normal_direction::AbstractVector, x, t,
                                            surface_flux_function, equations::CompressibleEulerEquations2D)
  # NOTE: Only for the supersonic outflow is this strategy valid
  # Calculate the boundary flux entirely from the internal solution state
  return flux(u_inner, normal_direction, equations)
end

# BCs must be passed as Dict
boundary_conditions = Dict(
  :x_neg => boundary_condition_outflow,
  :x_pos => boundary_condition_outflow,
  :y_neg => boundary_condition_outflow,
  :y_pos => boundary_condition_outflow
)

# Create DG solver with polynomial degree = 4 and (local) Lax-Friedrichs/Rusanov flux as surface flux
solver = DGSEM(polydeg=4, surface_flux=flux_lax_friedrichs)

# The initial condition is 2-periodic
coordinates_min = (-1.5, -1.5) # minimum coordinates (min(x), min(y))
coordinates_max = ( 2.5,  2.5) # maximum coordinates (max(x), max(y))

trees_per_dimension = (1024, 1024)

mesh = P4estMesh(trees_per_dimension, polydeg=4,
                 coordinates_min=coordinates_min, coordinates_max=coordinates_max,
                 periodicity=true)


volume_flux = flux_central
surface_flux = flux_lax_friedrichs

polydeg = 4
basis = LobattoLegendreBasis(polydeg)
shock_indicator = IndicatorHennemannGassner(equations, basis,
                                            alpha_max=0.5,
                                            alpha_min=0.001,
                                            alpha_smooth=true,
                                            variable=density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(shock_indicator;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)
solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux, volume_integral=volume_integral)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    # boundary_conditions=boundary_conditions
				   )

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.25)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     extra_analysis_integrals=(entropy,))

alive_callback = AliveCallback(analysis_interval=100)

save_solution = SaveSolutionCallback(interval=1000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)
stepsize_callback = StepsizeCallback(cfl=trixi2lw(0.98, solver))

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback
                       )

# positivity limiter necessary for this example with strong shocks
stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(5.0e-6, 5.0e-6),
                                                     variables=(Trixi.density, pressure))

###############################################################################
# run the simulation
sol = solve(ode, SSPRK54(stage_limiter!),
            dt = 1,
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary
