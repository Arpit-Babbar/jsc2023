
using OrdinaryDiffEq
using Trixi
include("../common/m2000_callback.jl")
include("../common/lw_cfl.jl")

###############################################################################

equations = CompressibleEulerEquations2D(5.0/3.0)

@inline function initial_condition_astro_jet(x, t, equations::CompressibleEulerEquations2D)
  rho  = 0.5
  v1 = 0.0
  v2 = 0.0
  p = 0.4127
  prim = SVector(rho, v1, v2, p)
  return prim2cons(prim, equations)
end

initial_condition = initial_condition_astro_jet

@inline function boundary_condition_astro_jet(x, t, equations::CompressibleEulerEquations2D)
  if t > 0 && x[2] >= -0.05 && x[2] <= 0.05 && x[1] ≈ 0.0
    rho  = 5.0
    v1 = 800
    v2 = 0.0
    p  = 0.4127
  else
    rho  = 0.5
    v1 = 0.0
    v2 = 0.0
    p  = 0.4127
  end
  prim = SVector(rho, v1, v2, p)
  return prim2cons(prim, equations)
end

boundary_condition_inflow = BoundaryConditionDirichlet(boundary_condition_astro_jet)

@inline function boundary_condition_outflow(u_inner,
                                            normal_direction::AbstractVector,
                                            x, t,
                                            surface_flux_function, equations::CompressibleEulerEquations2D)
  # NOTE: Only for the supersonic outflow is this strategy valid
  # Calculate the boundary flux entirely from the internal solution state
  return flux(u_inner, normal_direction, equations)
end

boundary_conditions = Dict( :Bottom => boundary_condition_outflow,
                            :Top    => boundary_condition_outflow,
                            :Right  => boundary_condition_outflow,
                            :Left   => boundary_condition_inflow   )
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

mesh_file = joinpath(@__DIR__, "M2000.inp")
isfile(mesh_file) || download("https://gist.githubusercontent.com/Arpit-Babbar/88216acbd3d6a257b12f3d03e03c3584/raw/db2b7a23767276a5afc4ca145ecf25aae6ea8e19/M2000.inp", mesh_file)
#                              mesh_file)isfile(mesh_file) || download("https://gist.githubusercontent.com/Arpit-Babbar/627f19ef40127b84624429b2f0b9e7f0/raw/c0f4adc777bd5eda7c38002f57f09de9277659f3/M2000_bigger.inp",
#                              mesh_file)
mesh = P4estMesh{2}(mesh_file, initial_refinement_level = 3)

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 1.0
ode = semidiscretize(semi, (0.0, 0.001));

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     extra_analysis_integrals=(entropy,))

alive_callback = AliveCallback(analysis_interval=50)

# Use LW-D2 CFL with a safety factor of 0.5
stepsize_callback = StepsizeCallbackM2000(cfl = trixi2lw(0.5, solver))

save_solution = SaveSolutionCallback(interval=10000000000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim,
				     output_directory = joinpath(@__DIR__, "out"))

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution, stepsize_callback
                       )

# positivity limiter necessary for this example with strong shocks
stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(5.0e-6, 5.0e-6),
                                                     variables=(Trixi.density, pressure))

###############################################################################
# run the simulation
sol = solve(ode, 
            # SSPRK43(stage_limiter!),
            SSPRK54(stage_limiter!), dt = 1e-8,
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary

# 1.3429597560927415
# 22.23376230393017
