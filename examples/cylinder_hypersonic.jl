using OrdinaryDiffEq
using H3PC

###############################################################################
# semidiscretization of the compressible Euler equations

mixture     = "air_5"
state_model = "ChemNonEq1T"
equations = HypersonicEulerEquations2D(mixture,state_model)

@inline function initial_condition_mach3_flow(x, t, equations::HypersonicEulerEquations2D)
  
  rho1 = 4.347779766802688e-114
  rho2 =  5.273650838884071e-60
  rho3 = 2.9208089786143016e-25
  rho4 = 0.0015474970136112847
  rho5 = 0.0004698831726471622
  rho = 0.002017380186258447
  v1 = 12*294.5166198856192
  rho_v1 = rho*v1
  rho_v2 = 0.0
  rho_e  = -145904.32485428682*rho+0.5*rho*v1*v1

    return SVector(rho_v1, rho_v2, rho_e, rho1, rho2, rho3, rho4, rho5)
end

initial_condition = initial_condition_mach3_flow


boundary_condition_inflow = BoundaryConditionDirichlet(initial_condition_mach3_flow)

# Outflow boundary condition.
# FIXME: For now only works for supersonic outflow where all values come from the internal state.
# The bones are here for the subsonic outflow as well. One simply needs to pass the reference pressure
# to set it from the outside, the rest comes from the internal solution.
# Once fixed this function should probably move to `compressible_euler_2d.jl`
# See the reference below for a discussion on inflow/outflow boundary conditions.
#
# - Jan-Reneé Carlson (2011)
#   Inflow/Outflow Boundary Conditions with Application to FUN3D.
#   [NASA TM 20110022658](https://ntrs.nasa.gov/citations/20110022658)
@inline function boundary_condition_outflow(u_inner, normal_direction::AbstractVector, x, t,
                                            surface_flux_function, equations::HypersonicEulerEquations2D)
  # # This would be for the general case where we need to check the magnitude of the local Mach number
  # norm_ = norm(normal_direction)
  # # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
  # normal = normal_direction / norm_

  # # Rotate the internal solution state
  # u_local = Trixi.rotate_to_x(u_inner, normal, equations)

  # # Compute the primitive variables
  # rho_local, v_normal, v_tangent, p_local = cons2prim(u_local, equations)

  # # Compute local Mach number
  # a_local = sqrt( equations.gamma * p_local / rho_local )
  # Mach_local = abs( v_normal / a_local )
  # if Mach_local <= 1.0
  #   p_local = # Set to the external reference pressure value (somehow? maybe stored in `equations`)
  # end

  # # Create the `u_surface` solution state where the local pressure is possibly set from an external value
  # prim = SVector(rho_local, v_normal, v_tangent, p_local)
  # u_boundary = prim2cons(prim, equations)
  # u_surface = Trixi.rotate_from_x(u_boundary, normal, equations)

  # Compute the flux using the appropriate mixture of internal / external solution states
  # flux = Trixi.flux(u_surface, normal_direction, equations)

  # NOTE: Only for the supersonic outflow is this strategy valid
  # Calculate the boundary flux entirely from the internal solution state
  flux = H3PC.flux(u_inner, normal_direction, equations)

  return flux
end

boundary_conditions = Dict( :Upper   => boundary_condition_inflow,
                            :Lower   => boundary_condition_inflow,
                            :Inlet   => boundary_condition_inflow,
                            :Cyl     => boundary_condition_slip_wall_riemann,
                            :Outlet  => boundary_condition_outflow)

volume_flux = flux_ducros
surface_flux = flux_lax_friedrichs

@inline function vel_mag(u, equations::HypersonicEulerEquations2D)
    rho_v1, rho_v2 = u
    rho = density(u,equations)
    return sqrt((rho_v1/rho)^2+(rho_v2/rho)^2)
end

polydeg = 3
basis = LobattoLegendreBasis(polydeg)
shock_indicator = IndicatorHennemannGassner(equations, basis,
                                            alpha_max=0.50,
                                            alpha_min=0.001,
                                            alpha_smooth=true,
                                            variable=H3PC.density)
volume_integral = VolumeIntegralShockCapturingHG(shock_indicator;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)
solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux, volume_integral=volume_integral)

# Get the unstructured quad mesh from a file (downloads the file if not available locally)
# default_mesh_file = joinpath(@__DIR__, "abaqus_forward_step.inp")
# isfile(default_mesh_file) || download("https://gist.githubusercontent.com/andrewwinters5000/b346ee6aa5446687f128eab8b37d52a7/raw/cd1e1d43bebd8d2631a07caec45585ec8456ca4c/abaqus_forward_step.inp",
#                                       default_mesh_file)
mesh_file = "cylinder_n.inp"

boundary_symbols = [:Upper, :Lower, :Cyl, :Inlet, :Outlet]
mesh = P4estMesh{2}(mesh_file,boundary_symbols = boundary_symbols)

#source_terms = source_terms_grossman
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions=boundary_conditions,source_terms = source_terms_grossman)
# source_terms = source_terms_grossman

# ###############################################################################
# # ODE solvers, callbacks etc.

tspan = (0.0, 10e-4)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 50
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     extra_analysis_integrals=(entropy,))

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim,
                                     output_directory="out_c")
# save_restart = SaveRestartCallback(interval = 4000,
#                                    save_final_restart = true)
# amr_indicator = IndicatorLöhner(semi, variable=density_pressure)#Trixi.density

# amr_indicator = IndicatorHennemannGassner(semi,
#                                           alpha_max=1.0,
#                                           alpha_min=0.001,
#                                           alpha_smooth=true,
#                                           variable=density_pressure)

# amr_controller = ControllerThreeLevel(semi, amr_indicator,
#                                       base_level=0,
#                                       med_level=2, med_threshold=0.05,
#                                       max_level=5, max_threshold=0.1)

# amr_controller = ControllerThreeLevel(semi, amr_indicator,
#                                       base_level=0,
#                                       max_level=4, max_threshold=0.1)
amr_indicator = IndicatorLöhner(semi, variable = vel_mag)

# amr_controller = ControllerThreeLevel(semi, amr_indicator,
#                                       base_level = 0,
#                                       med_level = 1, med_threshold = 0.1,
#                                       max_level = 2, max_threshold = 0.2)
amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level = 0,
                                      max_level = 3, max_threshold = 0.02)

# base_level: Minimum allowed refinement level
# med_level: Initial level of refinement for values above the med_threshold
# max_level: Maximum level of refinement for the indicator values above max_threshold
amr_callback = AMRCallback(semi, amr_controller,
                           interval=5,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)
stepsize_callback = StepsizeCallback(cfl=9.0e-4)

callbacks = CallbackSet(summary_callback,save_solution,stepsize_callback,
                        analysis_callback, alive_callback)

# positivity limiter necessary for this example with strong shocks
tresh = 0.0
stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds = (tresh,tresh,tresh,tresh,tresh,tresh),
                                                     variables = (H3PC.rho1,H3PC.rho2,H3PC.rho3,
                                                        H3PC.rho4,H3PC.rho5,pressure))

# ###############################################################################
# # run the simulation
# sol = solve(ode, CarpenterKennedy2N54(stage_limiter!),dt=1e-11,
#             save_everystep=false, callback=callbacks,
#             maxiters=999999);
sol = solve(ode, SSPRK43(stage_limiter!),
            save_everystep=false, callback=callbacks,
            maxiters=999999);
# sol = solve(ode, SSPRK932(stage_limiter!),
#             save_everystep=false, callback=callbacks,
#             maxiters=999999);
#CarpenterKennedy2N54,SSPRK43
summary_callback() # print the timer summary

