
using OrdinaryDiffEq
using H3PC


###############################################################################
# semidiscretization of the compressible Euler equations
mixture     = "air_5"
state_model = "ChemNonEq1T"
equations = HypersonicEulerEquations2D(mixture,state_model)

function initial_condition_blast_wave(x, t, equations::HypersonicEulerEquations2D)
    # Modified From Hennemann & Gassner JCP paper 2020 (Sec. 6.3) -> "medium blast wave"
    # Set up polar coordinates
    inicenter = SVector(0.0, 0.0)
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)
    
    if r <= 0.5
        rho1 = 0.027912209489708243
        rho2 =  0.008941625293754009
        rho3 = 3.4930594204054646e-5
        rho4 = 0.0015825170002801239
        rho5 = 5.017658307133176e-7
        rho_v1 = 0.0
        rho_v2 = 0.0
        rho_e  = 1.425938157389019e6
    else
        rho1 = 8.831713479861686e-81
        rho2 = 4.627828119933247e-42
        rho3 =  3.1197885320752415e-17
        rho4 = 0.08872316211371443
        rho5 = 0.026939968565103053
        rho_v1 = 0.0
        rho_v2 = 0.0
        rho_e  = -9783.724551855581
    end

    # U = [[0.0 0.0];
    #      [0.0 0.0];
    #      [1.425938157389019e6 -9783.724551855581];
    #      [0.027912209489708243 8.831713479861686e-81];
    #      [0.008941625293754009 4.627828119933247e-42];
    #      [3.4930594204054646e-5  3.1197885320752415e-17];
    #      [0.0015825170002801239 0.08872316211371443];
    #      [5.017658307133176e-7 0.026939968565103053]]
    # coef = 16.0;
    # rr = tanh(coef*(2.0*r-1.0))
    # uu = Vector{Float64}(undef, 8)
    # uu[:] = 0.5*(rr+1.0)*(U[:,2]-U[:,1])+U[:,1]

    return SVector(rho_v1, rho_v2, rho_e, rho1, rho2, rho3, rho4, rho5)
    # return SVector(uu[1], uu[2], uu[3], uu[4], uu[5], uu[6], uu[7], uu[8])
end
initial_condition = initial_condition_blast_wave
# for i=1:100
#     x= i*1.0/100.0
#     xx= initial_condition([x 0.0],0.0,equations)
#     println(xx[1],"\t",xx[2],"\t",xx[3],"\t",xx[4],"\t",xx[5],"\t",xx[6],"\t",xx[7],"\t",xx[8])
# end
boundary_conditions = Dict(:x_neg => BoundaryConditionDirichlet(initial_condition),
                           :y_neg => BoundaryConditionDirichlet(initial_condition),
                           :y_pos => BoundaryConditionDirichlet(initial_condition),
                           :x_pos => BoundaryConditionDirichlet(initial_condition))
# boundary_condition = BoundaryConditionDirichlet(initial_condition)

surface_flux = flux_lax_friedrichs
volume_flux = flux_ducros#flux_kennedy_gruber
basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

# coordinates_min = (-2.0, -2.0)
# coordinates_max = (2.0, 2.0)

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)

trees_per_dimension = (16, 16)
mesh = P4estMesh(trees_per_dimension,
                 polydeg = 3, initial_refinement_level = 1,
                 coordinates_min = coordinates_min, coordinates_max = coordinates_max,
                 periodicity = (false, false))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions,
                                    source_terms = source_terms_grossman)

###############################################################################
tspan = (0.0, 0.0002)
# tspan = (0.0, 0.0000)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt=0.0002/50.0,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

# amr_indicator = IndicatorHennemannGassner(semi,
#                                           alpha_max = 0.5,
#                                           alpha_min = 0.001,
#                                           alpha_smooth = true,
#                                           variable = Trixi.density)
# amr_controller = ControllerThreeLevel(semi, amr_indicator,
#                                       base_level = 0,
#                                       max_level = 3, max_threshold = 0.02)
# amr_callback = AMRCallback(semi, amr_controller,
#                            interval = 5,
#                            adapt_initial_condition = true,
#                            adapt_initial_condition_only_refine = true)

stepsize_callback = StepsizeCallback(cfl = 0.2)

callbacks = CallbackSet(summary_callback,stepsize_callback,
                        analysis_callback, alive_callback,save_solution)
tresh = 0.0
stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds = (tresh,tresh,tresh,tresh,tresh,0.0),
                                                     variables = (H3PC.rho1,H3PC.rho2,H3PC.rho3,
                                                        H3PC.rho4,H3PC.rho5,pressure))
###############################################################################
# run the simulation

# sol = solve(ode, CarpenterKennedy2N54(stage_limiter!,williamson_condition = false),
#             dt = 1e-6, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep = false, callback = callbacks);

sol = solve(ode, SSPRK43();
            ode_default_options()..., callback = callbacks);
# sol = solve(ode, RDPK3SpFSAL49(stage_limiter!); abstol = 1.0e-7, reltol = 1.0e-7,
#             ode_default_options()..., callback = callbacks);

# sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
#             dt = 5e-8, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep = false, callback = callbacks);

summary_callback() # print the timer summary
