using SafeTestsets

# Loading `CairoMakie` and `LaTeXStrings` activates the `ChargedParticlePlots` extension. Until
# these tests existed the extension was never loaded by the suite, which is how two breaks from the
# GeometricSolutions 0.6 restructure went unnoticed: `sol.t` became a `ScalarDataSeries` rather than
# a `TimeSeries`, so the `compute_*` convenience methods stopped dispatching, and
# `compute_relative_error` no longer accepts the plain `Vector` that `PoincareInvariants.compute!`
# returns.
#
# These are smoke tests: each function must build and return a Makie `Figure` without error.

@safetestset "Plotting extension and diagnostics                                                                  " begin

    using CairoMakie
    using LaTeXStrings
    using GeometricIntegrators
    using PoincareInvariants
    using ChargedParticleDynamics
    using ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical

    # The equilibrium is the module itself: the magnetic field code is injected into it by
    # `ElectromagneticFields.@code`, so `equ.R`, `equ.X`, … are its functions.
    equ = ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical

    prob = odeproblem(initial_conditions_deeply_passing())
    sol = integrate(prob, Gauss(2))

    @testset "Diagnostics" begin
        @test compute_energy(sol) isa ScalarDataSeries
        @test compute_energy_error(sol) isa Tuple
        @test compute_toroidal_momentum(sol) isa ScalarDataSeries
        @test compute_toroidal_momentum_error(sol) isa Tuple
        @test cartesian_solution(sol, equ) isa NamedTuple
    end

    @testset "Trajectory and field plots" begin
        cs = cartesian_solution(sol, equ)
        R = [q[1] for q in sol.q]
        Z = [q[2] for q in sol.q]
        φ = [q[3] for q in sol.q]
        u = [q[4] for q in sol.q]

        # `plot_invariant` has to be qualified: `PoincareInvariants` 0.5 exports a `plot_invariant`
        # of its own (for the error of a Poincaré invariant over an `EnsembleSolution`), and the
        # guiding centre Poincaré workflow requires both packages to be loaded together, so the two
        # names collide for every user of it.
        @test ChargedParticleDynamics.plot_invariant(sol.t, compute_energy(sol); label = "H") isa Figure
        @test plot_invariant_error(sol.t, compute_energy(sol); label = "H") isa Figure
        @test plot_trajectory_poloidal(R, Z) isa Tuple
        @test plot_trajectory_projection(cs.X, cs.Y) isa Tuple
        @test plot_trajectory_3d(cs.X, cs.Y, cs.Z) isa Tuple
        @test plot_trajectory_cylindrical(sol.t, R, Z, φ, u) isa Tuple
        @test plot_fieldlines(equ; xrange = (1.0, 2.5), yrange = (-0.8, 0.8), ngrid = (40, 30)) isa Tuple
    end

    # `plot_fieldlines` contours ψ = R A_φ, which is only meaningful in an axisymmetric cylindrical
    # chart. It used to accept any equilibrium and draw a plausible-looking wrong picture for the
    # cartesian and toroidal ones; it now refuses, and so does `plot_trajectory_poloidal(R, Z, equ)`,
    # which forwards its equilibrium here.
    @testset "Field lines reject charts they are not valid in" begin
        cylindrical = ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical
        cartesian   = ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian
        toroidal    = ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal

        @test is_axisymmetric_cylindrical(cylindrical)
        @test !is_axisymmetric_cylindrical(cartesian)
        @test !is_axisymmetric_cylindrical(toroidal)

        for bad in (cartesian, toroidal)
            @test_throws ArgumentError plot_fieldlines(bad; xrange = (1.0, 2.5), yrange = (-0.8, 0.8), ngrid = (4, 4))
            @test_throws ArgumentError plot_trajectory_poloidal([1.0, 1.1], [0.0, 0.1], bad)
        end

        # Without an equilibrium the trajectory still plots, which is the documented way out.
        @test plot_trajectory_poloidal([1.0, 1.1], [0.0, 0.1]) isa Tuple
    end

    @testset "Poincaré integral invariants" begin
        # The first invariant on the PoincareInvariants 0.5 interface: a noncanonical `FirstPI`
        # built from the guiding centre one-form, advected as an ensemble of loop points.
        pinv = poincare_invariant_1st(100)
        lprob = loop_iodeproblem(; timespan = (0.0, 10.0), timestep = 2.5)
        lsol = integrate(loop_ensemble(lprob, pinv), VPRKGauss(2))
        I₁ = compute!(pinv, lsol, parameters(lprob))

        @test I₁ isa Vector
        @test length(I₁) == 5
        # The dynamics is symplectic, so the invariant is conserved up to the integrator's error.
        @test maximum(abs, (I₁ .- I₁[begin]) ./ I₁[begin]) < 1E-2

        # The second invariant, on a finite-difference grid so that the surface can be plotted.
        pinv2 = poincare_invariant_2nd((8, 8); plan = PoincareInvariants.SecondFinDiffPlan)
        sprob = surface_iodeproblem(; timespan = (0.0, 10.0), timestep = 2.5)
        ssol = integrate(surface_ensemble(sprob, pinv2), VPRKGauss(2))
        I₂ = compute!(pinv2, ssol, parameters(sprob))

        @test I₂ isa Vector
        @test maximum(abs, (I₂ .- I₂[begin]) ./ I₂[begin]) < 1E-2

        # `plot_poincare_*` take plain coordinate vectors, one entry per plotted time.
        ts = [lsol[1].t[n] for n in 0:ntime(lsol[1])]
        cart(j, n) = to_cartesian(lsol[j].t[n], lsol[j].q[n])
        X = [[cart(j, n)[1] for j in 1:nsamples(lsol)] for n in 0:(length(ts)-1)]
        Y = [[cart(j, n)[2] for j in 1:nsamples(lsol)] for n in 0:(length(ts)-1)]
        Z = [[cart(j, n)[3] for j in 1:nsamples(lsol)] for n in 0:(length(ts)-1)]

        @test plot_poincare_invariant_error(ts, I₁) isa Figure
        @test plot_poincare_loop(ts, X, Y, Z; nplot = 3) isa Tuple
        @test plot_poincare_surface(ts, X, Y, Z; nplot = 3) isa Tuple
        @test plot_poincare_trajectories(X, Y, Z; nplot = 3) isa Tuple
    end

end
