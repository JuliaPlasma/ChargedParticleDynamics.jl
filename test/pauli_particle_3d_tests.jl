
using SafeTestsets

module PauliParticle3dTests

using GeometricIntegrators
using Test

export test_pauli_particle_3d

# ---------------------------------------------------------------------------------------------
# Time step and time span limits of the Pauli particle equilibria
# ---------------------------------------------------------------------------------------------
#
# Every call in this file runs at its module's declared thousand-step example; there are no
# exceptions. The steps at which these models stop integrating, worth knowing before changing one:
#
#   SolovevIter            Δt = 1    all three formulations throw `NaN detected in direction vector`
#   SolovevIterXpoint      Δt = 1    `pode`/`hode` throw on a NaN direction, `iode` reaches 2.7
#   TokamakSmallCartesian  Δt = 500  `pode`/`hode` reach 1.9E+37 with 441 capped steps, `iode` throws
#
# `TokamakSmallCylindrical` and `TokamakSmallToroidal` describe the same tokamak and keep Δt = 500,
# at which they are well behaved (2.5E-6 and 1.2E-6). The asymmetry is real: it is the cartesian
# chart of this equilibrium that is stiff, not the equilibrium.
#
# `SolovevIter` and `SolovevIterXpoint` are short for a reason that is not one formulation being
# catered to: over t ∈ [0, 1E3] `pode`, `hode` and `iode` alike lose `barely_passing`,
# `barely_trapped` and `deeply_passing` (2.7, 2.7, 5.6 in relative energy) while `deeply_trapped` and
# `trapped` survive. The orbit runs out, not a discretisation of it.
#
# The theta pinch keeps `MidpointExtrapolation(5)` for a different reason, recorded at that block:
# its momentum is constant along the orbit, so the default Hermite guess is degenerate.

# See `guiding_center_3d_tests.jl` for why `f_abstol` has to stay above the residual's round-off
# floor and why `f_reltol` is left at its default.
const options = (f_abstol=1E-12, max_iterations=50, warn_iterations=50)

# Asserts that the integration returns a solution rather than `@test_nowarn`; see the comment on
# the corresponding functions in `guiding_center_3d_tests.jl` for why.
function test_pauli_particle_3d(equ::Union{HODEProblem,PODEProblem}; kwargs...)
    @test integrate(equ, PartitionedGauss(2); options..., kwargs...) isa GeometricSolution
end

function test_pauli_particle_3d(equ::Union{IODEProblem,LODEProblem}; kwargs...)
    @test integrate(equ, VPRKGauss(2); options..., kwargs...) isa GeometricSolution
end

end


@safetestset "Pauli Particle in 3D with ITER-like Solov'ev Equilibrium                                            " begin

    using ChargedParticleDynamics.PauliParticle3d.SolovevIter
    using ..PauliParticle3dTests

    # 9.7E-3 in relative energy over the run. `Δt = 1` is not usable here: all three formulations
    # throw `NaN detected in direction vector`.
    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_barely_passing()))
    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_barely_trapped()))
    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_deeply_passing()))
    test_pauli_particle_3d(SolovevIter.hodeproblem(initial_conditions_deeply_trapped()))

    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_barely_passing()))
    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_barely_trapped()))
    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_deeply_passing()))
    test_pauli_particle_3d(SolovevIter.iodeproblem(initial_conditions_deeply_trapped()))

end


@safetestset "Pauli Particle in 3D with ITER-like Solov'ev Equilibrium with X-Point                               " begin

    using ChargedParticleDynamics.PauliParticle3d.SolovevIterXpoint
    using ..PauliParticle3dTests

    # As for `SolovevIter`: 1.1E-2 here, where `Δt = 1` throws on a `NaN` direction in the `pode` and
    # `hode` forms and reaches 2.7 in the `iode` one.
    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_barely_passing()))
    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_barely_trapped()))
    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_deeply_passing()))
    test_pauli_particle_3d(SolovevIterXpoint.hodeproblem(initial_conditions_deeply_trapped()))

    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_barely_passing()))
    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_barely_trapped()))
    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_deeply_passing()))
    test_pauli_particle_3d(SolovevIterXpoint.iodeproblem(initial_conditions_deeply_trapped()))

end


@safetestset "Pauli Particle in 3D in Symmetric Magnetic Field                                                    " begin

    using ChargedParticleDynamics.PauliParticle3d.SymmetricField
    using ..PauliParticle3dTests

    test_pauli_particle_3d(SymmetricField.podeproblem())
    test_pauli_particle_3d(SymmetricField.hodeproblem())
    test_pauli_particle_3d(SymmetricField.iodeproblem())

end


@safetestset "Pauli Particle in 3D in Theta Pinch                                                                 " begin

    using GeometricIntegrators: MidpointExtrapolation
    using ChargedParticleDynamics.PauliParticle3d.ThetaPinchField
    using ..PauliParticle3dTests

    # the momentum is constant along the initial trajectory, so the default Hermite
    # extrapolation of the initial guess is degenerate
    test_pauli_particle_3d(ThetaPinchField.podeproblem(); initialguess=MidpointExtrapolation(5))
    test_pauli_particle_3d(ThetaPinchField.hodeproblem(); initialguess=MidpointExtrapolation(5))
    test_pauli_particle_3d(ThetaPinchField.iodeproblem(); initialguess=MidpointExtrapolation(5))

end


@safetestset "Pauli Particle in 3D in ITER Equilibrium in Cylindrical Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakIterCylindrical.podeproblem())
    test_pauli_particle_3d(TokamakIterCylindrical.hodeproblem())
    test_pauli_particle_3d(TokamakIterCylindrical.iodeproblem())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cartesian Coordinates                                " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian
    using ..PauliParticle3dTests

    # All three formulations agree at 4.5E-8 here. `Δt = 500` — which the `Cylindrical` and
    # `Toroidal` charts of this same tokamak use happily — does not work in this one: `pode` and
    # `hode` reach 1.9E+37 with 441 capped steps and `iode` throws. It is the cartesian chart that is
    # stiff, not the equilibrium.
    test_pauli_particle_3d(TokamakSmallCartesian.podeproblem())
    test_pauli_particle_3d(TokamakSmallCartesian.hodeproblem())
    test_pauli_particle_3d(TokamakSmallCartesian.iodeproblem())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Cylindrical Coordinates                              " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallCylindrical.podeproblem())
    test_pauli_particle_3d(TokamakSmallCylindrical.hodeproblem())
    test_pauli_particle_3d(TokamakSmallCylindrical.iodeproblem())

end

@safetestset "Pauli Particle in 3D in Tokamak Equilibrium in Toroidal Coordinates                                 " begin

    using ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal
    using ..PauliParticle3dTests

    test_pauli_particle_3d(TokamakSmallToroidal.podeproblem())
    test_pauli_particle_3d(TokamakSmallToroidal.hodeproblem())
    test_pauli_particle_3d(TokamakSmallToroidal.iodeproblem())

end
