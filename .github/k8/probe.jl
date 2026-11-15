# Diagnostic for K8, not for merging: at the points where the failing structure tests evaluate the
# 3D guiding-centre fields, compare every entry of `DDb♭`, from ElectromagneticFields and from this
# package's field point, with a central difference of `Db♭`, and save the generated code of `Db♭`
# and `DDb♭`.
# Usage: julia --project=. probe.jl <output directory>; exit code 2 when an entry is wrong.
using ChargedParticleDynamics
using ElectromagneticFields
using ElectromagneticFields: RuntimeGeneratedFunctions

const G = ChargedParticleDynamics.GuidingCenter3d
const OUT = ARGS[1]
mkpath(OUT)

function findrgf(x, depth = 0)
    depth > 6 && return nothing
    x isa RuntimeGeneratedFunctions.RuntimeGeneratedFunction && return x
    for n in fieldnames(typeof(x))
        isdefined(x, n) || continue
        r = findrgf(getfield(x, n), depth + 1)
        r === nothing || return r
    end
    nothing
end

# The code of a generated function, as `get_expression` builds it; RuntimeGeneratedFunctions 0.5.9,
# the floor, has no `get_expression`.
function rgfexpr(::RuntimeGeneratedFunctions.RuntimeGeneratedFunction{
        A, C, X, I}) where {A, C, X, I}
    Expr(:->, Expr(:tuple, A...), RuntimeGeneratedFunctions._lookup_body(C, I))
end

cpu = string(Sys.cpu_info()[1].model, ", ", Sys.CPU_NAME, ", ", VERSION)
write(joinpath(OUT, "cpu.txt"), cpu)
println("[k8] ", cpu)

# The points of `test/integration/structure.jl`: the barely passing initial condition (field
# accessors, constraint derivatives), that point displaced by `0.01 .* [1, 2, 3]` (multiplier
# derivatives), and the point on the loop of each `b = e₃` fixture (constraint derivatives).
points = []
for M in (G.SolovevIterXpoint, G.TokamakMediumCartesian, G.TokamakSmallToroidal)
    q = M.initial_conditions_barely_passing().q
    push!(points, (M, q), (M, q .+ 0.01 .* [1.0, 2.0, 3.0]))
end
for M in (G.SymmetricField, G.ThetaPinchField)
    push!(points, (M, M.initial_conditions(0.0, M.f_loop(0.1)).q))
end

for M in unique(first.(points))
    label = string(nameof(M))
    for name in (:Db♭, :DDb♭)
        r = findrgf(getfield(M.FIELD.functions, name))
        body = r === nothing ? "no RuntimeGeneratedFunction found" : string(rgfexpr(r))
        h = string(hash(body); base = 16)
        write(joinpath(OUT, "$(label)_$(name)_$(h).jl"), body)
        println("[k8] ", label, " ", name, " code hash ", h)
    end
end

total = 0
for (M, q) in points
    label = string(nameof(M))
    t = 0.0
    D = ElectromagneticFields.DDb♭(M.FIELD, t, q)
    S = G.fieldpoint²(M.FIELD, t, q)
    bad = 0
    for i in 1:3, j in 1:3, k in 1:3
        δ = 1e-6 * max(1.0, abs(q[k]))
        qp = copy(q)
        qp[k] += δ
        qm = copy(q)
        qm[k] -= δ
        fd = (ElectromagneticFields.Db♭(M.FIELD, t, qp)[i, j] -
              ElectromagneticFields.Db♭(M.FIELD, t, qm)[i, j]) / 2δ
        ok = isapprox(D[i, j, k], fd; rtol = 1e-5, atol = 1e-8) &&
             isapprox(S.values.DDb♭[i, j, k], fd; rtol = 1e-5, atol = 1e-8)
        if !ok
            bad += 1
            println("[k8] BAD ", label, " q=", q, " ", (i, j, k), " EMF=", D[i, j, k],
                " CPD=", S.values.DDb♭[i, j, k], " fd=", fd)
        end
    end
    println("[k8] ", label, " q=", q, " wrong entries: ", bad)
    global total += bad
end
println("[k8] total wrong entries: ", total)
exit(total == 0 ? 0 : 2)
