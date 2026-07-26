using Documenter
using ChargedParticleDynamics

using Base.Docs: Binding, meta

const CP = ChargedParticleDynamics

# Every module in the package, collected recursively: the equilibria live one level down, in a module
# per equilibrium, and there are more than fifty of them, so a hand-written list is exactly what
# fails silently when the fifty-first is added.
function submodules(m::Module, seen = Set{Module}())
    m in seen && return Module[]
    push!(seen, m)
    mods = Module[m]
    for name in names(m; all = true, imported = false)
        isdefined(m, name) || continue
        sub = getfield(m, name)
        sub isa Module && parentmodule(sub) === m && sub !== m || continue
        append!(mods, submodules(sub, seen))
    end
    return mods
end

# Whether `m` carries a docstring of its own. A module's docstring is stored in the module's *own*
# metadata, keyed by the binding in its parent, which is why this is not just `Docs.doc(m)` — that
# has a fallback and never reports failure.
function hasdocstring(m::Module)
    md = meta(m; autoinit = false)
    md === nothing && return false
    return haskey(md, Binding(parentmodule(m), nameof(m)))
end

# Every fully qualified name appearing inside a `@docs` or `@autodocs` block, as a set of exact
# strings. Exact set membership rather than a substring search, so that `SolovevIter` is not
# considered documented merely because `SolovevIterXpoint` is listed.
function documented_names(dir)
    found = Set{String}()
    for (root, _, files) in walkdir(dir), file in files
        endswith(file, ".md") || continue
        text = read(joinpath(root, file), String)
        for block in eachmatch(r"```@(?:auto)?docs\r?\n(.*?)```"s, text)
            for name in eachmatch(r"ChargedParticleDynamics(?:\.\w+)*", block[1])
                push!(found, name.match)
            end
        end
    end
    return found
end

# Documenter's own missing-docs check is deliberately left off, and this check replaces it.
#
# The reason it cannot be used is structural rather than a matter of configuration. Enabling it means
# passing `modules`, and `Documenter.DocumentBlueprint` stores
# `submodules(modules; ignore = checkdocs_ignored_modules)` in a single field that serves two
# purposes at once: `docchecks.jl` uses it for the list of bindings that *must* be documented, and
# `expander_pipeline.jl` uses it to filter which docstrings a `@docs` or `@autodocs` block is
# *allowed* to render. One set, both jobs.
#
# That is irreconcilable with how this package is laid out. Each equilibrium module `include`s the
# shared model code, so all fifty-odd of them carry their own copy of every shared docstring — 491
# across the package — and the model pages render each once, from an anchor module, listing the rest
# in `@docs` blocks for their module docstring alone. Passing every module makes the pages work and
# reports 382 docstrings as "not included in the manual", every one of them excluded on purpose.
# Narrowing the set with `checkdocs_ignored_modules` silences those, and simultaneously makes the
# `@docs` blocks naming those modules fail with "no docs found", because it is the same set.
#
# So the check below does the job directly, and more precisely than the built-in would: every module
# in the package must have a docstring of its own, and must be named in some `@docs`/`@autodocs`
# block. Add an equilibrium and forget the page, or add one without a docstring, and the build says
# so. Warnings rather than errors, so the build still completes.
#
# What this does not cover is a *function* docstring that stops being rendered. The pages render the
# shared API through `@autodocs` on the anchor modules, which takes everything those modules own, so
# that can only regress by someone replacing an `@autodocs` block with a partial `@docs` list.
let documented = documented_names(joinpath(@__DIR__, "src")), undocumented = String[], unlisted = String[]
    for m in submodules(CP)
        m === CP && continue
        hasdocstring(m) || push!(undocumented, string(m))
        string(m) in documented || push!(unlisted, string(m))
    end
    isempty(undocumented) ||
        @warn "Modules without a docstring of their own:\n" * join("  " .* sort(undocumented), "\n")
    isempty(unlisted) ||
        @warn "Modules absent from every @docs/@autodocs block:\n" * join("  " .* sort(unlisted), "\n")
    isempty(undocumented) && isempty(unlisted) &&
        @info "Module docs check: all $(length(submodules(CP)) - 1) submodules have a docstring and appear in the manual."
end

makedocs(
    sitename = "ChargedParticleDynamics.jl",
    # No `modules` argument: see the comment on the module docs check above for why it cannot be
    # combined with rendering each shared docstring once.
    warnonly = Documenter.except(:autodocs_block, :cross_references, :docs_block, :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes, :linkcheck, :meta_block, :parse_error, :setup_block),
    format = Documenter.HTML(
                prettyurls = get(ENV, "CI", nothing) == "true",
                assets = [asset("assets/style.css", class=:css, islocal=true)],
                example_size_threshold = 32,
                ),
    pages = ["Overview"                      => "index.md",
             "Normalization"                 => "normalization.md",
             "Initialization"                => "initialization.md",
             "Charged Particles in 3D"       => "charged_particle_3d.md",
             "Pauli Particles in 3D"         => "pauli_particle_3d.md",
             "Guiding Center Dynamics in 3D" => "guiding_center_3d.md",
             "Guiding Center Dynamics in 4D" => "guiding_center_4d.md",
             "Gyrokinetics in 4D"            => "gyro_kinetics_4d.md",
             "Model Audit"                   => "audit.md",
             "Findings"                      => "findings.md",
             "Examples" => [
                    "ITER Equilibrium in Cylindrical Coordinates"   => "examples/iter_cylindrical.md",
                 ],
            ]
)

deploydocs(
    repo   = "github.com/JuliaPlasma/ChargedParticleDynamics.jl",
    devurl = "latest",
    devbranch = "main",
    push_preview = true,
)
