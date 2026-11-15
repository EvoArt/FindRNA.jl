# Until ProfileHMMs is registered: use a sibling checkout if there is one
# (../ProfileHMMs), otherwise the GitHub repository.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
sibling = joinpath(@__DIR__, "..", "..", "ProfileHMMs")
if isfile(joinpath(sibling, "Project.toml"))
    Pkg.develop(path=sibling)
else
    Pkg.add(url="https://github.com/EvoArt/ProfileHMMs.jl")
end
Pkg.instantiate()
