using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.develop(path=joinpath(@__DIR__, "..", "packages", "ProfileHMMs"))
Pkg.instantiate()
