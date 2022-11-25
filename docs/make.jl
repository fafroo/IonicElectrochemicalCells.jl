push!(LOAD_PATH,"../src/")
using Documenter, Literate
using IonicElectrochemicalCells


makedocs(
    sitename = "IonicElectrochemicalCells",
    format = Documenter.HTML(),
    modules = [IonicElectrochemicalCells],
    pages = Any[
        "Physics"=>"lattice_gas.md",
        "Types, Constructors and Methods"=>"allindex.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/fafroo/IonicElectrochemicalCells.jl>"
)
