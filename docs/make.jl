push!(LOAD_PATH,"../src/")
using Documenter
using IonicElectrochemicalCells

makedocs(
    sitename = "IonicElectrochemicalCells",
    format = Documenter.HTML(),
    modules = [IonicElectrochemicalCells]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/fafroo/IonicElectrochemicalCells.jl>"
)
