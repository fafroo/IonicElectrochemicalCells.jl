module IonicElectrochemicalCells

using Base: @kwdef
using CompositeStructs
using DataFrames            # TODO bump to examples
using DocStringExtensions   
using ExtendableGrids
using Markdown
using ModelParameters
using SamplerIO
using VoronoiFVM            # https://github.com/j-fu/VoronoiFVM.jl

include("auxiliary.jl");            export set_parameters!, ishless, row2dict
include("cell.jl");                 export AbstractCell, CommonCell #, inival, stationary_update!, update_parameters!
include("grids.jl");                export half_electrolyte_1D, full_electrolyte_1D, cell1D, simplefullcell1D
include("finite_volume_fluxes.jl");   export LGS_flux, Fick_flux, PB_flux
include("units.jl");                export kB, Planck_constant, mₑ, e0, ε0, N_Avo, T0, a0, nm

include("./cells/AuYSZAu.jl")       
# cells
export YHElectrolyte, AYABoltzmann, AYALGBoltzmann, AYALG1iBoltzmann
# methods 
export inival, biasshow!, biastest!, stationary_update!, update_parameters!, AuL_charge

end
