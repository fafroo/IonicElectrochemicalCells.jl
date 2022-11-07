abstract type AbstractCell end

mutable struct CommonCell <: AbstractCell
    parameters#::ModelParameters.Model
    system#::VoronoiFVM.AbstractSystem
    U#::VoronoiFVM.AbstractUnknown
    Uold#::VoronoiFVM.AbstractUnknown
end

# function inival(cell :: AbstractCell)
#     cell.U = VoronoiFVM.unknowns(cell.system, inival=0.0)
#     cell.Uold = VoronoiFVM.unknowns(cell.system, inival=0.0)
# end
# 
# function update_sys_parameters!(cell::AbstractCell, params_dict::Dict)
# 	set_parameters!(cell.system.physics.data, params_dict)
# end

# function stationary_update!(cell::AbstractCell, params_dict::Dict; tend=1e-3)
# 	update_sys_parameters!(cell, params_dict)
#     cell.Uold .= cell.U
#     tsol = VoronoiFVM.solve(
#             cell.U,
#             cell.system, 
#             [0.0, tend]
#             , control=control
#         )
#     cell.U .= tsol[end]
# end


function biastest!(cell :: AbstractCell)
    for bias ∈ collect(0.0:-0.1:-1.0)
        stationary_update!(cell , Dict(:bias => bias))
    end
end

