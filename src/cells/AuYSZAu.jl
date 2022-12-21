#=
Domain annotation
=#
const Ω_Aul, Ω_YSZ, Ω_Aur = 1, 2, 3
const bulk_domains = (Ω_Aul, Ω_YSZ, Ω_Aur)
const Γ_Aul, Γ_YSZl, Γ_YSZr, Γ_Aur = 1, 3, 4, 2
const interfaces = (Γ_Aul, Γ_YSZl, Γ_YSZr, Γ_Aur)
const ISR = (Γ_YSZl, Γ_YSZr)
const ipsi, iyV = 1, 2
const bulk_species = (ipsi, iyV)
const iyVs = 3
const surface_species = (iyVs)
const species_names = ("ipsi", "iyV")

#=
Grid parameters
=#
const electrolyte_length = 1.0e-8 # YHElectrolyte
const electrolyte_thickness = 1e-3
const electrode_thickness = 1e-5
const dmin = 1.0e-11
#
const nu = 0.93 # backward compatibility for the original model

#=
YSZ parameters
=#
# lattice
const nLYSZ = 7.46268656716418e27
const x_frac = 0.13
const mC = 4.0
const mA = 8.0
# charges and electrostatics
const epsYSZ = 28.0 * ε0
const zV = 2.0
const zi = -zV
""" lattice charge density in YSZ """
zCf(x) = 4.0 * (1.0 - x) / (1.0 + x) + 2.0 * 3.0 * x / (1.0 + x)
const Nv = (1.0 - nu) * mA * nLYSZ # upper bound of vacancy concentration 
const zC = zCf(x_frac)
const zL = mC * zC - 2 * mA

YSZ_charge_density(nV, x) = e0 * nLYSZ * (mC * zCf(x) + mA * zi) + e0 * zV * nV # TODO promote x_frac to YSZ parameters
YSZ_charge_density(nV) = YSZ_charge_density(nV, x_frac)#e0*nLYSZ*(mC*zC + mA*zi + zV*nV/nLYSZ)
electroneutral_nV_YSZ(x) = -(mC * zCf(x) + mA * zi) * nLYSZ / zV
const electroneutral_nV = electroneutral_nV_YSZ(x_frac)
X(dpsi, sys) = (electroneutral_nV / (nVmax(sys) - electroneutral_nV)) * exp(zV * e0 / kB / sys.physics.data.T * dpsi)
yVeq(x) = x ./ (1 .+ x)
nVeq(dpsi, sys) = nVmax(sys) * yVeq(X(dpsi, sys))

yVmax(a::Float64) = (-zL / zV * (1 - a) + mA * a)
yVmax(sys::VoronoiFVM.AbstractSystem) = yVmax(sys.physics.data.alpha)
nVmax(a::Float64) = nLYSZ * yVmax(a)
nVmax(sys::VoronoiFVM.AbstractSystem) = nVmax(sys.physics.data.alpha)

@composite mutable struct YHElectrolyte <: AbstractCell
    CommonCell...
end
function YHElectrolyte(grid=half_electrolyte_1D(electrolyte_length, dmin))
    new = YHElectrolyte(Model((SharedParams, YSZparams)), bulk_lattice_half_cell(grid=grid), 0.0, 0.0)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end

@kwdef struct SharedParams{A,B}
    T::A = Param(T0 + 800, bounds=(600, 900), name = "temperature")
    bias::B = Param(0.0, bounds=(-3.0, 3.0))
end
@kwdef struct YSZparams{A,B,C}
    AYSZ::A = Param(0.0 * 1.0e-3, bounds=(1e-5, 1e2))
    DYSZ::B = Param(1.0e-10, bounds=(1e-5, 1e2))
    alpha::C = Param(0.05, bounds=(0.0, 1.0))
end
@kwdef struct Auparams{A}
    DAu::A = Param(1.0e-11, bounds=(1e-5, 1e2))
end
@kwdef mutable struct parametersY
    AYSZ::Float64 = 0.0 * 1.0e-3 # rename to AYSZ in YHElectrolyte
    DYSZ::Float64 = 1.0e-10 # rename to DYSZ in YHElectrolyte
    Nv::Float64 = Nv# 1.0e-3 # TODO introduce the alpha reparametrization
    #
    T::Float64 = T0 + 800
    #	
    bias::Float64 = 0.0
end

function update_sys_parameters!(sys, params_dict::Dict)
    set_parameters!(sys.physics.data, params_dict)
    boundary_dirichlet!(sys, ipsi, 1, sys.physics.data.bias)
end

update_parameters!(cell::YHElectrolyte, params_dict::Dict) = update_sys_parameters!(cell.system, params_dict)

function flux_test(sys)
    edge = nothing
    data = sys.physics.data
    flux_wrapper!(f, u) = sys.physics.flux!(f, u, edge, data)
    u_kl = [1.0 2.0; 3.0 2.0; 2.0 1.0]
    f_k = [0.8, 0.8, 0.8]
    @show ForwardDiff.jacobian(flux_wrapper!, f_k, u_kl)
end

function YSZfluxes!(f, u, edge, data)
    f[iyV] = LGS_flux(u[iyV, 1], u[iyV, 2], u[ipsi, 1], u[ipsi, 2], data.DYSZ, e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T) # interaction energy A enters here
    f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsYSZ)
end

function YSZstorage!(f, u, node, data)
    f[iyV] = u[iyV]
    f[ipsi] = 0.0
end

function YSZreaction!(f, u, node, data)
    f[iyV] = 0.0
    f[ipsi] = -YSZ_charge_density(data.Nv * u[iyV])
end


function bulk_lattice_half_cell(; grid=nothing)
    grid = half_electrolyte_1D(electrolyte_length, dmin)  # TODO verify the indexing of the grid
    data = parametersY()
    physics = VoronoiFVM.Physics(
        data=data,
        num_species=2,
        flux=YSZfluxes!,
        #
        storage=YSZstorage!,
        ##
        reaction=YSZreaction!
    )
    sys = VoronoiFVM.System(grid, physics, unknown_storage=:sparse)


    enable_species!(sys, iyV, [Ω_YSZ])
    enable_species!(sys, ipsi, [Ω_YSZ])
    # electroneutral Dirichlet and zero Neumann for iyV
    boundary_dirichlet!(sys, iyV, 2, electroneutral_nV / Nv)
    # Dirichlet BCs for ipsi
    boundary_dirichlet!(sys, ipsi, 1, data.bias)
    boundary_dirichlet!(sys, ipsi, 2, 0.0)

    return sys
end

function inival(cell::YHElectrolyte)
    inival = VoronoiFVM.unknowns(cell.system)
    inival[ipsi, :] .= 0.0
    inival[iyV, :] .= electroneutral_nV / nLYSZ
    return inival
end

function equilibrium_solution(sys)
    inival = VoronoiFVM.unknowns(sys)
    inival[ipsi, :] .= 0.0
    inival[iyV, :] .= electroneutral_nV / nLYSZ
    #
    update_sys_parameters!(sys, Dict(:bias => 0.0))
    VoronoiFVM.solve!(inival, inival, sys, control=control)
    return inival
end
function stationary_update!(cell::YHElectrolyte, params_dict; tend=1e-3)
    reload = false
    for par in [:alpha, :alphas]
        if (par ∈ keys(params_dict))
            if (getfield(cell.system.physics.data, par) != params_dict[par])
                reload = true
            end
        end
    end
    update_parameters!(cell, params_dict)
    if reload
        cell.U = inival(cell)
    end
    cell.Uold = cell.U
    VoronoiFVM.solve!(cell.U, cell.U, cell.system, control=control)
end


function stationary_update!(inival, sys, params_dict)
    update_sys_parameters!(sys, params_dict)
    VoronoiFVM.solve!(inival, inival, sys, control=control)
end

function current_integrator(sys, U)
    factory = VoronoiFVM.TestFunctionFactory(sys)
    tfL = testfunction(factory, 2, 1)
    return VoronoiFVM.integrate_stdy(sys, tfL, U)
end

function impedance_current(sys, U)
    return current_integrator(sys, U)[ipsi]
end





#=
Au (gold) domain parameters Electrostatic parameters
=#
const epsAu = 6.9 * ε0 # Separation of the contribution of free and bound electrons into real and imaginary parts of the dielectric constant of gold.     Shklyarevskii, I. N.; Pakhmov, P. L. USSR. Optika i Spektroskopiya  (1973),  34(1),  163-6.
const ze = -1.0

#=
Lattice densities
=#
# const nLYSZ = nL
const nLAu = 1 / (4 * pi / 3 * (3.01 * a0)^3) # RM Martin, Electronic structure (eq 5.1, )


Au_charge_density(ne) = e0 * (nLAu * 1.0 + ze * ne)


# electroneutral concentrations
const neutral_eAu = 1.0
const electroneutral_nAu = 1.0 * nLAu


# stable configuration of the solver
const testing = true
const control = VoronoiFVM.NewtonControl()
#control.tol_absolute = 1e-4
#control.tol_relative = 1e-5
control.max_iterations = 100
# control.damp_initial = 1e-6
# control.damp_growth = 1.15
# control.verbose=true
# control.force_first_step = true
# control.handle_exceptions= true
control.Δt = 1.0e-4
control.Δt_min = 1.0e-50
control.Δt_max = 1.0
control.Δu_opt = 0.1 # smaller if imprecise
control.verbose = false#testing ? true : false
control.in_memory = false
control.store_all = false

# @kwdef mutable struct parameters
#     AYSZ::Float64 = 0.0 * 1.0e-3
#     DYSZ::Float64 = 1.0e-10
#     DAu::Float64 = 1.0e-11
#     alpha::Float64 = 0.05
#     #
#     T::Float64 = T0 + 800
#     #	
#     bias::Float64 = 0.0
# end


#=

The full transient gPNP problem for the blocking Au|YSZ|Au cell

=#
@composite mutable struct AYABoltzmann <: AbstractCell
    CommonCell...
end

parameters2VoronoiData(Model((SharedParams(), YSZparams(), Auparams())), :AYABoltzmannData)()

function AYABoltzmann(grid=cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-12))
    new = AYABoltzmann(Model((SharedParams(), YSZparams(), Auparams())), ayasystem(grid)[1], 0.0, 0.0)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end


function ayasystem(grid)
    system = VoronoiFVM.System(grid, unknown_storage=:sparse)
    psi = VoronoiFVM.ContinuousQuantity(system, collect(bulk_domains))
    dspec = VoronoiFVM.DiscontinuousQuantity(system, collect(bulk_domains))
    function flux(f, u, edge, data)
        if edge.region == 1
            f[psi] = Fick_flux(u[psi, 1], u[psi, 2], epsAu)
            f[dspec] = PB_flux(u[dspec, 1], u[dspec, 2], u[psi, 1], u[psi, 2], data.DAu, e0 * ze / kB / data.T)
        elseif edge.region == 2
            f[psi] = Fick_flux(u[psi, 1], u[psi, 2], epsYSZ)
            f[dspec] = LGS_flux(u[dspec, 1], u[dspec, 2], u[psi, 1], u[psi, 2], data.DYSZ, e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T)
        elseif edge.region == 3
            f[psi] = Fick_flux(u[psi, 1], u[psi, 2], epsAu)
            f[dspec] = PB_flux(u[dspec, 1], u[dspec, 2], u[psi, 1], u[psi, 2], data.DAu, e0 * ze / kB / data.T)
        end
    end

    function storage(f, u, node, data)
        f[psi] = 0.0
        f[dspec] = u[dspec]
    end

    function reaction(f, u, node, data)
        if node.region == 1
            f[psi] = -Au_charge_density(nLAu * u[dspec])
        elseif node.region == 2
            f[psi] = -YSZ_charge_density(nVmax(data.alpha) * u[dspec])
        elseif node.region == 3
            f[psi] = -Au_charge_density(nLAu * u[dspec])
        end
    end

    VoronoiFVM.physics!(system,
        VoronoiFVM.Physics(
            # data=parameters(),
            data=AYABoltzmannData(),
            flux=flux,
            reaction=reaction,
            storage=storage,
        )
    )
    boundary_dirichlet!(system, psi, 1, 0.0)
    boundary_dirichlet!(system, psi, 2, 0.0)
    #
    boundary_dirichlet!(system, dspec, 1, neutral_eAu)
    boundary_dirichlet!(system, dspec, 2, neutral_eAu)
    #
    check_allocs!(system, true)
    return system, psi, dspec
end

AYAsystem() = ayasystem()[1]

function aya_inival(system)
    inival = VoronoiFVM.unknowns(system, inival=nVmax(0.0) / nVmax(system))
    inival[1, :] .= 0.0
    if size(inival)[1] > 2
        inival[2, :] .= 1.0
        inival[4, :] .= 1.0
    end
    return inival
end
inival(cell::AYABoltzmann) = aya_inival(cell.system)
#=
Reduced problem: 
- YSZ -- transient gPNP for ionic vacancies
- Au -- nonlinear Poisson problem

=#
# FIXME solve the math in #= ... =# comments
@composite mutable struct AYALGBoltzmann <: AbstractCell
    CommonCell...
end

function AYALGBoltzmann(grid=cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-12))
    new = AYALGBoltzmann(Model((SharedParams(), YSZparams(), Auparams())), ayasystemLG(grid, AueDensity=BoltzmannAue), 0.0, 0.0)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end


ayaLGphys = function (; AueDensity=BoltzmannAue)
    return VoronoiFVM.Physics(
        # data=parameters(),
        data=AYABoltzmannData(),
        #
        flux=function (f, u, edge, data)
            if edge.region == 1
                f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsAu)
            elseif edge.region == 2
                #f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsYSZ)
                #f[iyV] = LGS_flux(u[iyV, 1], u[iyV, 2], u[ipsi, 1], u[ipsi, 2], data.DYSZ, e0 * zV / kB / data.T, data.AYSZ * e0 / kB / data.T)
                YSZfluxes!(f,u,edge,data)
            elseif edge.region == 3
                f[ipsi] = Fick_flux(u[ipsi, 1], u[ipsi, 2], epsAu)
            end
        end,
        #
        reaction=function (f, u, node, data)
            if node.region == 1
                f[ipsi] = -Au_charge_density(AueDensity(data, data.bias - u[ipsi]))
            elseif node.region == 2
                f[ipsi] = -YSZ_charge_density(nVmax(data.alpha) * u[iyV])
            elseif node.region == 3
                f[ipsi] = -Au_charge_density(AueDensity(data, -u[ipsi]))
            end
        end,
        #
        storage=function (f, u, node, data)
            f[ipsi] = 0.0
            if node.region == 2
                #f[iyV] = u[iyV]
                YSZstorage!(f, u, node, data)
            end
        end,
    )
end

BoltzmannAue(data, V) = nLAu * exp(ze * e0 / kB / data.T * V)

# Thomas-Fermi-Dirac model of electrons in metal 
# TODO check -> correct -> test
# const Ck = 3*Planck_constant^2/40/mₑ/a0^2*(3/pi)^(2/3) 
# TFDAue(data, V) = (nLAu^(2/3) + V/Ck)^(3/2) 

function ayasystemLG(grid; AueDensity=BoltzmannAue)
    system = VoronoiFVM.System(grid, ayaLGphys(AueDensity=AueDensity), unknown_storage=:sparse)
    #
    enable_species!(system, ipsi, collect(bulk_domains))
    enable_species!(system, iyV, [Ω_YSZ])
    #
    boundary_dirichlet!(system, ipsi, 1, 0.0)
    boundary_dirichlet!(system, ipsi, 2, 0.0)
    #
    check_allocs!(system, true)
    return system
end

function ayaLGinival(system)
    inival = VoronoiFVM.unknowns(system, inival=nVmax(0.0) / nVmax(system))
    inival[ipsi, :] .= 0.0
    return inival
end

inival(cell::AYALGBoltzmann) = ayaLGinival(cell.system)

#=

the previous case extended with interface species
chemical potential of the surface vacancies

=#
# TODO governing eqns for the interface species
@composite mutable struct AYALG1iBoltzmann <: AbstractCell
    CommonCell...
end


function AYALG1iBoltzmann(grid=cell1D(electrode_thickness, electrolyte_thickness, electrode_thickness, dmin=1e-12))
    new = AYALG1iBoltzmann(0.0, 0.0, 0.0, 0.0)
    new.parameters = Model((SharedParams(), YSZparams(), Auparams(), ISRparameters()))
    new.system = ayasystemLG1i(grid, AueDensity=BoltzmannAue)
    new.U = inival(new)
    new.Uold = copy(new.U)
    return new
end


const nLYSZs = nLYSZ^(2 / 3)  # [# YSZ faces/m^2]
const nLAus = nLAu^(2 / 3)    # [# Au faces/m^2]
const mCs = 2               # [# cation sites/YSZ face]
const mAs = 4               # [# anion sites/YSZ face]
const zLYSZs = mCs * zC - 2 * mAs # [static elementary charge/YSZ face] <=> fully occupied oxide vacancies
const ISR_staticcharge = nLYSZs * zLYSZs + 1 * nLAus                          # per ISR area [C/m^2]
# ISR_chargedensity(nes, nVs, areaRatio::Float64) = e0 * areaRatio * (nes * (-1.0) + nVs * zV + ISR_staticcharge) # per cell cross-section [C/m^2]
ISR_chargedensity(nes, nVs, areaRatio::Float64) = e0 * (-nes + nVs * zV +  areaRatio * ISR_staticcharge) # per cell cross-section [C/m^2]

yVmaxs(a::Float64)::Float64 = (-zLYSZs / zV * (1 - a) + mAs * a) # = nVmaxs / S_l nLYSZs = (mAs - (1-a)*mCs*zC/zV) # checked
yVmaxs(sys::VoronoiFVM.AbstractSystem)::Float64 = yVmaxs(sys.physics.data.alpha)
nVmaxs(a::Float64, S::Float64)::Float64 = nLYSZs * S * yVmaxs(a)                                       # per ISR area
#nVmaxs(sys::VoronoiFVM.AbstractSystem)::Float64 = nVmaxs(sys.physics.data.alpha)     # per ISR area

@kwdef struct ISRparameters{A,B,C,D,E,F,G,H,I}
    alphas::A = Param(0.05, bounds=(0.0,1.0)) # [1] ratio of admissible vacancies at ISR
    AYSZs::B = Param(0.0, bounds=(0.0,1.0)) # [eV] interaction energy of vacancies at ISR
    GA::C = Param(0.0 * e0, bounds=(0.0,1.0)) # [eV] Gibbs energy of vacancy adsorption
    Ge::D = Param(0.0 * e0, bounds=(0.0,1.0)) # [eV] Gibbs energy of electron adsorption
    kA::E = Param(1.0e15, bounds=(0.0,1.0)) # [1/m^2/s] rate of vacancy adsorption
    dL::F = Param(1.0 * nm, bounds=(0.0,1.0)) # [nm] thickness of the left ISR
    dR::G = Param(1.0 * nm, bounds=(0.0,1.0)) # [nm] thickness of the right ISR
    SL::H = Param(1.0, bounds=(0.0,1.0)) # [1] (ISR area)/(cross section area) -- left
    SR::I = Param(1.0, bounds=(0.0,1.0)) # [1] (ISR area)/(cross section area) -- right
end
parameters2VoronoiData(Model((SharedParams(), YSZparams(), Auparams(), ISRparameters())), :AYALG1iBoltzmannData)()
# @composite @kwdef mutable struct parameters1i
#     parameters...
#     alphas::Float64 = 0.05 # [1] ratio of admissible vacancies at ISR
#     AYSZs::Float64 = 0.0 # [eV] interaction energy of vacancies at ISR
#     GA::Float64 = 0.0 * e0 # [eV] Gibbs energy of vacancy adsorption
#     Ge::Float64 = 0.0 * e0 # [eV] Gibbs energy of electron adsorption
#     kA::Float64 = 1.0e15 # [1/m^2/s] rate of vacancy adsorption
#     dL::Float64 = 1.0 * nm # [nm] thickness of the left ISR
#     dR::Float64 = 1.0 * nm # [nm] thickness of the right ISR
#     SL::Float64 = 1.0 # [1] (ISR area)/(cross section area) -- left
#     SR::Float64 = 1.0 # [1] (ISR area)/(cross section area) -- right
# end

ISR_arearatio(bnode,data) = (bnode.region == Γ_YSZl ? data.SL : data.SR)
reduced_voltage(u, bnode, data) = (bnode.region == Γ_YSZl ? data.bias : 0.0) - u[ipsi]
ISR_electrondensity(u, bnode, data) = nLAus / nLAu * exp(-data.Ge / kB / data.T) * BoltzmannAue(data, reduced_voltage(u,bnode,data)) # [# electrons/ ISR area]
ISR_chargedensity(u, bnode, data) = ISR_chargedensity(ISR_electrondensity(u, bnode, data), u[iyVs], ISR_arearatio(bnode,data))  # per cell cross-section [C/m^2]

ayaLGphys1i = function (; AueDensity=BoltzmannAue)
    ayaLGp = ayaLGphys(AueDensity=AueDensity)
    return VoronoiFVM.Physics(
        #data=parameters1i(),
        data=AYALG1iBoltzmannData(),
        # bulk properties inherited from ayaLGphysics
        flux=ayaLGp.flux,
        reaction=ayaLGp.reaction,
        storage=ayaLGp.storage,
        breaction=function (f, u, bnode, data)
            if bnode.region in ISR
                # !! notation: ν is an outer normal vector
                # Adsorption of vacancies 
                # TODO add the difference of the electrostatic potential or its derivative*thickness of the ISR to the "equilibrium constant for vacancies"
                KVsq = exp(-data.GA / data.T / kB / 2 + (data.AYSZ * u[iyV] - data.AYSZs * u[iyVs]) / 2) # sqrt(KV)
                ReducedRateA = KVsq * u[iyV] * (1.0 - u[iyVs]) - 1 / KVsq * u[iyVs] * (1.0 - u[iyV])
                # INFO the reaction rate is calculated in [# vacancies/cross section area]
                # and its contribution to both coverages isn't thus far scaled appropriately
                # However, this is no problem for a blocking electrode in equilibrium...
                # the fix should look something like
                # ... a check is needed though
                # AreaRatio = (bnode.region == Γ_YSZl ? data.SL : data.SR)
                AreaRatio = ISR_arearatio(bnode, data)
                # implementation for bspecies 
                # bstorage + breaction = 0
                f[iyVs] = -data.kA * ReducedRateA / nVmaxs(data.alpha, AreaRatio)
                # implementation for species
                # - j ̇ν + breaction = 0
                f[iyV] = data.kA * ReducedRateA / nVmax(data.alpha)
                # equilibrium of electrons
                # V = (bnode.region == Γ_YSZl ? data.bias : 0.0) - u[ipsi]
                V = reduced_voltage(u, bnode, data) 
                # TODO add the difference of the electrostatic potential or its derivative*thickness of the ISR to the "equilibrium constant for electrons"
                # FIXME implicitly assuming the Boltzmann statistics for the ISR electrons
                # nes = nLAus / nLAu * exp(-data.Ge / kB / data.T) * AueDensity(data, V) # [# electrons/ ISR area]
                nes = ISR_electrondensity(u, bnode, data) # [# electrons/ ISR area]
                # TODO ISRthickness = (bnode.region  == Γ_YSZl ? data.dL : data.dR)
                ###
                # The surface Poisson equation is consistent with [BSE2018]
                # (-ε_+ ∇ψ_+)⋅ν_+ (-ε_- ∇ψ_-)⋅ν_- + ISR_chargedensity = 0
                # However, the equation for normal fluxes and the breaction for an internal node between two REVs is implemented as
                # - (j_+ ̇ν_+ + j_- ̇ν_-) + breaction = 0 
                # thus ISR_chargedensity below, correctly, enters with a negative sign !!!
                f[ipsi] = -ISR_chargedensity(nes, nVmaxs(data.alpha, AreaRatio ) * u[iyVs], AreaRatio)
                #
            end
        end,
        #
        bstorage=function (f, u, bnode, data)
            if bnode.region in ISR #[Γ_YSZl,Γ_YSZr]
                f[iyVs] = u[iyVs]
            end
        end
    )
end

# ayasystemLGi1 constructor
function ayasystemLG1i(grid; AueDensity=BoltzmannAue)
    system = VoronoiFVM.System(grid, ayaLGphys1i(AueDensity=AueDensity), unknown_storage=:sparse)
    #
    enable_species!(system, ipsi, collect(bulk_domains))
    enable_species!(system, iyV, [Ω_YSZ])
    enable_boundary_species!(system, iyVs, [Γ_YSZl, Γ_YSZr])
    #
    boundary_dirichlet!(system, ipsi, 1, 0.0)
    boundary_dirichlet!(system, ipsi, 2, 0.0)
    #
    check_allocs!(system, true)
    return system
end
# initial values for ayasystemLGi1
function ayaLG1iinival(system)
    inival = VoronoiFVM.unknowns(system, inival=nVmax(0.0) / nVmax(system))
    inival[ipsi, :] .= 0.0
    inival[iyVs, :] .= yVmaxs(0.0) / yVmaxs(system) # TODO adjust for each boundary
    return inival
end

inival(cell::AYALG1iBoltzmann) = ayaLG1iinival(cell.system)

function biasshow!(cell::AbstractCell; bend=1.0, bstep=0.1, tend=1e-3)
    cell.U = inival(cell)
    cell.Uold = copy(cell.U)
    #results = DataFrame(bias=Float64[], U=typeof(cell.U))
    results = DataFrame(bias=Float64[], U=Any[])
    for bias ∈ collect(0.0:bstep:bend)
        stationary_update!(cell, Dict(:bias => bias), tend=tend)
        push!(results, Dict(:bias => bias, :U => copy(cell.U)))
    end
    cell.U = copy(cell.Uold)
    stationary_update!(cell, Dict(:bias => 0.0), tend=tend)
    for bias ∈ collect(-bstep:-bstep:-bend)
        stationary_update!(cell, Dict(:bias => bias), tend=tend)
        push!(results, Dict(:bias => bias, :U => copy(cell.U)))
    end
    return sort!(results, [:bias])
end



# function biassweep(system::VoronoiFVM.System, biasend) 
#     biasstart = system.physics.data.bias
#     U = VoronoiFVM.unknowns(system, inival=biasend)
#     inival= VoronoiFVM.unknowns(system, inival=0.5*biasend)
#     for bias ∈ collect(biasstart:0.1*(biasend-biasstart):biasend)
#         @show bias
#         update_parameters!(system, Dict(:bias => bias))
#         VoronoiFVM.solve!(U, U, system, control=control)
#     end
#     println("sys swept to ",biasend," V")
#     return U
# end

AuL_charge(cell::AbstractCell) = AuL_charge(cell.system, cell.U)

function AuL_charge(system::VoronoiFVM.AbstractSystem, state)
    QrAuL = 0.0
    QrAuL = VoronoiFVM.integrate(system, system.physics.reaction, state)
    return -QrAuL[ipsi, Ω_Aul]
end

function AuL_charge(cell::AYALG1iBoltzmann)
    alternative = true
    QrAuL = VoronoiFVM.integrate(cell.system, cell.system.physics.reaction, cell.U) # nspec x nregion
    # bulk, bnd = bulk_domains[1], ISR[1]  # Γ_YSZl 
    bulk, bnd = bulk_domains[3], ISR[2]  # Γ_YSZr 
    bnd_index = cell.system.grid.components[BFaceNodes][bnd]
    data = cell.system.physics.data
    if alternative
        AreaRatio = (bnd == Γ_YSZl ? data.SL : data.SR)
        V = (bnd == Γ_YSZl ? data.bias : 0.0) - cell.U[ipsi, bnd_index]
        nes = nLAus / nLAu * exp(-data.Ge / kB / data.T) * BoltzmannAue(data, V) # [# electrons/ ISR area]
        ISRcontribution = e0*nes
        #ISRcontribution = -ISR_chargedensity(nes, nVmaxs(data.alpha, AreaRatio) * cell.U[iyVs, bnd_index], AreaRatio)
    else
        QrAuLb = VoronoiFVM.integrate(cell.system, cell.system.physics.breaction, cell.U, boundary=true) # nspec x nregion
        ISRcontribution = QrAuLb[ipsi, bnd]
        # TODO try to integrate breaction instead of reaction !!
        # FIXME the result of the boundary integration is really large for Γ_YSZl, but 0.0 for Γ_YSZr. This needs to be resolved.
        # @show QrAuLb[ipsi, Γ_YSZl:Γ_YSZr] 
    end
    return -QrAuL[ipsi, bulk] + ISRcontribution # the sign of the surface charge is opposite, see surface Poisson
end


function stationary_update!(cell::AbstractCell, params_dict; tend=1e-3)
    reload = false
    for par in [:alpha, :alphas]
        if (par ∈ keys(params_dict))
            if (getfield(cell.system.physics.data, par) != params_dict[par])
                reload = true
            end
        end
    end
    update_parameters!(cell, params_dict)
    if reload
        cell.U = inival(cell) # no re-allocations
    end
    tsol = VoronoiFVM.solve(
        cell.U,
        cell.system,
        [0.0, tend],
        control=control
    )
    cell.U .= tsol[end]
end

function update_parameters!(sys::VoronoiFVM.AbstractSystem, params_dict::Dict)
    set_parameters!(sys.physics.data, params_dict)
    boundary_dirichlet!(sys, ipsi, Γ_Aul, sys.physics.data.bias)
end

update_parameters!(cell, params_dict::Dict) = update_parameters!(cell.system, params_dict::Dict)