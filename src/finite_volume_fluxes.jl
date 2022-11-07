"""
Q1: what is a reasonable scaling for fluxes?
- per the crystalline lattice density
- as an electric current density 
- molar fluxes
"""
# TODO adjust the multiline comments
function Fick_flux(psi1, psi2, eps)
    return -eps*(psi2 - psi1)
end  

function sedan_flux(u1, u2, D, g)
	# flux: j = - ( u' + u*q)
	# dirichlet BCs u1        u2 
	# line segment  |---------|
	# colloc. pts  x1        x2
	# h = |x2-x1|
	# g = q*h
	bp,bm=VoronoiFVM.fbernoulli_pm(g) # Bernoulli function of positive and negative argument
	return -1.0*D*(u2*bm - u1*bp) 
end

function PB_flux(u1, u2, psi1, psi2, D=1.0, UT=1.0)
	g = UT*(psi2 - psi1) # UT = zV*e0/kB/T
    sedan_flux(u1, u2, D, g)
end

#=
## Free energy
```math
\frac{F}{k_B T} = n\mu^{ref} + n\log n + (nMax-n)\log(nMax-n) + A/2 \frac{n^2/ nMax}
```
## Chemical potential
```math
\frac{mu}{k_B T}= \mu^{ref} + \log \frac{n}{nMax-n} + A \frac{n/ nMax}
```

`` y = \frac{n}{nMax} ``
## continuous Flux
```math
\textbf{j}{nMax} 
	&= - D y(1-y)\left(
			\frac{\nabla y}{y(1-y)}
			+ U_T \nabla\psi 
			+ A \nabla y
		\right)
    &= - D \left(
			\nabla y
			+ y\left[ 
				(1-y)\left(
					U_T \nabla\psi 
					+ A \nabla y
				\right)
			\right]
		\right)
```
## Argument `g` of the Sedan flux in `LGS_flux``
```math
g = \textbf{Qh}
	=\left(
		1 - \frac{y_2 + y_1}{2} 
	\right)
	\left(
		U_T (\psi_2 - \psi_1)
		+ A (y_2 - y_1)
	\right)```
=#
function LGS_flux(u1, u2, psi1, psi2, D=1.0, UT=1.0, A=0.0)
	#  g = Q*h
	g = (1.0 - 0.5*(u1 + u2))*(
		  A*(u2 - u1)
		 +UT*(psi2 - psi1) # UT = zV*e0/kB/T
	)
	#return -1.0*D*(u2*bm - u1*bp)# y_V^R B(-gh) - y_V^L B(gh)
    return sedan_flux(u1, u2, D, g)
end # interactive_vac_flux


function legacy_ion_flux(u1, u2, psi1, psi2, D, UT) 
    mu1=log(1-u1)
    mu2=log(1-u2)
    g = +(mu2-mu1) - UT*(psi2-psi1)# UT = e0*za/T/kB
    return D*sedan(g, u1, u2)
end # legacy_ion_flux



# export LGS_flux, Fick_flux, PB_flux, # sym_lat_gas_flux, legacy_ion_flux, , PB_flux2
