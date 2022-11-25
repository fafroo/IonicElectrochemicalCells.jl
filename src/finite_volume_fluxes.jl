"""

Q1: what is a reasonable scaling for fluxes?
- per the crystalline lattice density
- as an electric current density 
- molar fluxes
"""

"""
$(SIGNATURES)

Return integrated Fick's diffusion flux
"""
function Fick_flux(psi1, psi2, eps)
    return -eps*(psi2 - psi1)
end  

"""
$(SIGNATURES)

Return Sedan approximation of flux ``\v j``



"""
function sedan_flux(u1, u2, D, g)

	bp,bm=VoronoiFVM.fbernoulli_pm(g) # Bernoulli function of positive and negative argument
	return -1.0*D*(u2*bm - u1*bp) 
end

"""
$(SIGNATURES)

Return integrated Poisson-Boltzmann drift-diffusion flux

# Sedan-scheme for
```math
\\int_0^\\infty = 
```
"""
function PB_flux(u1, u2, psi1, psi2, D=1.0, UT=1.0)
	g = UT*(psi2 - psi1) # UT = zV*e0/kB/T
    sedan_flux(u1, u2, D, g)
end



#=
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
	\right)
```
=#
"""
$(SIGNATURES)

Return integrated Sedan approximation of lattice-gas-drift-diffusion flux

## continuous Flux
```math
\\begin{aligned}
\\textbf{j}{nMax} 
	&= - D y(1-y)\\left(
			\\frac{\\nabla y}{y(1-y)}
			+ U_T \\nabla\\psi 
			+ A \\nabla y
		\right)
    &= - D \\left(
			\\nabla y
			+ y\\left[ 
				(1-y)\\left(
					U_T \\nabla\\psi 
					+ A \\nabla y
				\\right)
			\\right]
		\\right)
\\end{aligned}
```

"""
function LGS_flux(u1, u2, psi1, psi2, D=1.0, UT=1.0, A=0.0)
	#  g = Q*h
	g = (1.0 - 0.5*(u1 + u2))*(
		  A*(u2 - u1)
		 +UT*(psi2 - psi1) # UT = zV*e0/kB/T
	)
	#return -1.0*D*(u2*bm - u1*bp)# y_V^R B(-gh) - y_V^L B(gh)
    return sedan_flux(u1, u2, D, g)
end # LGS_flux

function legacy_ion_flux(u1, u2, psi1, psi2, D, UT) 
    mu1=log(1-u1)
    mu2=log(1-u2)
    g = +(mu2-mu1) - UT*(psi2-psi1)# UT = e0*za/T/kB
    return -sedan_flux(u1, u2, D, g)
end # legacy_ion_flux