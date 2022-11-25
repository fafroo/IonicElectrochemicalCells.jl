# Capacitance of blocking YSZ|Au electrode
The package contains four model instances:  `YHElectrolyte`,`AYABoltzmann`, `AYALGBoltzmann` and `AYALG1iBoltzmann`.

- TODO hierarchy diagramm 


# Vacancies in Yttria-stabilized Zirconia electrolyte
#### Model of crystalline structure and charge density
Yttria-stabilized zirconia (YSZ) is a solid mixture of yttrium oxide and zirconium oxide: ``x^{\#}\mathrm{Y}_2\mathrm{O}_3 + (1-x^{\#})\mathrm{Zr}\mathrm{O}_2`` for ``x^{\#} \in (0,1)``. 
```@raw html
<center>
<img src="../res/ysz_m.png" width="50%">
</center>
```
In face-centered cubic YSZ, each unit cell contains ``m_\text{C}=4`` cation sites and ``m_\text{A} = 8`` anion sites.
The cation sites are occupied either by ``\mathrm{Zr}^{\textrm{IV}}`` or ``\mathrm{Y}^{\textrm{III}}`` ions and formally carry 
```math
z_\text{C}(x^{\#}) = 4\frac{1 - x^{\#}}{1 + x^{\#}} + 2\cdot 3\frac{x^{\#}}{1 + x^{\#}}
```
elementary charge per cation site.

The anion sites can be occupied by doubly-negatively charged oxide ions ``\mathrm{O}^{\textrm{II-}}`` or be vacant which compensates the lower valency of ``\mathrm{Y}^{\textrm{III}}``.
For a *volume-averaged concentration of oxide vacancies* ``n_\text{V}``, the volumetric charge density reads 
```math
\rho_\text{YSZ}(n_\text{V}, x^{\#}) = n_\text{YSZ}^{\#}\Big(z_\text{C}(x^{\#}) m_\text{C}  - 2 m_\text{A}\Big) + z_\text{V} n_\text{V}~,
```
where ``z_\text{V}=2`` denotes elementary charge of oxide vacancy and ``n_\text{YSZ}^{\#}`` stands for the volumetric density of elementary cells.
For a given ``x^{\#}``, there's an *electroneutral oxide vacancy concentration* ``n_\text{V}^\text{eq}`` such that ``0 = \rho_\text{YSZ}(n_\text{V}^\text{eq}, x^{\#})``.
Explicitely,
```math 
\begin{aligned}
n_\text{V}^\text{eq}(x^{\#}) = n_\text{YSZ}^{\#}\Big(m_\text{A} - \frac{z_\text{C}(x^{\#})}{z_\text{V}} m_\text{C}  \Big)~.
\end{aligned}
```

The concentration of the oxide vacancies is bound from above by the number of anion lattice sites. However, reaching this bound is unphysical, because 
it would imply that the entire crystal contains only cations. To avoid this, we introduce an additional bound, *maximum conecentration of vacancies* ``n_\text{V}^\text{max}``.
Moreover, we assume ``n_\text{V}^\text{eq} \leq n_\text{V}^\text{max}`` to allow for electroneutrality. In total, we have
```math
0 
\leq n_\text{V}^\text{eq}
\leq n_\text{V}^\text{max}
\leq m_\text{A} n_\text{YSZ}^{\#}
~.
```
#### Free energy
The free energy model of the oxide vacancies contains the interacting lattice gas
```math
\frac{F}{k_B T} = n\frac{\mu^\text{ref}_\text{V}(T)}{k_\text{B}T} + \underbrace{n\log n + (n^\text{max}-n)\log(n^{max}-n)}_{\text{ideal lattice gas}} + \underbrace{\frac{A}{2} \frac{n^2}{n^\text{max}}}_{\text{interaction}}
~,
```
the chemical potential of the oxide vacancies, defined as ``\mu_\text{V}=\frac{\partial F}{\partial n_\text{V}}``, reads
```math
\frac{\mu_\text{V}}{k_B T}= \frac{\mu^\text{ref}_\text{V}(T)}{k_\text{B}T} + \log\frac{n}{n^\text{max}-n} + A \frac{n}{n^\text{max}}
```
#### Lattice-gas Poisson-Nernst-Planck
Given the *maximum concentration of vacancies* ``n_\text{V}^\text{max}``, we define *vacancy coverage* ``y_\text{V}= \frac{n_\text{V}}{n_\text{V}^\text{max}}``.
Moreover, parameter ``\alpha\in[0,1]`` interpolates the maximum concentration of vacancies ``n_\text{V}^\text{max}`` between its bounds, 
```math
n_\text{V}^\text{max}(\alpha, x^{\#}) := (1-\alpha) n_\text{V}^\text{eq}(x^{\#}) + \alpha m_\text{A} n_\text{YSZ}^{\#}
```
Hence the vacancy concentration now reads,
```math
\begin{aligned}
\text{vacancy coverage}\quad n_\text{V}(y_\text{V}, \alpha, x^{\#}) 
&= n_\text{V}^\text{max}(\alpha, x^{\#}) y_\text{V} 
\\
&= y_\text{V}n_\text{YSZ}^{\#}\Big(m_\text{A} - (1-\alpha) m_\text{C} \frac{z_\text{C}(x^{\#})}{z_\text{V}}  \Big)~. 
\end{aligned}
```
The drift-diffusion system governing the oxide vacancies in bulk YSZ, expressed in ``(y_\text{V}, \psi)`` reads
```math
\begin{aligned}
-\nabla\cdot \left(
    \nabla \varepsilon_\text{YSZ} \psi 
\right)
%- e_0 \left(
%    n^{\#}_\text{YSZ} m_\text{C} z_\text{C}(x^{\#}) 
%    - 2 n^\text{\#}_\text{YSZ} m_\text{A} 
%    + z_\text{V} n_\text{V}(y_\text{V}, \alpha, x^{\#})
%\right) = 0\\
- e_0 n^{\#}_\text{YSZ}\bigg(
     m_\text{C} z_\text{C}(x^{\#}) 
    -  m_\text{A} z_\text{V}
    + y_\text{V} \Big(
        m_\text{A} z_\text{V} 
        - (1-\alpha)m_\text{C} z_\text{C}(x^{\#})
    \Big)
\bigg) = 0\\
%
\partial_t y_\text{V} 
+ \nabla\cdot \left(
    \underbrace{
        -D_\text{YSZ}\left(
            \nabla y_\text{V}
            +y_\text{V}(1-y_\text{V})\frac{e_0 z_\text{V}}{k_\text{B} T} \nabla\psi 
            + A_\text{YSZ} \frac{e_0}{k_\text{B} T} y_\text{V}(1-y_\text{V}) \nabla y_\text{V}
        \right) 
    }_{=:\,\vec{\mathbf{j}}_\text{V}\ \text{oxide vacancy flux}}
\right)    
= 0
~.
\end{aligned}
```
#### `YHElectrolyte`
The `YHElectrolyte` implements the above-stated Lattice-gas Poisson-Nernst-Planck system on an interval ``[0,L]``, endowed with the following boundary conditions,
```math
\psi(t, 0) = V_\text{bias}, 
\
\psi(t, L) = 0,
```
and 
```math
\
\vec{\mathbf{j}}_\text{V}(t, 0) = 0, 
\
y_\text{V}(t, L) = \frac{n^\text{eq}_\text{V}}{n^\text{max}_\text{V}}
~.
```
The condition at ``x=0`` simulate a blocking interface, whereas at ``x=L`` they correspond to electronuetral bulk.

# Symmetric Au|YSZ|Au cell
## Valence electrons in Au lattice
In `AYABoltzmann`, we model the drift-diffusion of Au valence electrons using Poisson-Nernst-Planck  
```math
\begin{aligned}
\nabla\cdot\left(-\varepsilon_\text{Au}\nabla \psi \right)
- e_0 n^{\#}_\text{Au}\left(
    1  -  y_\text{e}
\right) 
&= 0
\\
\partial_t y_\text{e} 
- \nabla\cdot D_\text{Au}\left(
    \nabla y_\text{e}
    + y_\text{e}\frac{e_0 z_\text{e}}{k_\text{B} T} \nabla\psi 
\right) 
&= 0~,
\end{aligned}
```
where ``y_\textrm{e} = \frac{n_\textrm{e}}{n_\textrm{Au}^{\#}}``.

## `AYABoltzmann`
The `AYABoltzmann` represents a symmetric Au|YSZ|Au cell. It is defined on three connected intervals ``\Omega = [] \cup [] \cup []``.
The middle interval ``[]`` represents the YSZ electrolyte and hosts the Lattice-gas PNP equations.
Since the oxide vacancies cannot leave YSZ, we apply zero Neumann boundary conditions,
```math
\vec{\mathbf{j}}_\text{V}(t, ) = 0, 
\
\vec{\mathbf{j}}_\text{V}(t, ) = 0, 
~.
```
The two remaining intervals, ``[]`` and ``[]`` represent the Au electrodes and the Au PNP is deployed there. 
Since the electrons cannot enter the YSZ at the internal interfaces, we apply zero Neumann boundary conditions, 
```math
\vec{\mathbf{j}}_\text{e}(t, ) = 0, 
\,
\vec{\mathbf{j}}_\text{e}(t, ) = 0, 
~.
```
At the outer boundaries, the charge electroneutrality and voltage is applied
```math
y_\text{e}(t, ) = \frac{n^\text{eq}_\text{e}}{n^\text{\#}_\text{Au}} = 1.0,
\
y_\text{e}(t, ) = \frac{n^\text{eq}_\text{e}}{n^\text{\#}_\text{Au}} = 1.0,
\
\psi(t, 0) = V_\text{bias}, 
\
\psi(t, L) = 0,
~.
```
## `AYALGBoltzmann` -- drift-diffusion equilibrium of electrons
The `AYALGBoltzmann` model is a simplification of the `AYABoltzmann` model.
Assuming ``D_\text{Au}\to\infty`` renders the electrochemical potential of the valence Au electrons to be in quasi-equilibrium, that is, ``\nabla y_\text{e} + y_\text{e}\frac{e_0 z_\text{e}}{k_\text{B} T} \nabla\psi = 0``.
Thus in 1D, the electrostatic potential parametrizes the density of electrons using its boundary values,
```math
y_e(x) = \exp{\frac{z_e e_0}{k_B T}(\psi_{l,r} - \psi(x)) }~.
```
This we introduce into the Poisson equation and solve non-linear Poisson equation in the Au domains instead,
```math
\begin{aligned}
\nabla\cdot\left(
    -\varepsilon_\text{Au}\nabla \psi
\right) 
- e_0 n^{\#}_\text{Au}\left(
    1  - \exp{\frac{z_e e_0}{k_B T}(\psi_{l,r} - \psi(x)) }
\right) 
= 0
~.
\end{aligned}
```

# Symmetric Au|YSZ|Au cell with interface dynamics
## Au|YSZ interface specific region (ISR)

Notation: quantities related to the ISR are underlined, that is, ``\underline{\cdot}``

The interface specific region is a thin shell between YSZ and Au, the ISR
- is treated as a manifold of codimension 1
- has thickness ``\underline{d}``
- contains both YSZ and Au cations as a fixed background charge density
- contains Au electrons and YSZ oxide vacancies
- its charge density causes discontinuity of the displacement field
- the electrostatic potential remains continuous

We define 
- lattice densities ``\underline{n}_\text{YSZ}^{\#} := S_l\left(n_\text{YSZ}^{\#}\right)^\frac{2}{3}`` and ``\underline{n}_\text{Au}^{\#} := S_l\left(n_\text{Au}^{\#}\right)^\frac{2}{3}``
- ``\underline{m}_\text{C}^\text{YSZ}=2`` and ``\underline{m}_\text{A}^\text{YSZ}=4`` 

### ISR charge density -- discontinuity of displacement field
```math
\left(
    \varepsilon_\text{YSZ}\nabla\psi\cdot\vec{\nu}_\text{YSZ}
    +\varepsilon_\text{Au}\nabla\psi\cdot\vec{\nu}_\text{Au}
\right)|_{x=x^l}
- \underbrace{e_0 \left(
   \underline{n}_\text{Au}^{\#}
   -\underline{n}_\text{e} 
   + z_\text{V}\underline{n}_\text{V}
   +\underline{n}_\text{YSZ}^{\#}\left(
       \underline{m}_\text{C}^\text{YSZ} z_\text{C} - 2 \underline{m}^\text{YSZ}_\text{A}
   \right)
\right)}_{\text{ISR charge density}}
=
0
```
Analogously to the bulk, we define the *electroneutral* density of ISR vacancies as
```math
\underline{n}_\text{V}^\text{eq} 
=
\underline{S}_l (n_\text{YSZ}^{\#})^\frac{2}{3}\left(
  \underline{m}_\text{A} 
  - \underline{m}_\text{C} \frac{z_\text{C}(x^{\#})}{z_\text{V}} 
\right)
```
and again, we assume that
```math
0 
\leq \underline{n}_\text{V}^\text{eq}
\leq \underline{n}_\text{V}^\text{max}
\leq S_l\underline{m}_\text{A}(n_\text{YSZ}^{\#})^\frac{2}{3}~.
```
```math
\underline{n}_\text{V}^\text{max}(\alpha, x^{\#}) 
:= 
(1-\underline{\alpha}) \underline{n}_\text{V}^\text{eq}(x^{\#}) 
+ \underline{\alpha} \underline{m}_\text{A} \underline{n}_\text{YSZ}^{\#}
```
Hence the vacancy concentration now reads,
```math
\begin{aligned}
\text{ISR vacancy coverage}\quad \underline{n}_\text{V}(\underline{y}_\text{V}, \underline{\alpha}, x^{\#}, \underline{S}_l) 
&= \underline{n}_\text{V}^\text{max}(\underline{\alpha}, x^{\#},\underline{S}_l )y_\text{V} 
\\
&= \underline{y}_\text{V} ({n}_\text{YSZ}^{\#})^\frac{2}{3} \underline{S}_l\Big(
    \underline{m}_\text{A} 
     - (1-\underline{\alpha}) \underline{m}_\text{C} \frac{z_\text{C}(x^{\#})}{z_\text{V}}  
  \Big)~. 
\end{aligned}
```
where ``\underline{y}_\text{V}`` is the *ISR vacancy coverage*.
### Normal transport (adsorption) of oxide vacancies
The adsorption proceeds in positive direction from the YSZ bulk to the ISR,
```math
\text{V}^{..}_\text{O}  + \text{O}_\text{O}(s) 
\rightarrow 
\text{V}^{..}_\text{O}(s)  + \text{O}_\text{O} 
~.
```
This means, that rate 
```math
\underline{R}_\text{V} 
= 
\underline{k}_\text{V} 
\left(
    K_\text{V}^{\frac{1}{2}}  y_\text{V} (1 - \underline{y}_\text{V}) 
    - 
    K_\text{V}^{-\frac{1}{2}}  \underline{y}_\text{V}(1 - y_\text{V})
\right)
\quad\text{where}\quad
\underline{K}_\text{V}^{\frac{1}{2}}
= 
\exp
\left(
    -\frac{\Delta \underline{G}_\text{A}}{ 2 k_\text{B} T} 
    + \frac{
        A_\text{YSZ} y_\text{V} 
        - \underline{A}_\text{YSZ} \underline{y}_\text{V}
    }{2}
\right)
```
enters the LHS of the ISR balance with negative sign. 
```math
\partial_t \underline{y}_\text{V}
-
\frac{1}{\underline{S}_l \underline{n}_\text{V}^\text{max}}\underline{R}_\text{V}
=0
```
Finally, the rate has to be equal to the outflux of bulk vacancies, which induces the following Robin boundary condition,
```math
-\vec{\mathbf{j}}_\text{V}\cdot\nu_\text{YSZ} 
+ 
\frac{1}{n_\text{V}^\text{max}}\underline{R}_\text{V} 
= 0
```
### Drift-diffusion equilibrium of ISR electrons
```math
\underline{n}_\text{e} 
=
\frac{\underline{n}_\text{Au}^{\#}}{{n}_\text{Au}^{\#}}
\exp\left(
    -\frac{\Delta \underline{G}_\text{e}}{k_\text{B} T}
\right)
\exp{\frac{z_\text{e} e_0}{k_\text{B}T}(\psi_{l,r} - \psi(x))}
```

## Sedan scheme 
We discretize the drift-diffusion fluxes using the Sedan scheme.
Roughly, the main idea of the Sedan scheme splits the normal projection of the flux to the diffusion part and the rest:
```math
\vec{j} {\cdot}  = - ( \nabla u\cdot  + uq), u(x_1) = u_1, u(x_2) = u_2
```

dirichlet BCs u1        u2 
line segment  |---------|
colloc. pts  x1        x2
h = |x2-x1|
g = q*h

### Lattice-gas flux 

### Nernst-Planck flux
