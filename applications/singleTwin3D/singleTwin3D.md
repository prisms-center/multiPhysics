## singleTwin3D: Evolution of a single $\{10\bar{1}2\}$ twin within a single crystal

We employ a Phase-Field / Crystal Plasticity models coupling strategy similar to the one introduced by by Li et al. in Ref. [1]. 

## Phase Field Model

The total free energy of the system can be expressed as

$$
\mathcal{F} = \int_{\Omega} \left( f_{tw} + f_{grad}+ f_{el} \right)~dV, 
$$

where $f_{tw}$ is a bulk energy density contribution defined as a double-well function with minima at values of the order parameter corresponding to the parent grain ($\eta=0$) and twin ($\eta=1$) phases,

$$
f_{tw}=\Delta f_{tw}~\phi^2(1-\phi)^2.
$$

The constant is $\Delta f_{tw}$ the energy barrier height between the parent crystal phase and the twin. The term $f_{grad}$ is an anisotropic contribution from the gradient of $\eta=1$ at twin boundaries, defined as

$$
f_{grad}=\frac{1}{2} \nabla\phi\cdot\boldsymbol{\kappa}\cdot\nabla\phi,
$$

where $\boldsymbol{\kappa}$ is a second-order anisotropic tensor proportional to the direction and magnitude of the twin boundary energy. Finally, $f_{el}$ is an elastic energy density contribution expressed as

$$
f_{el}=\frac{1}{2}\mathbf{E}:\mathbb{C}:\mathbf{E}.
$$

### Kinetics
The evolution of the twin boundary is described by Allen-Cahn dynamics:

$$
\begin{align}
\frac{\partial \phi}{\partial t} = -\mathbf{M}\frac{\partial \mathcal{F}}{\partial \phi}
\end{align}
$$

where $\mathbf{M}$ is an anisotropic mobility tensor mobility.

==Insert derivation of $\partial \mathcal{F} / \partial \phi$==

### Time discretization

Considering forward Euler explicit time stepping, we have the time discretized kinetics equation:

==Insert derivation==

### Weak formulation

In the weak formulation, considering an arbitrary variation $w$, the above equation can be expressed as a residual equation:

==Insert derivation==

The above values of $r_{\phi}$ and $r_{\phi x}$ are used to define the residuals in `applications/singleTwin3D/equations.cc`.

### CPFE model 

[1] G. Liu, H. Mo, J. Wang, and Y. Shen, Coupled crystal plasticity finite element-phase field model with kinetics-controlled twinning mechanism for hexagonal metals, Acta Mater. **202**, 399-416 (2021).
