# singleTwin3D

This application simulates the evolution of a single $\{10\bar{1}2\}$ twin within a single crystal.

We employ a Phase-Field / Crystal Plasticity models coupling strategy similar to the one introduced by by Li et al. in Ref. [1]. 

## Phase Field Model

The total free energy of the system can be expressed as

$$
\mathcal{F} = \int_{\Omega} \left( f_{tw} + f_{grad}+ f_{el} \right)~dV, 
$$

where $f_{tw}$ is a bulk energy density contribution defined as a double-well function with minima at values of the order parameter corresponding to the parent grain ($\phi=0$) and twin ($\phi=1$) phases,

$$
f_{tw}=\Delta f_{tw}~\phi^2(1-\phi)^2.
$$

The constant $\Delta f_{tw}$ is the energy barrier height between the parent crystal phase and the twin. The term $f_{grad}$ is an anisotropic contribution from the gradient of $\phi$ at twin boundaries, defined as

$$
f_{grad}=\frac{1}{2} \nabla\phi\cdot\pmb{\kappa}\cdot\nabla\phi,
$$

where $\pmb{\kappa}$ is a second-order anisotropic tensor proportional to the direction and magnitude of the twin boundary energy. Finally, $f_{el}$ is an elastic energy density contribution expressed as

$$
f_{el}=\frac{1}{2}\pmb{E^e}:\mathbb{C}:\pmb{E^e},
$$

where, $\mathbb{C}$, is the local elastic tensor, obtained via a linear interpolation between the value in the parent grain, $\mathbb{C}^0$, and value in the twin, $\mathbb{C}^{tw}$,

$$
\mathbb{C}=(1-\phi)\mathbb{C}^0 + \phi~\mathbb{C}^{tw}.
$$

and $\pmb{E^e}$ is the elastic Green-Lagrange strain tensor

### Kinetics
The evolution of the twin boundary is described by Allen-Cahn dynamics:

$$
\begin{align}
\frac{\partial \phi}{\partial t} = -M\frac{\delta \mathcal{F}}{\delta \phi}
\end{align}
$$

where $M$ is an anisotropic mobility that depends on the orientation of the twin boundary. The variational derivative $\delta \mathcal{F} / \delta \phi$ is given by

$$
\begin{align}
\frac{\delta \mathcal{F}}{\delta \phi} =\mu_\phi= 4 \Delta f_{tw}~\phi (\phi-1) ( \phi-0.5)
-\nabla\cdot(\pmb{\kappa}\nabla\phi) -\gamma_{tw}\tau_{tw},
\end{align}
$$

where $\gamma_{tw}$ is characteristic twin shear strain and $\tau_{tw}$ is the twinning resolved shear stress (TRSS) for the twin variant, obtained by solving the constitutive equations of the CPFE model.

### Time discretization

Considering forward Euler explicit time stepping, we have the discretized kinetics equation to aproximate the field $\phi$ at time increment $n+1$ from the fields $\phi$ and $\mu_\phi$  at increment $n$:

$$
\phi^{n+1}=\phi^{n}-M\Delta t~\mu_\phi^{n},
$$

where $\Delta t$ is the time increment. 

### Weak formulation

In the weak formulation, considering an arbitrary variation $\omega$, the above equation can be expressed as a residual equation:

$$
\begin{align}
\int_{\Omega} \omega \phi^{n+1} ~dV &= \int_{\Omega} \omega \left[ \phi^{n} - \Delta t M(4 \Delta f_{tw}~\phi^n (\phi^n-1) ( \phi^n-0.5) -\nabla\cdot(\pmb{\kappa}\nabla\phi^n) -\gamma_{tw}\tau^n_{tw})\right] ~dV \\
&=\int_{\Omega}\omega r_\phi~dV + \int_{\Omega}\nabla \omega\cdot \pmb{r}_{\phi x} ~dV,
\end{align}
$$

where

$$
\begin{align}
r_{\phi}= \phi^{n} - \Delta t M(4 \Delta f_{tw}~\phi^n (\phi^n-1) ( \phi^n-0.5) 
-\gamma_{tw}\tau^n_{tw})
\end{align}
$$

and

$$
\begin{align}
\pmb{r}_{\phi x} = -\Delta t M \pmb{\kappa}\nabla \phi^{n}
\end{align}
$$

are the right-hand side (RHS) scalar and gradient residual terms, respectively, that need to be computed in `applications/singleTwin3D/equations.cc`.

### CPFE model 

[1] G. Liu, H. Mo, J. Wang, and Y. Shen, Coupled crystal plasticity finite element-phase field model with kinetics-controlled twinning mechanism for hexagonal metals, Acta Mater. **202**, 399-416 (2021).
