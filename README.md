# Simulation of the Penetration of thermal Neutrons through material

## Abstract
In this report, Monte Carlo methods were investigated and then used to build a simulation of thermal neutrons moving through different types of shielding via recording the particle’s history. The simulation was then used to calculate the characteristic attenuation of water, lead and graphite and found to be λ_c= 1.950±0.001,9.34±0.01 and 11.22±0.02. However, the model used to obtain them could not describe all the behaviour in the survivability of neutrons throughout the material. Further complexity in the simulation was built by using multiply materials; this demonstrated the power and versatility of the method.

## 1. Introduction
The modern Monte Carlo method was first developed by Stanislaw Ulam when working on the Manhattan project to solve neutron diffusion in fissionable material. All the factors were known and understood, but the problem was too complex to solve by analytical methods. The name “Monte Carlo” was suggested as a code name by one of the team members after a casino in Monaco [1]. 
A similar technique was used in 18th century by Georges-Louis Leclerc and Comte de Buffon to estimate π, by scattering needless along an evenly striped floor and working out the probability that the needles stand between the strips [2]. This method is extremely powerful across a wide arrangement of problems and is used frequency in modern-day physics, so it is vital to understand the basics.

## 2. Theory of Neutron scattering
This report will look at estimating the penetration of neutrons through different types of shielding by simulating the random motion of each neutron through the material. When a neutron enters the shielding, it can fly through or collide with an atom and so will travel an average distance through the material before it collides. The average number of particles in a given volume, n, is,

### <img src="https://latex.codecogs.com/gif.latex?n=&space;\frac{\rho&space;N_A}{M_A}" title="n= \frac{\rho N_A}{M_A}" />  (1)

where ρ is the mass density, M_A is the atomic mass and N_A=6.02214076×10^23 Avogadro constant. Now assume that each atom has some cross-section in which a collision will occur, σ, and that the nucleus of the material is much greater than a neutron, σ≫σ_n. The probability of no collision occurring in some distance x+dx, is,

### <a href="https://www.codecogs.com/eqnedit.php?latex=P(x&plus;dx)=P(x)[1-P(\overline{dx})]," target="_blank"><img src="https://latex.codecogs.com/gif.latex?P(x&plus;dx)=P(x)[1-P(\overline{dx})]," title="P(x+dx)=P(x)[1-P(\overline{dx})]," /></a>

where P(x) is the probability of no collision occurring in a distance of x and P((dx) ̅) is the probability of a collision in distance dx. If the neutron is travelling through dx a collision will occur if a nucleus and neutron occupy the same cylinder shown in fig.1a. Therefore, P((dx) ̅ )=nσdx, making eq.2 become,

### <img src="https://latex.codecogs.com/gif.latex?P(x&plus;dx)=P(x)[1-n\sigma&space;dx]" title="P(x+dx)=P(x)[1-n\sigma dx]" /> (2)

Taylor expanding P(x+dx) around x, ignoring O(〖dx〗^2 ) terms, and cancelling the P(x) and dx gives,

### <img src="https://latex.codecogs.com/gif.latex?\frac{\mathrm{d}&space;P(x)}{\mathrm{d}&space;x}=-P(x)n&space;\sigma" title="\frac{\mathrm{d} P(x)}{\mathrm{d} x}=-P(x)n \sigma" /> (3)

which can be solved and normalized, removing the constant of integration, forming,

### <img src="https://latex.codecogs.com/gif.latex?P(x)=e^{\frac{-x&space;}{\lambda}}" title="P(x)=e^{\frac{-x }{\lambda}}" /> (4)

where λ=1/nσ and λ is the mean free path [4]. The area σ has no physical or geometric significance due to the collisions being purely quantum mechanical interactions, so σ is determined experimentally. It is also noted that σ is not a constant and will depend on many factors such as the energy of the neutron and the temperature of the shielding [5]. This report will only look at thermal neutrons; neutrons that have similar energy to the internal energy of the material. There are also fast neutrons, which have much higher energy than the material they move into and are important to consider in more complex systems.

### <img>

Once a collision occurs many different processes can occur such as fission, elastic scattering, inelastic scattering, absorption and even more; each of these interactions has a corresponding area [5]. This report will keep things simple and only look at elastic scattering and absorption, therefore the total collision area is,

### <img src="https://latex.codecogs.com/gif.latex?\sigma_{total}=\sigma_{scattering}&plus;\sigma_{absorbtion}" title="\sigma_{total}=\sigma_{scattering}+\sigma_{absorbtion}" /> (5)

There is some probability that one process occurs when colliding, for example, absorption occurring if the collision is more direct. This can be thought of that concentric circles of areas each corresponding to a different process, as shown in fig.1b, so the probability that one process occurs is,

### <img src="https://latex.codecogs.com/gif.latex?P(proccess)=\frac{\sigma_{proccess}}{\sigma_{total}}" title="P(proccess)=\frac{\sigma_{proccess}}{\sigma_{total}}" /> (6)

If scattering occurs the neutron, it will be assumed they will travel off in a uniformly random direction, which is a safe assumption to make with thermal neutrons. So, it now becomes clear that even from all these assumptions calculating the characteristic length a neutron will travel is a complex problem.

## 3. Simulation approach

The basis of this simulation is to take a statistical approach, so a way of producing random numbers is required. Linear-congruential generators have been historically used, which operate by using the simple sequence, 

