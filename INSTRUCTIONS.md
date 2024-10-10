# Modelling mechanobiological processes

## Origin

This tutorial was authored by Claire Pritchard, Guillaume Stirnemann, and Glen M. Hocky in October 2024

## Aims

This Masterclass explains how mechanical forces can be modeled using PLUMED, and how to compute their result using free energy and kinetics calculations

## Objectives

The objectives of this Masterclass are:
- Learn how to apply a constant force in PLUMED
- See how transition rates change with different forces

## Prerequisites

We assume that the person that will follow this tutorial is familiar with the Linux terminal, Gromacs and basic functionality of PLUMED.

Familiarity with python and matplotlib is recommended for plotting

## Setting up PLUMED

We will use GROMACS, PLUMED, and PLUMED's pesmd function to perform the calculations.
Conda packages with the software required for this class have been prepared and you can install them following the instructions in [this link](https://github.com/plumed/masterclass-2022).

The data needed to run the exercises of this Masterclass can be found on [GitHub](https://github.com/hockyg/plumed-tutorial-force2).
You can clone this repository locally on your machine using the following command:

````
git clone https://github.com/hockyg/plumed-tutorial-force2
````

## Background

A force along some direction here is defined as the negative gradient of the potential energy along that direction. 

A constant force $F$ on a scalar collective variable $Q(\vec{X})$ therefore is a simple addition to the system's energy function.

$$
 U(\vec{X},F) = U(\vec{X}) - F Q(\vec{X}) 
$$

Notice that, because of the negative sign, a postive value of $F$ results in a lower energy for large values of $Q$, meaning $F\gt0$ corresponds to a "pulling" force.

A mechanical force would often in reality would correspond to pulling apart two atoms, and so $Q$ would often be a simple distance coordinate.

Note however, that other quantities could be used, such as an area which would mean $F$ corresponds to a pressure.

Dimensional analysis implies that the units of $F$ must be [Energy]/[Q].

The effect of constant force can be assessed using any enhanced sampling method.

The work done by the bias in this case is very simple

$$
W = \int_a^b F \cdot d Q = F ( Q_b - Q_a )
$$

Constant forces can be applied in PLUMED with the SLOPE keyword of the RESTRAINT bias.

## Steered MD

Steered molecular dynamics (SMD) is one way of pulling on a molecular coordinate, and has a connection to experiments done where a molecule is attached via a "spring" to an object such as an optical tweezer or AFM tip. To represent this in simulation, instead of applying a constant force, we impose a Harmonic restraint on $Q$ with a center that moves:

$$
 U(\vec{X},F) = U(\vec{X}) + \frac{1}{2} k (Q-Q_0(t))^2
$$

Typically $Q_0(t)$ would move linearly, with $ Q_0(t) = Q_0(0)-\lambda t $ although that is not a requirement. 

At any given time, the force along $Q$ from the moving bias is given as:

$$
 F(t) = -\frac{\partial U}{\partial Q} = -k(Q-Q_0(t))
$$

This force is positive (pulling) when $Q_0(t)$ is bigger than $Q$, and it can get very large if the spring moves quickly to larger values.

SMD is implemented in PLUMED using the MOVINGRESTRAINT bias, and the work is computed automatically.

$$
W = \int_a^B F dQ \approx \sum_{i=1}^{N_{steps}} \bar{F}_i (Q_0(t_i)-Q_0(t_{i-1})) = \lambda dt \sum_{i=1}^{N_{steps}} \bar{F}_i, 
$$

$$
\bar{F_i}=\frac{1}{2}( F_{i}+F_{i-1} ) = -\frac{k}{2}( Q(t_i)-Q_0(t_i) -Q(t_{i-1})+Q_0(t_{i-1})) = -\frac{k}{2} (\Delta Q_i - \lambda d t)
$$

## Exercises

The exercises are presented below.

### Effect of force on a 2-dimensional potential

Use the RESTRAINT function to add a constant force of different magnitudes (e.g. -2.5 to 2.5 in these units) and look at how the force changes the resulting free energy surface.

```plumed
#SOLUTIONFILE=work/plumed_ex1.dat
UNITS ENERGY=kcal/mol

d1: DISTANCE ATOMS=1,2
ff: MATHEVAL ARG=d1 PERIODIC=NO FUNC=0.2*(((x-10)^2)*((x-20)^2))
metad: METAD ARG=d1 PACE=500 HEIGHT=0.1 SIGMA=2.5 FILE=__FILL__ BIASFACTOR=10 TEMP=300.0 GRID_WFILE=__FILL__ GRID_MIN=0 GRID_MAX=30 GRID_BIN=251 GRID_WSTRIDE=10000

RESTRAINT __FILL__

PRINT ARG=* FILE=__FILL__ STRIDE=100
```

Then run the simulation using the command:

````
plumed pesmd < doublewell_prod.pesmd.input
````

Plot the free energy surface from the GRID or after using sum_hills to compute the surface, and zero the potential at the left minimum. What do you notice about the other minimum and barrier?

![Pulling on a double well, sampled by metadynamics](figs/masterclass-22-15-doublewell_metad.jpg)

### Ligand dissociation from a membrane protein system 

#### System and context 

In the following, we illustrate how Plumed can be used to probe ligand-protein interactions in a constant force set-up mimicking single-molecule force spectroscopy experiments that can produce such forces from the pN to the nN regime. On purpuse, we apply it here to a large system comprising the ghrelin receptor GHSR bound to a model peptide GHRP-6 embedded in a model infinite lipid bilayer and solvated by water molecules. All input files necessary to run these simulations are available here XX. Given the system size, and required simulation lengths, it is highly recommended that you perform these simulations on a cluster, preferably on GPUs. Note that these simulations can be run on a laptop or desktop computer, but you will be limited to quite short simulation times which would require the use of much larger forces to see dissociation within a short time window. 

We start by giving a little bit of context on the system chosen for this example, which is shown in the figure below: 

![tutorial_figure](https://github.com/user-attachments/assets/4f587ac5-34a1-4abc-9afc-3669a01c96ad)

 Panel a gives an overview of the system {GHSR+ligand+membrane+water} which is periodic in all directions. In panel b and c, we show the geometric definition of two collective variables, $d_S$ and $d_O$, corresponding respectively to the distance between an anchor in the protein (center of mass between three residues), and either the N-terminal alpha carbon ($d_S$) or the center of mass of the ligand ($d_O$). In plumed, this is achieved as follows: 

```plumed
GHRP6_COM: COM ATOMS=4822-4942
GHSR_COM: COM ATOMS=769,1490,3760
dO: DISTANCE ATOMS=GHRP6_COM,GHSR_COM COMPONENTS NOPBC
dS: DISTANCE ATOMS=4826,GHSR_COM COMPONENTS NOPBC
```
Please note the NOPBC option that is important here as we don't want to calculate the minimal distance between these two ligand/protein sites. 

#### Applying a constant, vertical force 

Say now that we want to apply a vertical force of a given amplitude to one of these CVs. As detailed in the preambule, this can be achieved by playing with the SLOPE keyword of the RESTRAINT bias. The fact that we apply a vertical force is controlled by the fact that bias is acting on the z component of the CV: 

For example, 
```plumed
restraint: RESTRAINT ARG=dS.z AT=0 SLOPE=-270
```

Since the energy and distance units are respectively XX and XX, can you guess to which force this corresponds? Why is the SLOPE defined with a negative sign?

The amplitude of the force can be tuned by changing the SLOPE, and this can be applied to any CV you like. Steered and constant force MD are by essence out-of-equilibrium set-ups, therefore, if the force is high enough, you don't expect the unbinding event to be reversible. This means that to calculate converged quantities (e.g. unbinding times), you need to repeat such simulations multiple times. 

For example, the following figures shows the results (CV vs time) of 5 trajectories at two different forces acting on $d_S$. How does the unbinding time depends on force? Is this expected?

![different_forces](https://github.com/user-attachments/assets/e458c453-1ae3-4192-95ca-21048878b3c4)

As an exercise, you can either try to play with force further and re-run some simulations to see how this affects the unbinding time, or try to apply force to a different CV, $d_O$ for example. 

#### Pulling at an angle

Finally, it could be interesting and relevant for the comparision with experimental results to change the directionnality of the force. For example, providing you with the following plumed file, and based on what we did above, can you guess what is done exactly to the system? Is is expected to have a significant influnce on the unbinding kinetics?

For example, 
```plumed
GHSR_COM: COM ATOMS=769,1490,3760
dS: DISTANCE ATOMS=4826,GHSR_COM COMPONENTS NOPBC
restraint_x: RESTRAINT ARG=dS.x AT=0 SLOPE=-173
restraint_z: RESTRAINT ARG=dS.z AT=0 SLOPE=-207
PRINT FMT=%g ARG=dS.*,restraint_x.*,restraint_z.* STRIDE=500 FILE=colvars.data

```

In the following figure, we show the results of the application of a force of 450 pN to the N-terminal residue of the ligand peptide, pulling at different angles. The first angle is that corresponds to the tilt along x, and the second one, if shown, to that along y. 

![different_angles](https://github.com/user-attachments/assets/e19dce85-7a76-407a-9af7-cf17e2cfbf0b)

#### Conclusion

You should now be able to apply constant forces to any CV in a given system, and to tune the direction of this force. 
