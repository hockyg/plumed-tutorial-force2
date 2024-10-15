# Modelling mechanobiological processes

## Origin

This tutorial was authored by Claire Pritchard, Guillaume Stirnemann, and Glen M. Hocky in October 2024.

## Aims

This Masterclass explains how mechanical forces can be modeled using PLUMED, and how to compute their result using free energy and kinetics calculations.

## Objectives

The objectives of this Masterclass are:
- Learn how to apply a constant force in PLUMED;
- See how transition rates change with different forces.

## Prerequisites

We assume that the person that will follow this tutorial is familiar with the Linux terminal, Gromacs and basic functionality of PLUMED.

Familiarity with python and matplotlib is recommended for plotting

## Setting up PLUMED

We will use GROMACS, PLUMED, and PLUMED's pesmd function to perform the calculations.
Conda packages with the software required for this class have been prepared and you can install them following the instructions in [this link](https://github.com/plumed/masterclass-2022).

The data needed to run the exercises of this Masterclass can be found on [GitHub](https://github.com/hockyg/plumed-tutorial-force2).

You can clone this repository locally on your machine if you want, but specific data will be linked to below.

## Background

### Force basics

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

The work done by the bias in this case is very simple:

$$
W = \int_a^b F \cdot d Q = F ( Q_b - Q_a )
$$

Constant forces can be applied in PLUMED with the SLOPE keyword of the RESTRAINT bias.

### Steered MD

Steered molecular dynamics (SMD) is one way of pulling on a molecular coordinate, and has a connection to experiments done where a molecule is attached via a "spring" to an object such as an optical tweezer or AFM tip. To represent this in simulation, instead of applying a constant force, we impose a Harmonic restraint on $Q$ with a center that moves:

$$
 U(\vec{X},F) = U(\vec{X}) + \frac{1}{2} k (Q-Q_0(t))^2
$$

Typically $Q_0(t)$ would move linearly, with $Q_0(t) = Q_0(0)-\lambda\ t$ although that is not a requirement. 

At any given time, the force along $Q$ from the moving bias is given as:

$$
 F(t) = -\frac{\partial U}{\partial Q} = -k(Q-Q_0(t))
$$

This force is positive (pulling) when $Q_0(t)$ is bigger than $Q$, and it can get very large if the spring moves quickly to larger values.

SMD is implemented in PLUMED using the MOVINGRESTRAINT bias, and the work is computed automatically:

$$
 W = \int_a^b F dQ \approx \sum_{i=1}^{N_{steps}}  \bar{F_i} (Q_0(t_i)-Q_0(t_{i-1})) = \lambda dt \sum_{i=1}^{N_{steps}} \bar{F_i} 
$$

where,

$$
 \bar{F_i}=\frac{1}{2}( F_{i}+F_{i-1} ) = -\frac{k}{2}( Q(t_i)-Q_0(t_i) -Q(t_{i-1})+Q_0(t_{i-1})) = -\frac{k}{2} (\Delta Q_i - \lambda d t)
$$

### Dynamics from Infrequent Metadynamics

To compute the rate constant of unbinding in a standard molecular dynamics simulation, one would run many simulations starting from the bound pose, and then compute the mean time of escape, by averaging the times at which unbinding is deemed to have occured, $\tau=\frac{1}{N} \sum_{i=1}^N t^{unbind}_i$ (this assumes every simulation ends in an unbinding event). 

In this case, the unbinding rate constant is $k_{off}=1/\tau$. 


One way to compute dynamics from biased simulations is [Infrequent Metadynamics](https://doi.org/10.1103/PhysRevLett.111.230602). In this approach, we apply bias with a slow PACE so that the transition over the barrier to escape is not likely to be biased. In this case with some other assumptions, $\tau = \frac{1}{N} \sum_{i=1}^{N} \alpha_i \t_i$, where 

$$
\alpha_i = \frac{1}{t_i} \int_0^{t_i} e^{\beta V_i(t)} dt
$$

Where $V_i(t)$ is the Metadynamics bias at time $t$.

$\alpha_i$ is called the _acceleration factor_ and can be computed by PLUMED for METAD automatically.

Note - there are other approaches to take the same data and compute $\tau$ by fitting the survival distribution, but we will not discuss that here. 

### Dependence of unbinding rate constants with force

For a double well potential, the rate of crossing a high barrier with a small pulling force follows what is sometimes called Bell's law, and is a simple consequence of the linear response of the barrier to pulling:

$$
k_{off}(F)=k_{off}(F=0) e^{\beta F \Delta Q^*}
$$

This means that the unbinding accelerates exponentially with force relative to zero force, and this speedup depends on $\Delta Q^*$ which is the "distance to the transition state", assumed to be constant (in reality, this distance moves even for a double well potential, and it is easy to correct for this motion for a quadratic barrier. Feel free to tackle this as an exercise!).

## Exercises

The exercises are presented below.


### Part 1 - Ligand dissociation from a membrane protein system 

#### System and context 

In the following, we illustrate how Plumed can be used to probe ligand-protein interactions in a constant force set-up mimicking single-molecule force spectroscopy experiments that can produce such forces from the pN to the nN regime. 

On purpose, we apply it here to a large system comprising the ghrelin receptor GHSR bound to a model peptide GHRP-6 embedded in a model infinite lipid bilayer and solvated by water molecules. 

All input files necessary to run these simulations are available [in this zip file](https://github.com/hockyg/plumed-tutorial-force2/raw/refs/heads/main/plumed-tutorial-force2-data-part1.zip).

Given the system size, and required simulation lengths, it is highly recommended that you perform these simulations on a cluster, preferably on GPUs. Note that these simulations can be run on a laptop or desktop computer, but you will be limited to quite short simulation times which would require the use of much larger forces to see dissociation within a short time window. 

We start by giving a little bit of context on the system chosen for this example, which is shown in the figure below: 

![tutorial_figure](images/GHRP-setup.png).

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
#SOLUTIONFILE=work/pulling_units.plumed.dat
UNITS __FILL__

GHRP6_COM: COM ATOMS=4822-4942
GHSR_COM: COM ATOMS=769,1490,3760
dO: DISTANCE ATOMS=GHRP6_COM,GHSR_COM COMPONENTS NOPBC
dS: DISTANCE ATOMS=4826,GHSR_COM COMPONENTS NOPBC

restraint: RESTRAINT ARG=dS.z AT=0 SLOPE=-270

PRINT __FILL__ STRIDE=100 FILE=pull_270.colvar.dat
```

Set the energy and distance units to kJ/mol and nm respectively. In these units, what does a SLOPE of -270 correspond to in Piconewtons? (And remember from above,  why is the SLOPE defined with a negative sign?)

Make sure to print out the distances you are interested in, as well as the work done by the bias.

The amplitude of the force can be tuned by changing the SLOPE, and this can be applied to any CV you like. Steered and constant force MD are by essence out-of-equilibrium set-ups, therefore, if the force is high enough, you don't expect the unbinding event to be reversible. This means that to calculate converged quantities (e.g. unbinding times), you need to repeat such simulations multiple times. 

For example, the following figures shows the results (CV vs time) of 5 trajectories at two different forces acting on $d_S$. How does the unbinding time depend on force? Is this expected?

![different_forces](images/GHRP_distance_time.png).

As an exercise, set slope to match these different conditions. can either try to play with force further and re-run some simulations to see how this affects the unbinding time, or try to apply force to a different CV, $d_O$ for example. 

#### Pulling at an angle

Finally, it could be interesting and relevant for the comparison with experimental results to change the directionality of the force. 

For example, providing you with the following plumed file, and based on what we did above, can you guess what is done exactly to the system? Is it expected to have a significant influnce on the unbinding kinetics?

For example, 
```plumed
#SOLUTIONFILE=work/pull_angle.plumed.dat
UNITS __FILL__

GHSR_COM: COM ATOMS=769,1490,3760
dS: DISTANCE ATOMS=4826,GHSR_COM COMPONENTS NOPBC
restraint_x: RESTRAINT ARG=dS.x AT=0 SLOPE=__FILL__
restraint_z: RESTRAINT ARG=dS.z AT=0 SLOPE=__FILL__

PRINT FMT=%g ARG=dS.*,restraint_x.*,restraint_z.* STRIDE=500 FILE=pull_40degrees.dat
```

Fill in this PLUMED file so that constant force pulling is performed with a total magnitude of 450 pN at an angle of 40 degrees from the vertical.

In the following figure, we show the results of the application of a force of 450 pN to the N-terminal residue of the ligand peptide, pulling at different angles. The first angle is that corresponds to the tilt along x, and the second one, if shown, to that along y. 

![different_angles](images/GHRP_pulling_angles.png).

#### Conclusion

You should now be able to apply constant forces to any CV in a given system, and to tune the direction of this force. 

### Part 2 - Effect of force on unbinding rates

In this section, we will demonstrate how unbinding rate constants can be computed via Infrequent Metadynamics. This will be demonstrated on a double-well potential, and simulations will be performed using plumed's `pesmd` function.

Note, the first step is a recap from the earlier [FISST module masterclass/tutorial](https://www.plumed-tutorials.org/lessons/22/015/data/INSTRUCTIONS.html). 

#### Reminder - Effect of force on a 2-dimensional potential

Use the RESTRAINT function to add a constant force of different magnitudes (e.g. -2.5 to 2.5 in these units) and look at how the force changes the resulting free energy surface.

```plumed
#SOLUTIONFILE=work/pulling_example_doublewell.plumed.dat
UNITS ENERGY=kcal/mol

d1: DISTANCE ATOMS=1,2
ff: MATHEVAL ARG=d1 PERIODIC=NO FUNC=0.2*(((x-10)^2)*((x-20)^2))

metad: METAD ARG=d1 PACE=500 HEIGHT=0.1 SIGMA=2.5 FILE=__FILL__ BIASFACTOR=10 TEMP=300.0 GRID_WFILE=__FILL__ GRID_MIN=0 GRID_MAX=30 GRID_BIN=251 GRID_WSTRIDE=10000

RESTRAINT __FILL__

PRINT ARG=* FILE=__FILL__ STRIDE=100
```

Then download the [pesmd input file](inputs/doublewell_prod.pesmd.input) run the simulation using the command:

````
plumed pesmd < doublewell_prod.pesmd.input
````

Plot the free energy surface from the GRID or after using sum_hills to compute the surface, and zero the potential at the left minimum. What do you notice about the other minimum and barrier?

![Pulling on a double well, sampled by metadynamics](images/doublewell_force_metad.jpg)

### Getting transition times versus force for a double well

Perform MetaD using the options below to predict kinetics at different forces. In the figure, we show results from computing the unbinding time using 40 runs with independent seeds for each _pulling_ force from 0 to 1 (be careful with the sign). As you can see, the time goes down exponentially, following Bell's law. Can you extract the transition distance from this figure? 

Note that here we used reduced units. An example pesmd input file is given at [this link](inputs/doublewell_rates.pesmd.input). 

Fill in the blanks below to turn on calculation of the acceleration factor, and use the COMMITTOR function to stop simulations when the system crosses the barrier (e.g. at x=15.5)

```plumed
#SOLUTIONFILE=work/pulling_rate_example.plumed.dat
UNITS NATURAL

d1: DISTANCE ATOMS=1,2
ff: MATHEVAL ARG=d1 PERIODIC=NO FUNC=0.02*(((x-10)^2)*((x-20)^2))
bb: BIASVALUE ARG=ff

metad: METAD ARG=d1 PACE=200 HEIGHT=0.1 SIGMA=2.5 FILE=run_metad_F-1.0.hills BIASFACTOR=10 TEMP=1.0 __FILL__

COMMITTOR ...
    ARG=d1
    STRIDE=10
    BASIN_LL1=__FILL__
    BASIN_UL1=10000
    LABEL=c
...

F: RESTRAINT ARG=d1 AT=0.0 KAPPA=0 SLOPE=__FILL__

PRINT ARG=* FILE=run_metad_F-1.0.metad.dat STRIDE=10
```

Your results should look like the following:

![Taus from double well, from infrequent metadynamics](images/potential_and_taus.png)
