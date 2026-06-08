# Interplanetary Trajectory Planning
This repo contains the code for handling interplanetary trajectory planning for rendez-vous mission in the solar system.
The scenarios are generated randomly from the Solar System and we use the [`pykep`](https://github.com/esa/pykep/) package from the ESA to compute trajectories' parameters. This package relies on data from the Jet Propulsion Laboratory of NASA which can also be parsed and utilized using the `jplephem` package. Description of the package and the data source are detailed [here](https://pypi.org/project/jplephem/).

The scripts are implementations of **Monte Carlo Tree Search** based algorithms.

We also use data from the Global Trajectory Optimization Problem [portal](https://sophia.estec.esa.int/gtoc_portal/) previous competitions.

## Multi-Gravity Assist model (MGA)

We use the MGA model described to calculate velocities of the spacecraft by solving the Lambert equation [1]. The code for this model was generated with the help of Claude from Anthropic.

The sequential decision making nature of this problem explains the use of MCTS algorithms, and the decision variables are the following:
* The planets' sequence, starting from Earth and ending at the target planet. This can later be modified to included other trajectories (as for the scenarios from GTOC).
* The launch parameters from the start planet, defined by the launch window and a vector (velocity and angles)
* Times of flight between planets
* Arrival parameters at each planet: angle and radius.

In general, the goal is to minimize the variation of velocity through the voyage and respect the feasability constraints.


## References
[1] Song, Yuqi, et al. "Gravity assist space pruning and global optimization of spacecraft trajectories for solar system boundary exploration." Complex & Intelligent Systems 10.1 (2024): 323-341. [Link](https://link.springer.com/content/pdf/10.1007/s40747-023-01123-2.pdf)