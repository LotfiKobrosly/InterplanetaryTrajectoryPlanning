# Interplanetary Trajectory Planning
This repo contains the code for handling interplanetary trajectory planning for rendez-vous mission in the solar system.
The scenarios are generated using data from the Jet Propulsion Laboratory of NASA and handled using the `jplephem` package. Description of the package and the data source are detailed [here](https://pypi.org/project/jplephem/).

The scripts are implementations of **Monte Carlo Tree Search** based algorithms.

We also use data from the Global Trajectory Optimization Problem [portal](https://sophia.estec.esa.int/gtoc_portal/) previous competitions.

## Multi-Gravity Assist model (MGA)

We use the MGA model described to calculate velocities of the spacecraft by solving the Lambert equation [1]. The code for this model was generated with the help of a LLM.






## References
[1] Song, Yuqi, et al. "Gravity assist space pruning and global optimization of spacecraft trajectories for solar system boundary exploration." Complex & Intelligent Systems 10.1 (2024): 323-341. [Link](https://link.springer.com/content/pdf/10.1007/s40747-023-01123-2.pdf)