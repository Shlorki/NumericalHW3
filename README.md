<h3 align="center">Advanced Numerical Analysis HW 3</h3>

  <p align="center">
    Code for solving a periodic one dimensional linear hyperbolic problem over two spatial periods.
  </p>
</div>

<!-- GETTING STARTED -->
## Getting Started

The code requires Python to be installed on the machine.

### Prerequisites

The modules required in Python are
* tabulate
  ```sh
  pip install tabulate
  ```
* numpy
  ```sh
  pip install numpy
  ```
* scipy
  ```sh
  pip install scipy
  ```
* matplotlib
  ```sh
  pip install matplotlib
  ```

### Installation

1. Simply download HW3.py or HW3nb.pynb and run in your favorite Python environment
2. Install required modules

## Examples
We solve the baby wave equation $u_t=u_x$ with periodic boundary conditions over the interval $[0,1]$ and display the solution at time $t=0,0.5,1$. 

#### Backward Time and Centered Space (BTCS):

<a href="https://github.com/Shlorki/NumericalHW2">
  <img src="Images/Dirichlet.png" alt="helmpt" width="400" height="300">
  <img src="Images/DirichletCvg.png" width="150" height="300">
</a>

#### Crank_Nicholson Time and Centered Space (CNCS):

<a href="https://github.com/Shlorki/NumericalHW2">
  <img src="Images/Neumann-Dirichlet.png" alt="helmpt" width="400" height="300">
  <img src="Images/NeuDirCvg.png" width="150" height="300">
</a>

#### Lax-Friedrichs:

<a href="https://github.com/Shlorki/NumericalHW2">
  <img src="Images/Dirichlet-Robin.png" alt="helmpt" width="400" height="300">
  <img src="Images/DirRobCvg.png" width="150" height="300">
</a>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

#### Lax-Wendroff:

asdfsdf

#### RK4 in Time and Compact Differences in Space (RK4/CD4):

asdfsadf 
