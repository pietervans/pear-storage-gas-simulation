# Pear project (team 14)

# Description
At present, no good methods are available to measure internal gas concentrations in fruit.
Therefore, in recent years,a scientific computing approach has been adopted to simulate
and predict internal gas concentrations/distributions.
Furthermore, this approach allows to study the effect of fruit geometry (shape and size) or controlled storage conditions on local oxygen and carbon dioxide concentrations, while reducing experimental
costs.

The goal of this project is to simulate the oxygen and carbon dioxide concentrations inside a pear. A respiration–diffusion system for metabolic gas exchange is solved using the Finite Element Method. We refer to [this file](pear_assignment.pdf) for more information. For more information about the Finite Element method, see [References](#the-finite-element-method).

## Dependencies
- To run the code, you will need to have the C++ library [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) pre-installed.
- To use the provided python scripts such as [demo.py](demo.py), the following modules are required:
    - [numpy](https://pypi.org/project/numpy/)
    - [matplotlib](https://pypi.org/project/matplotlib/)

## Installation

- Make sure all project [dependencies](#dependencies) are installed.
- Clone the repository:
```
git clone --depth=1 https://gitlab.kuleuven.be/math-eng/h0t46a/2023/team14.git
```
- Compile the code:
```
cd team14 & chmod +x build & ./build
```

## Usage

To simulate the gas concentrations in the pear, simply run the demo.py file along with a mesh specification and a storage description. For example:

```
python3 demo.py uniform_1mm orchard
```

Results are stored in a file called c_solution.txt. The first half of values listed are the oxygen concentrations in the mesh vertices. The latter half represent the carbon dioxide concentrations.

The different meshes are stored in [this folder](matlab/pear_meshes/data/). To add a mesh, two files should be added. A points_{name}.txt file that specifies the coordinates of the vertices of your mesh. The triangle_{name}.txt file specifies where the triangles of your mesh lie by indexing the points. (Note: indexing starts at 1)

We have also included a script that shows how results improve for finer meshes:
```
python3 convergence_test.py orchard
```

## Storage descriptions

| Storage description | $T_{cel} (^\circ C)$ | $\eta_{O_2} (\%)$ | $\eta_{CO_2} (\%)$ |
| ------------------- | --------- | ------------ | ------------- |
| orchard             | 25        | 20.8         | 0.04          |
| shelf_life          | 20        | 20.8         | 0             |
| refrigerator        | 7         | 20.8         | 0             |
| precooling          | -1        | 20.8         | 0             |
| disorder_inducing   | -1        | 2            | 5             |
| optimal_ca          | -1        | 2            | 0.7           |

## Project status & Links
Our planning can be found in this [Google Sheet](https://docs.google.com/spreadsheets/d/1Fhpz5enRTB4le8ANImo4l1pByl8DFYfVXlU_I1k14jk/edit?usp=sharing).

[This Google document](https://docs.google.com/document/d/1FK4v7cVcI_l2MSTxBLiBNz6RzeOSgufXm_0JWumz0ao/edit?usp=sharing) contains some questions and remarks for internal use by the team.

The final presentation can be found [here](https://docs.google.com/presentation/d/1S1AQRmjl5wM7BUBt8en8uUPmlFcIurLW7UBM3sGOk1A/edit?usp=sharing).

## References

### The Finite Element Method
- Huebner, K.H., Thornton, E.A., Byrom, T.G., 1995. The finite element method for engineers.
John Wiley and Sons, Inc.
- [Nikishkov, G. P. Introduction to the Finite Element Method. 2004 Lecture Notes](http://nliebeaux.free.fr/ressources/introfem.pdf)
- [Gagandeep Singh, Short Introduction to Finite Element Method.](https://folk.idi.ntnu.no/elster/tdt24/tdt24-f09/gagan-fem.pdf)
- [Roache, Patrick J. Code verification by the method of manufactured solutions. J. Fluids Eng.
124.1 (2002): 4-10.](https://asmedigitalcollection.asme.org/fluidsengineering/article-abstract/124/1/4/462791/Code-Verification-by-the-Method-of-Manufactured)


