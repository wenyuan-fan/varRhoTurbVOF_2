# varRhoTurbVOF

varRhoTurbVOF contains a set of OpenFOAM volume of fluid (VOF) solvers for turbulent isothermal multiphase flows, which are variable-density incompressible. Unlike their official counterparts, where Favre-averaged and Reynolds-averaged velocities coexist in different equations, new solvers use Favre-averaged velocities consistently in all equations. This major update introduces three main improvements to the previous version of varRhoTurbVOF. First, the implementation is extended to VOF solvers for isothermal and non-isothermal phase change two-phase flows, where the flow is no longer incompressible. Second, in order to introduce backward compatibility and to avoid code duplication, the turbulence model construction procedure is redesigned such that solvers can determine whether the variable-density effect is considered or not in the turbulence modeling part based on the input file at run time. Third, the Egorov turbulence damping model for \omega-based turbulence models is implemented with its most recent developments. Plus, an extension to \epsilon-based turbulence models is developed and implemented.

Reference:

Wenyuan Fan and Henryk Anglart. "varRhoTurbVOF: A new set of volume of fluid solvers for turbulent isothermal multiphase flows in OpenFOAM." Computer Physics Communications (2020), 247, 106876

Wenyuan Fan and Henryk Anglart. "varRhoTurbVOF 2: Modified OpenFOAM volume of fluid solvers with advanced turbulence modeling capability." Computer Physics Communications (2020), 256, 107467.

## File structure
  * OpenFOAM-7: folder for code to be compiled with OpenFOAM v7;
  * OpenFOAM-8: folder for code to be compiled with OpenFOAM v8;
  * OpenFOAM-1912: folder for code to be compiled with OpenFOAM v1912;
  * OpenFOAM-2006: folder for code to be compiled with OpenFOAM v2006;
  * bruteForceExamples: folder for three examples of the brute-force approach;
  * turbulenceDamping: folder for the turbulence damping fvOptions;
  * tutorials: folder for tutorials;
  * Manual.pdf: a manual focusing on how to use new features in varRhoTurbVOF 2;
  * LICENSE: the license file.


## Installation
The installation of a supported version of OpenFOAM is a prerequisite for using the newly designed solvers and the turbulence damping fvOptions.
In order to use the code, one needs to enter the desired version folder and load the corresponding environment variable for OpenFOAM.
Then run "./Allwmake" to compile the code.

## Using solvers
As for the usage for a specific solver, e.g. varRhoInterFoam, it is almost the same with the corresponding existing solver interFoam. For instance, it could be executed on 1024 processes by simply typing "mpirun -np 1024 varRhoInterFoam -parallel" in the terminal. When the full-form turbulence models are used in the new solvers, the corresponding discretization schemes should be provided to solve the governing equations numerically. Other than this, the users could reuse all their interFoam input files for varRhoInterFoam.

## Running tutorials
Several tutorials are provided in the "tutorials" folder, the "Allrun" file in each tutorial is used to run the simulation.
