# AvionicOptimizer

## 1 - Project Description
The goal of this optimizer is to find the best avionic architecture for a satellite dedicated to Active Debris Removal. It uses a MIQCQP ( Mixed Integer Quadratically-Constrained Quadratic Programming ) approach to implement the structure of the avionic. More details on the implementation can be found in part 3 of this readme.

## 2 - How to Use

This tool requires [Matlab](https://ch.mathworks.com/products/matlab.html) to run. In addition, it needs the [Gurobi Matlab libarary](https://www.gurobi.com/) that can be obtained for free with an academic license. The link to the library should be updated in the first line of the main.m file.

The "main.m" file contains the implementation of the MIQCQP into the library. It is also the one running the whole optimization. The "loadParamXML.m" load the various parameter of the architecture from the XML file stored in the "parameter" folder.

The "display" folder contains two standalone scripts to plot the data generated by the tool. 

The "parameter" folder contains the description of the avionic architecture that should be optimized. The new architecture can be added easily and should use the structure of the "param_test.xml" file.

The "display" folder contains the models of the system. New components can be added easily and should use the structure of the "testXXX.xml" file.

The "parametricAnalyses" folder has the Matlab code to load the parametric analyses. It is the same code as the main.m with a few additional linked to the parameter used. (In a future version, it will simply call the main.m file.)

The "results" folder will store the results of the optimizer.

## 3 - Research linked
Complement information the tool can be found here
- "Simulation Tool: Resources Management in High Performance Avionic for ADR Missions" by Michaël Juillard & Jean-Paul Kneib (2021 IEEE Aerospace Conference) [PDF] (https://www.researchgate.net/publication/352234234_Simulation_Tool_Resources_Management_in_High_Performance_Avionic_for_ADR_Missions)
- "Optimal Control Approach for Dedicated On-Board Computer in Active Debris Removal Mission" by Michaël Juillard & Jean-Paul Kneib (35th Proceedings of the AIAA/USU Conference on Small Satellites) [PDF}(https://www.researchgate.net/publication/354784632_Optimal_Control_Approach_for_Dedicated_On-Board_Computer_in_Active_Debris_Removal_Mission)


## 4 - Credits
Michael Juillard
