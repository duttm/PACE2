# Lipid Example

This example uses PACE<sup>2</sup> to generate lipid-based vesicles. 

## Workflow 

MD simulation -> Analysis -> Extended MD Simulation 

## Analysis Kernel

The analysis kernel for this example hosts a k-means clustering code. The code analyzes the output of the MD simulation and provides a value of either 1 or 0 in a text file. If the value is 1, the simulation is extended. If the value is 0, the simulation is terminated. 

## Candidates

This example has 2 candidates in the `DV` folder. 
The first candidate (`DV.0`)  will result in simulation termination, and the second candidate (`DV.1`) will result in simulation extension. 
