## Amber Example

Here, we test a short energy minimization simulation of a protein. We use the amber MD engine for this task.

We have used 3 input files for this workflow: min.in (MD input parameters) , prmtop and prmcrd (coordinate and topology files). These files are stores in `amber_example/DV/DV.0/`

The PACE workflow controls the simulation with the help of `simconfig_local.json` (input files, executables, etc) and `resconfig_local.json` (walltime, target resource for running the simulation)

Here are the relevant keys (from `simconfig_local.json`) for this example:


```
    "pre_md_executable" : "sander", ## Amber execuatable
    "pre_md_args"       : "-i min.in -o 6pti.min1.out -c prmcrd -r 6pti.min1.xyz", ## Amber/sander commands to run an amber simulation
    "chead_files"       : "min.in prmcrd prmtop", ## list of simulation input files
    "md_binary"         : "6pti.min1.xyz", ## simulation output file (amber restart file)
    "structure_in"     :  "prmcrd", ## simulation input file 
    "structure_out"     : "6pti.min1.xyz" ## simulation output file (amber restart file)
```

We note that PACE has been designed for separate pre-md and md stage. This is applicable for MD engines such as GROMACS. However, the amber MD engine does
not require 2 separate stages to run an MD simulation. Hence, we have run the amber simulation in the pre-md stage, and the md-stage is left empty. In future versions 
of PACE, we will fix this design issue. 
