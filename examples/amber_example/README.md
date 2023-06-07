## Amber Example

Here, we test a short energy minimization simulation of a protein using the Amber MD engine. Please visit http://ambermd.org/ for details on obtaining and using the Amber software.

We have used 3 input files for this workflow: min.in (MD input parameters), prmtop, and prmcrd (coordinate and topology files). These files are stored in `amber_example/DV/DV.0/`

The PACE<sup>2</sup> workflow controls the simulation with the help of `simconfig_local.json` (input files, executables, etc) and `resconfig_local.json` (walltime and  target resource for running the simulation).

Here are the relevant keys (from `simconfig_local.json`) for this example:


```
    "pre_md_executable" : "sander", ## Amber execuatable
    "pre_md_args"       : "-i min.in -o 6pti.min1.out -c prmcrd -r 6pti.min1.xyz", ## Amber/sander commands to run an amber simulation
    "chead_files"       : "min.in prmcrd prmtop", ## list of simulation input files
    "md_binary"         : "6pti.min1.xyz", ## simulation output file (amber restart file)
    "structure_in"     :  "prmcrd", ## simulation input file 
    "structure_out"     : "6pti.min1.xyz" ## simulation output file (amber restart file)
```

We note that PACE<sup>2</sup> has been designed for separate pre-MD and MD stages. This is applicable for MD engines such as GROMACS. However, the Amber MD engine does
not require two separate stages to run an MD simulation. In this implementation, the Amber simulation runs in the pre-MD stage, and the MD stage is left empty. Future versions 
of PACE<sup>2</sup> will incorporate a single stage MD simulation feature. 
