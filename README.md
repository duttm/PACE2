# PACE<sup>2</sup>: Pipeline for Automating Compliance-based Elimination and Extension

PACE<sup>2</sup> (sometimes also written PACE2) is a framework written in the RADICAL-EnTK API (EnTK) (https://github.com/radical-cybertools/radical.entk). PACE<sup>2</sup> has been designed to create and run various flavours of simulation-analysis loops, convergence checks and output-targeted simulations. 

Candidate systems are fed to PACE<sup>2</sup> as directories containing all inputs needed for an MD run. PACE<sup>2</sup> itself invokes the MD engine to run the simulations. After the simulations complete, PACE<sup>2</sup> calls a user-defined script to check compliance of each MD system. If a system meets the user-specified compliance requirements, the system is extended within the pipeline where it can continue processing. If the compliance requirements are not met, the system is eliminated and its resources are returned to the pool so that they can be reassigned to a different system. In this way PACE<sup>2</sup> manages resources to adapt to the needs of the user without needing a user in the loop.

See the Examples directory for examples running GROMACS and using custom analysis scripts.

## Installation

We recommend executing in a virtual environment using e.g. venv. Create and activate a Python 3.9+ environment using venv as follows. (This is just an example. Check your environment manager's docs.)

```
python3.9 -m virtualenv my_pace2_env
source my_pace2_env/bin/activate 
```

Clone and install the repository using:

```
git clone https://github.com/duttm/PACE2.git
cd PACE2/
pip install -e .
```

### Other runtime considerations

EnTK is installed automatically during the PACE<sup>2</sup> installation process, but EnTK still has some dependencies which must be configured. We recommend reviewing the EnTK Github documentation.

PACE<sup>2</sup> is tested targeting execution on remote hpc clusters using SSH authentication, using the RCT Radial-Pilot (RP) system (https://github.com/radical-cybertools/radical.pilot). RP is installed automatically during the PACE<sup>2</sup> installation process.

## Usage

To use PACE<sup>2</sup>:
1. Prepare the Candidate Pool (directories)
2. Prepare simulation configuration and resource configuration (JSON files)
3. Execute via the Linux command line

### Step 1
The Candidate is the fundamental unit of the PACE<sup>2</sup> adaptivity model. The full Candidate Pool is organized as a directory structure, with the main directory representing the Pool, and subdirectories representing individual Candidates. Each Candidate directory includes all of the files necessary to execute an MD simulation. 
    The Candidate Pool directory is given a user-defined "basename", and the subdirectories are named "basename.0", "basename.1", etc. The basename string must obey the character requirements of several upstream pieces of software, so we recommend using a simple name. See the /examples directory for a better understanding of how this should look.

### Step 2
Simulation Configuration and Resource Configuration parameters are organized in a pair of JSON files conventionally called simconfig.json and resconfig.json. 

See associated files in the /examples directory. Here is a brief description of keys in each of the config files.

**simconfig.json keys and descriptions** 
* "basename" : A name for the Candidate Pool. See Step 1 instructions.
* "candidates" : Number of Candidates in the Candidate Pool.
* "cycle_max" : Maximum number of pre-MD/MD/analysis cycles.
* "pipeline_cores" : Number of cores on the compute resource to dedicate to the run. 
* "pre_md_pre_exec" : A command to be executed immediately before pre-MD, e.g. "module load gromacs/2018".
* "pre_md_executable" : The command line program to be used for pre-MD, e.g. "/usr/bin/gmx". Use fully qualified paths. N.B. that these paths are to be parsed on the compute resource.
* "pre_md_args" : Arguments for the pre-MD program, e.g. "grompp -f FFFNF.mdp -c FFFNF.gro -o sys.tpr -p FFFNF.top"
* "md_pre_exec" : A command to be executed immediately before MD, e.g. "module load openmpi5.5".
* "md_executable" : The command line program to be used for pre-MD, e.g. "/usr/bin/gmx". Use fully qualified paths. N.B. that these paths are to be parsed on the compute resource.
* "md_args" : Arguments for the MD program, e.g. "mdrun -s sys.tpr -deffnm sys -c outcrd.gro"
* "an_pre_exec" : A command to be executed immediately before analysis, e.g. "module load python3.9.x"
* "an_executable" : The command line program to be used for analysis, e.g. "/usr/local/bin/python3.9". Use fully qualified paths. N.B. that these paths are to be parsed on the compute resource. 
* "an_args" : Arguments for the analysis program, e.g. "/$HOME/PACE2/examples/dummy_example/force_outcome.py 1",
* "chead_files" : A list of all user-furnished files in each Candidate, e.g. "FFFNF.mdp FFFNF.top FF.itp FNF.itp WF.itp FFFNF.gro martini_v2.2.itp",
* "md_binary" : An intermediate file which constitutes the standalone MD binary (this is a GROMACS thing), e.g. "sys.tpr",
* "structure_in" :  The filename of the initial structure in each Candidate. This file gets special handling. e.g. "FFFNF.gro",
* "structure_out" : The filename of the final structure output during each MD phase. This file gets special handling. e.g. "outcrd.gro"

**resconfig.json keys and descriptions** 
* "resource" : A Radical-Pilot resource key.
* "walltime" : Maximum number of minutes to run.
* "cpus" : Number of cpus to be used on the compute resource.

### Step 3
With candidate directories set up, run PACE<sup>2</sup> (**remember to activate your virtual environment!**) using:

```
pace2.py simconfig.json resconfig.json
```

`simconfig.json` is the *simulation configuration* file. This is where simulation parameters are specified.
`resconfig.json` is the *resource configuration* file. This is where the resource, i.e. localhost/Bridges2/Expanse etc. are specified. 

## Issues and Support

If you believe you have found a bug, please open a ticket at https://github.com/duttm/PACE2.

For other feedback, please send an email to the authors (copy all of them for better response):

Srinivas Mushnoori srinivas.mushnoori@gmail.com

Akash Banerjee akashbanerji@gmail.com

Mason Hooten mason.simulation@gmail.com

## Acknowledgements

The developers gratefully acknowledge financial support from the National Science Foundation CAREER award DMR-1654325, awards OAC-1547580, OAC-1835449 and DMR-2118860.

Thanks for using PACE<sup>2</sup>!
