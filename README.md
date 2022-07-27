# PACE<sup>2</sup>: Pipeline for Automating Compliance-based Elimination and Extension

PACE<sup>2</sup> (sometimes also written PACE2) is a framework written in the RADICAL-EnTK API (https://github.com/radical-cybertools/radical.entk). PACE<sup>2</sup> has been designed to create and run various flavours of simulation-analysis loops, convergence checks and output-targeted simulations. 

Candidate systems are fed to PACE<sup>2</sup> as directories containing all inputs needed for an MD run. PACE<sup>2</sup> itself invokes the MD engine to run the simulations. After the simulations complete, PACE<sup>2</sup> calls applies a user-defined script to check compliance of each MD system. If a system meets the user-specified compliance requirements, the system is extended within the pipeline where it can continue processing. If the compliance requirements are not met, the system is eliminated and its resources are returned to the pool so that they can be reassigned to a different system. In this way PACE<sup>2</sup> manages resources to adapt to the needs of the user without needing a user in the loop.

See the Examples directory for examples running GROMACS and using custom analysis scripts.

## Installation

We recommend using a virtual environment using e.g. venv. Create and activate a Python 3.9+ environment using venv as follows. (This is just an example. Check your environment manager's docs.)

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

### RabbitMQ

PACE<sup>2</sup> uses the RADICAL Cybertools (RCT) suite of tools. RCT depends on a program called RabbitMQ, a message broker service. See the Support section below if you have trouble getting RMQ set up.  

## Usage

Run PACE<sup>2</sup> (inside the environment) using:

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

Mason Hooten mason.hooten@gmail.com

## Acknowledgements

The developers gratefully acknowledge financial support from the National Science Foundation CAREER award DMR-1654325, awards OAC-1547580 and DMR-2118860.

Thanks for using PACE<sup>2</sup>!
