#!/usr/bin/env python


import os, sys, json

from candidate import Candidate
from candidate_manager import CandidateManager
from typing import Dict, Tuple


def read_jsons() -> Tuple[Dict, Dict]:
    """
    Read from specified resource and simulation json configuration files.
    Returns:
        Tuple[Dict, Dict] specifying the simulation and resource configurations, respectively
    """
    simconfig = sys.argv[1]  # simulation configuration
    resconfig = sys.argv[2]  # resource configuration
    sim_dict = {}

    # Load simulation and resource json files
    with open(simconfig) as simconf:
        simdata = json.load(simconf)

    with open(resconfig) as resconf:
        resdata = json.load(resconf)

    # Assign simulation configs
    sim_dict['basename'] = simdata["basename"]
    sim_dict['candidates'] = simdata["candidates"]
    sim_dict['cycle_max'] = simdata["cycle_max"]
    sim_dict['pipeline_cores'] = simdata["pipeline_cores"]

    # Executables, args and pre exec for all the tasks
    
    sim_dict['md_executable']     = simdata["md_executable"]
    sim_dict['md_args']           = simdata["md_args"].split()
    sim_dict['md_pre_exec']       = simdata["md_pre_exec"]

    sim_dict['pre_md_executable'] = simdata["pre_md_executable"]
    sim_dict['pre_md_args']       = simdata["pre_md_args"].split()
    sim_dict['pre_md_pre_exec']   = simdata["pre_md_pre_exec"]
        
    sim_dict['an_pre_exec']       = simdata["an_pre_exec"]
    sim_dict['an_executable']     = simdata["an_executable"]
    sim_dict['an_args']           = simdata["an_args"].split()




    sim_dict['pilot_cores'] = resdata["cpus"]

    # Assign resource configs when all fields are specified
    try:
        res_dict = {
            "resource": str(resdata["resource"]),
            "walltime": int(resdata["walltime"]),
            "cpus": sim_dict['pilot_cores'],
            "gpus_per_node": int(resdata["gpus_per_node"]),
            "access_schema": str(resdata["access_schema"]),
            "queue": str(resdata["queue"]),
            "project": str(resdata["project"]),
        }

    # Assign resource configs when minimum fields specified
    except:
        res_dict = {
            "resource": str(resdata["resource"]),
            "walltime": int(resdata["walltime"]),
            "cpus": sim_dict['pilot_cores'],

        }
    return sim_dict, res_dict


def main() -> None:
    """
    Main pace driver.
    """
    # Read from simconfig.json
    candidate_specifications_dict, resource_dict = read_jsons()

    # Obtain pre-run environment variables
    if 'RADICAL_ENTK_VERBOSE' in os.environ:
        os.environ['RADICAL_ENTK_REPORT'] = 'True'
    hostname = os.environ.get('RMQ_HOSTNAME', 'localhost')
    port = os.environ.get('RMQ_PORT', 32769)
    username = os.environ.get('RMQ_USERNAME')
    password = os.environ.get('RMQ_PASSWORD')

    # Generate all pipelines
    pipelines = []
    for cid in range(candidate_specifications_dict['candidates']):
        # Create candidate pipeline
        candidate = Candidate(candidate_specifications_dict, cid)
        candidate.create_candidate_pipeline()
        pipelines.append(candidate.pipeline)

    # Create candidate manager and run
    candidate_manager = CandidateManager(hostname, port, username, password, resource_dict, pipelines)
    candidate_manager.run()


if __name__ == '__main__':
    main()
