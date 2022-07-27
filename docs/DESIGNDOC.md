## Dev Notes

### The following section is meant for developers of PACE. 


![Pace Schematic](/Images/PACE.png) 

### Objectives, broad design

- Modular, OOP based PACE code
  - Candidate class: 
    - Overview: Candidate class defines the candidate object, which is a unique set of phases (pre-md, md, analysis), each of which is composed of stages (currently just 1 stage for each phase). This unique ordered list of stages can also be abstracted as a "cycle." Each Candidate object is also exactly one radical.entk pipeline (note: if there are multiple cycles for this candidate, the pipeline contains all cycles).
    - Input: Dictionary of candidate phase specifications and candidate ID. 
    - Output: EnTK pipeline, containing the phases (which contains stages) which defines the candidate.
    - Uses the information from `simconfig.json` and `resconfig.json` to generate the required stages as an ordered list of EnTK stage objects, and pushes them to the pipeline.
    - Only `create_candidate_pipeline` is public for the purposes of proper object abstraction. The remaining functions to create md stages, analysis stages, and to extend pipeline are privately invoked by `create_candidate_pipeline`.
    - The "Phase" is an abstraction that may consist of one or more EnTK stages. This is to accommodate pre-exec needs of specific MD kernels.
    - Checks whether each pipeline needs to terminate / extend based on selection criteria (i.e. checks if cycles are needed).
        - If analysis condition fails: extends pipeline with new cycle.
        - If analysis condition succeeds: Nothing happens.
    - Repeats the above iteratively until either 1. Max number of iterations has been reached or 2. All candidates have been exhausted
  - Candidate manager class:
    - Overview: Manages execution of candidate objects and handles resource allocation to candidates.
    - Input: Candidates (where each candidate is a unique pipeline), elimination criteria (which defines when the candidates terminate vs. extend based on the analysis result) 
    - Inherits from EnTK AppManager
  - Top level driver script
    - Reads `simconfig.json` and `resconfig.json`.
    - Creates Candidate objects for each of the candidates specified, and adds them to a pipelines list.
    - Creates an instance of the Candidate Manager and runs all of the candidates.

## Dev Notes

### The following section is meant for developers of PACE. 


![Pace Schematic](/Images/PACE.png) 

### Objectives, broad design

- Modular, OOP based PACE code
  - Candidate class: 
    - Overview: Candidate class defines the candidate object, which is a unique set of phases (pre-md, md, analysis), each of which is composed of stages (currently just 1 stage for each phase). This unique ordered list of stages can also be abstracted as a "cycle." Each Candidate object is also exactly one radical.entk pipeline (note: if there are multiple cycles for this candidate, the pipeline contains all cycles).
    - Input: Dictionary of candidate phase specifications and candidate ID. 
    - Output: EnTK pipeline, containing the phases (which contains stages) which defines the candidate.
    - Uses the information from `simconfig.json` and `resconfig.json` to generate the required stages as an ordered list of EnTK stage objects, and pushes them to the pipeline.
    - Only `create_candidate_pipeline` is public for the purposes of proper object abstraction. The remaining functions to create md stages, analysis stages, and to extend pipeline are privately invoked by `create_candidate_pipeline`.
    - The "Phase" is an abstraction that may consist of one or more EnTK stages. This is to accommodate pre-exec needs of specific MD kernels.
    - Checks whether each pipeline needs to terminate / extend based on selection criteria (i.e. checks if cycles are needed).
        - If analysis condition fails: extends pipeline with new cycle.
        - If analysis condition succeeds: Nothing happens.
    - Repeats the above iteratively until either 1. Max number of iterations has been reached or 2. All candidates have been exhausted
  - Candidate manager class:
    - Overview: Manages execution of candidate objects and handles resource allocation to candidates.
    - Input: Candidates (where each candidate is a unique pipeline), elimination criteria (which defines when the candidates terminate vs. extend based on the analysis result) 
    - Inherits from EnTK AppManager
  - Top level driver script
    - Reads `simconfig.json` and `resconfig.json`.
    - Creates Candidate objects for each of the candidates specified, and adds them to a pipelines list.
    - Creates an instance of the Candidate Manager and runs all of the candidates.

