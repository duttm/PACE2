import os

from radical.entk import Pipeline, Stage, Task, AppManager
from typing import Dict

class Candidate:
    """
    Defines class to manage candidates.
    """
    def __init__(self, candidate_specifications: Dict, cid: int):
        """
        Initializes pipeline, candidate specifications dictionary, and candidate id
        Args:
            candidate_specifications: dictionary containing simulation configurations
            cid: unique candidate id
        """
        self.pipeline = Pipeline()
        self.candidate_specifications = candidate_specifications
        self.cid = cid
        self.cycle_count = 0

    def create_candidate_pipeline(self) -> Pipeline:
        """
        Creates candidate pipeline with appropriate stages
        Returns:
            Pipeline: pipeline containing appropriate stages (pre md, md, analysis, etc.) for the candidate
        """
        basename = self.candidate_specifications['basename']
        sysname = basename + "." + str(self.cid)
        self.pipeline.name = sysname

        # Create md stage
        self._create_md_stage()
        # Create analysis stage
        analysis_stage = self._create_analysis_stage()
        # Extend pipeline in post exec as necessary
        analysis_stage.post_exec = self._extend_pipeline

        return self.pipeline

    def _extend_pipeline(self) -> None:
        """
        Extends or terminates pipeline according to analysis output and max cycles threshold.
        """
        basename = self.candidate_specifications['basename']
        cycle_max = int(self.candidate_specifications['cycle_max'])
        sysname = basename + "." + str(self.cid)

        file_number = 0
        temp_filename = 'out{}.{}.txt'.format(sysname,file_number)

        # Find the most recently generated output analysis file
        while os.path.isfile(temp_filename):
            analysis_output_filename = temp_filename
            file_number += 1
            temp_filename = 'out{}.{}.txt'.format(sysname,file_number)

        # Open most recently generated output analysis file
        with open(analysis_output_filename, "r") as f:
            analysis_output = f.read()
            print('Output has been read as', analysis_output,'from file',f.name)

            # Extend if analysis indicates to do so and max cycles has not been reached
            if analysis_output == "1" and self.cycle_count < cycle_max:
                self.cycle_count += 1

                print('Extending pipeline.')
                self._create_md_stage()
                analysis_stage = self._create_analysis_stage()
                analysis_stage.post_exec = self._extend_pipeline

    def _create_pre_md_stage(self) -> None:
        """
        Creates pre md stage
        """
        basename = self.candidate_specifications['basename']
        sysname = basename + "." + str(self.cid)
        
        c_head_extension = self.candidate_specifications['chead_files']    
       
        c_head_initial = basename + "/" + sysname + "/"
        
        c_head = []
        
        for i in range(len(c_head_extension)):
            c_head.append(c_head_initial + c_head_extension[i])
            ##print(c_head[i])
        
#        c_head = [basename + "/" + sysname + "/" + basename + ".gro",
#                  basename + "/" + sysname + "/" + basename + ".itp",
#                  basename + "/" + sysname + "/" + basename + ".mdp",
#                  basename + "/" + sysname + "/" + basename + ".top",
#                  basename + "/" + sysname + "/" + "martini_v2.2.itp"]
                        

        # Create pre md stage
        pre_md_stage = Stage()
        pre_md_stage.name = 'premdstage' + str(self.cycle_count)
        # Create pre md task
        pre_md_task = Task()
        pre_md_task.name = 'premdtask'  ####
        pre_md_task.executable = self.candidate_specifications['pre_md_executable']
        pre_md_task.arguments = self.candidate_specifications['pre_md_args']
        pre_md_task.pre_exec = self.candidate_specifications['pre_md_pre_exec']
        if self.cycle_count == 0:
            pre_md_task.upload_input_data = c_head
        else:
            # Take the same file names from c_head, no sysname/basename
            
            
            c_head_2 = []
            for i in range(len(c_head_extension)):
                ##print(c_head_extension[i])
                if c_head_extension[i] != self.candidate_specifications['structure_in']:
                    c_head_2.append('$Pipline_%s_Stage_%s_Task_%s/%s'% (sysname, 'premdstage0', 'premdtask', c_head_extension[i]))    
                
               
                else:     
                    c_head_2.append('$Pipline_%s_Stage_%s_Task_%s/%s > %s'% (sysname, 'mdstage'+str(self.cycle_count-1), 'mdtask', self.candidate_specifications['structure_out'], self.candidate_specifications['structure_in']))
                    
            pre_md_task.link_input_data = c_head_2        
            
#            pre_md_task.link_input_data = ['$Pipline_%s_Stage_%s_Task_%s/outcrd.gro > %s.gro' % (sysname, 'mdstage'+str(self.cycle_count-1), 'mdtask', basename),
#                                           '$Pipline_%s_Stage_%s_Task_%s/%s.itp'              % (sysname, 'premdstage0', 'premdtask', basename),
#                                           '$Pipline_%s_Stage_%s_Task_%s/%s.mdp'              % (sysname, 'premdstage0', 'premdtask', basename),
#                                           '$Pipline_%s_Stage_%s_Task_%s/%s.top'              % (sysname, 'premdstage0', 'premdtask', basename),
#                                           '$Pipline_%s_Stage_%s_Task_%s/martini_v2.2.itp'    % (sysname, 'premdstage0', 'premdtask')]
        pre_md_stage.add_tasks(pre_md_task)
        # Add pre md stage to pipeline
        self.pipeline.add_stages(pre_md_stage)

    def _create_md_stage(self) -> None:
        """
        Creates md stage
        """

        # Create pre md stage
        self._create_pre_md_stage()
        # Create md stage
        md_stage = Stage()
        md_stage.name = 'mdstage' + str(self.cycle_count)
        # Create md task
        md_task = Task()
        md_task.name = 'mdtask'
        md_task.executable = self.candidate_specifications['md_executable']
        md_task.arguments = self.candidate_specifications['md_args']
        md_task.pre_exec = self.candidate_specifications['md_pre_exec']
        sysname = self.candidate_specifications['basename'] + "." + str(self.cid)
        
        md_task.link_input_data = ['$Pipline_%s_Stage_%s_Task_%s/%s' % (sysname, "premdstage" + str(self.cycle_count), "premdtask", self.candidate_specifications['md_binary'])]

        md_stage.add_tasks(md_task)
        # Add md stage to pipeline
        self.pipeline.add_stages(md_stage)

    def _create_analysis_stage(self) -> Stage:
        """
        Creates analysis stage
        """
        
        basename = self.candidate_specifications['basename']
        sysname = basename + "." + str(self.cid)
        
        # Create analysis stage
        an_stg = Stage()
        an_stg.name = 'analysisstage' + str(self.cycle_count)
        # Create the analysis task
        an = Task()
        an.name = 'analysistask'
        an.executable = self.candidate_specifications['an_executable']
        an.pre_exec = self.candidate_specifications['an_pre_exec']
        an.link_input_data = ['$Pipline_%s_Stage_%s_Task_%s/%s > %s' % (sysname, 'mdstage'+str(self.cycle_count), 'mdtask', self.candidate_specifications['structure_out'], self.candidate_specifications['structure_in'])]
        
        an_args = []
        # Add all analysis arguments in candidate specifications json
        for arg in self.candidate_specifications['an_args']:
            an_args.append(arg)
        an.arguments = an_args

        # Create file with success/failure output of analysis stage
        file_number = 0
        filename = 'out{}.{}.txt'.format(sysname,file_number)
        while os.path.isfile(filename):
            filename = 'out{}.{}.txt'.format(sysname,file_number)
            file_number += 1

        f = open(filename, 'x')
        f.close()
        
        # Copy output of analysis out.txt to unique file
        an.download_output_data = ['out.txt > {}'.format(filename)]
        # an.copy_input_data = ['$Pipline_%s_Stage_%s_Task_%s/file.tar.gz' % (p.name, md_stg.name, md.name),]
        an_stg.add_tasks(an)

        # Add analysis stage to pipeline
        self.pipeline.add_stages(an_stg)
        return an_stg
