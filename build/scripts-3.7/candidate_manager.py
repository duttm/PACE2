from radical.entk import Pipeline, Stage, Task, AppManager
from typing import Dict, List

class CandidateManager:
    """
    Defines class to manage and run candidates in specified pipeline using radical.entk.
    """
    def __init__(self, hostname: str, port: int, username: str, password: str, resource_dict: Dict, pipelines: List[Pipeline]):
        """
        Initializes hostname, port, resource dictionary, and pipelines to run.
        Args:
            hostname: hostname (e.g. localhost)
            port: port number
            resource_dict: dictionary specifying resource configurations
            pipelines: list of candidates to run
        """
        self.hostname = hostname
        self.port = port
        self.username = username
        self.password = password
        self.app_manager = AppManager(hostname=hostname, port=port, username=username, password=password)
        self.app_manager.resource_desc = resource_dict
        self.app_manager.workflow = set(pipelines)

    def run(self) -> None:
        """
        Runs pipelines in workflow.
        """
        self.app_manager.run()
