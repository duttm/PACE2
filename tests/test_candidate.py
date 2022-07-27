import os
import sys
import unittest

from unittest.mock import Mock, patch

sys.modules['radical.entk'] = Mock()
from src import candidate


class TestCandidate(unittest.TestCase):
    """
    Unit test cases for candidate.py constructor.
    """

    def test_init(self):
        # Define testing parameters
        test_sim_config = {'basename': 'basename',
                           'candidates': 1,
                           'cycle_max': 1,
                           'pipeline_cores': 1,
                           'md_executable': 'md_executable',
                           'md_pre_exec': 'md_pre_exec',
                           'md_args': 'md_args',
                           'an_executable': 'an_executable',
                           'an_args': 'an_args',
                           'pilot_cores': 'pilot_cores'
                           }
        test_cid = Mock()
        test_candidate = candidate.Candidate(test_sim_config, test_cid)

        # Test attributes of candidate object
        self.assertDictEqual(test_candidate.candidate_specifications, test_sim_config)
        self.assertEqual(test_cid, test_candidate.cid)
        self.assertEqual(0, test_candidate.cycle_count)

    @patch('src.candidate.Candidate._create_analysis_stage')
    @patch('src.candidate.Candidate._create_md_stage')
    def test_create_candidate_pipeline(self, mock_create_md_stage, mock_create_analysis_stage):
        # Define testing parameters
        test_sim_config = {'basename': 'basename',
                           'candidates': 1,
                           'cycle_max': 1,
                           'pipeline_cores': 1,
                           'md_executable': 'md_executable',
                           'md_pre_exec': 'md_pre_exec',
                           'md_args': 'md_args',
                           'an_executable': 'an_executable',
                           'an_args': 'an_args',
                           'pilot_cores': 'pilot_cores'
                           }
        test_cid = Mock()
        test_candidate = candidate.Candidate(test_sim_config, test_cid)
        test_create_candidate_pipeline_return = test_candidate.create_candidate_pipeline()

        # Check that _create_md_stage, and _create_analysis_stage were each called once
        mock_create_md_stage.assert_called_once()
        mock_create_analysis_stage.assert_called_once()

        # Check that the return type is of type mock, which was used to patch radical.entk.Pipeline()
        self.assertIsInstance(test_create_candidate_pipeline_return, unittest.mock.Mock)


if __name__ == '__main__':
    unittest.main()
