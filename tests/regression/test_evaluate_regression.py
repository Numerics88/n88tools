from __future__ import division
import os
import unittest
from config_regression import cfg
import shutil, tempfile
import subprocess
import re


class TestEvaluateRegression(unittest.TestCase):
    filename = 'test25a_uniaxial_solved.n88model'

    def setUp(self):
        # Create temporary directory to work in
        self.test_dir = tempfile.mkdtemp()

        # Fetch the data
        input_uncompress_model = cfg['DOWNLOAD_TESTING_DATA'](self.filename)
        self.assertNotEqual(input_uncompress_model, '', 'Unable to download file ' + self.filename)

        # Copy to temporary directory
        shutil.copy(input_uncompress_model, self.test_dir)
        self.model_filename = os.path.join(self.test_dir, self.filename)
        self.assertTrue(os.path.isfile(self.model_filename))

        # Run the command
        command = ['n88evaluate', self.model_filename]
        self.output = subprocess.check_output(command)

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_evaluate(self):
        '''Can run `n88evaluate` on a file'''
        self.assertNotEqual(self.output, '')

    def test_solution_displacement_max_error(self):
        '''`n88evaluate` gives correct solution displacement max error'''
        expression=r'\s+max err\s*:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.00E+00)

    def test_solution_displacement_rms_error(self):
        '''`n88evaluate` gives correct solution displacement rms error'''
        expression=r'\s+rms err\s*:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.00E+00)

    def test_force_max_error(self):
        '''`n88evaluate` gives correct force max error'''
        expression=r'\s+max err\s*:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[1][0]), 5.02E-05)

    def test_force_rms_error(self):
        '''`n88evaluate` gives correct force rms error'''
        expression=r'\s+rms err\s*:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[1][0]), 9.92E-06)

    def test_force_max_relative_error(self):
        '''`n88evaluate` gives correct relative force max error'''
        expression=r'\s+max err/max force\s*:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 5.86E-04)

    def test_force_rms_relative_error(self):
        '''`n88evaluate` gives correct relative force rms error'''
        expression=r'\s+rms err/max force\s*:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 1.16E-04)


if __name__ == '__main__':
    unittest.main()
