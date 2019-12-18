from __future__ import division
import os
import unittest
import numpy as np
from .config_regression import cfg
import shutil, tempfile
import subprocess
import re


class TestPistoiaRegressionCritical(unittest.TestCase):
    filename = 'test25a_uniaxial_solved.n88model'

    def setUp(self):
        # Create temporary directory to work in
        self.test_dir = tempfile.mkdtemp()

        # Fetch the data
        download_location = cfg['DOWNLOAD_TESTING_DATA'](self.filename)
        self.assertNotEqual(download_location, '', 'Unable to download file ' + self.filename)

        # Copy to temporary directory
        shutil.copy(download_location, self.test_dir)
        self.assertTrue(os.path.isfile(os.path.join(self.test_dir, self.filename)))

        # Run pistoia
        command = ['n88pistoia', '-v', '3.0', '-s', '0.01', os.path.join(self.test_dir, self.filename)]
        self.output = subprocess.check_output(command)

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_pistoia(self):
        '''Can run `n88pistoia` on a file'''
        self.assertNotEqual(self.output, '')

    def test_pistoia_failure_load(self):
        '''`n88pistoia` returns correct failure load'''
        expected = np.array([-2.0027E-04, 1.0511E-04, -7.5904E+00])

        expression=r'\s+Failure load \(RF \* factor\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = np.array([float(result[0][i]) for i in [0,4,8]])
        self.assertTrue(np.allclose(result, expected), '{}{}{}'.format(result, os.linesep, expected))

    def test_pistoia_critical_volume(self):
        '''`n88pistoia` returns correct critical volume'''
        expected = 3.0

        expression=r'\s+Critical volume \(\%\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_critical_ees(self):
        '''`n88pistoia` returns correct critical ees'''
        expected = 1.0000E-02

        expression=r'\s+Critical EES:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_at_vol_crit(self):
        '''`n88pistoia` returns correct ess at vol_crit'''
        expected = 1.3425E-02

        expression=r'\s+EES at vol_crit:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_factor(self):
        '''`n88pistoia` returns correct factor'''
        expected = 7.4489E-01

        expression=r'\s+Factor \(from table\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
