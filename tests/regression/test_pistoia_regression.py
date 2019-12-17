from __future__ import division
import os
import unittest
from n88tools import transformations
import numpy as np
import math
from config_regression import cfg
import shutil, tempfile
import subprocess
import vtkbone
import re


class TestPistoiaRegression(unittest.TestCase):
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
        command = ['n88pistoia', os.path.join(self.test_dir, self.filename)]
        self.output = subprocess.check_output(command)

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_pistoia(self):
        '''Can run `n88pistoia` on a file'''
        self.assertNotEqual(self.output, '')

    def test_pistoia_failure_load(self):
        '''`n88pistoia` returns correct failure load'''
        expected = np.array([-1.3131E-04, 6.8915E-05, -4.9765E+00])

        expression=r'\s+Failure load \(RF \* factor\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = np.array([float(result[0][i]) for i in [0,4,8]])
        self.assertTrue(np.allclose(result, expected), '{}{}{}'.format(result, os.linesep, expected))

    def test_pistoia_critical_volume(self):
        '''`n88pistoia` returns correct critical volume'''
        expected = 2.0

        expression=r'\s+Critical volume \(\%\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_critical_ees(self):
        '''`n88pistoia` returns correct critical ees'''
        expected = 7.0000E-03

        expression=r'\s+Critical EES:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_at_vol_crit(self):
        '''`n88pistoia` returns correct ess at vol_crit'''
        expected = 1.4333E-02

        expression=r'\s+EES at vol_crit:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_factor(self):
        '''`n88pistoia` returns correct factor'''
        expected = 4.8838E-01

        expression=r'\s+Factor \(from table\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_RF(self):
        '''`n88pistoia` returns correct RF'''
        expected = np.array([-2.6886E-04, 1.4111E-04, -1.0190E+01])

        expression=r'\s+RF \(node set 1\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = np.array([float(result[0][i]) for i in [0,4,8]])
        self.assertTrue(np.allclose(result, expected), '{}{}{}'.format(result, os.linesep, expected))

    def test_pistoia_U(self):
        '''`n88pistoia` returns correct U'''
        expected = np.array([-6.5332E-04, 2.0151E-04, -8.5000E-03])

        expression=r'\s+U \(node set 1\):\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = np.array([float(result[0][i]) for i in [0,4,8]])
        self.assertTrue(np.allclose(result, expected), '{}{}{}'.format(result, os.linesep, expected))

    def test_pistoia_axial_stiffness(self):
        '''`n88pistoia` returns correct axial stiffness'''
        expected = np.array([4.1154E-01, 7.0026E-01, 1.1988E+03])

        expression=r'\s+Axial stiffness:\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) \s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = np.array([float(result[0][i]) for i in [0,4,8]])
        self.assertTrue(np.allclose(result, expected), '{}{}{}'.format(result, os.linesep, expected))

    def test_pistoia_ess_distribution_average(self):
        '''`n88pistoia` returns correct ESS distribution average'''
        expected = 5.446E-03

        expression=r'\s+average\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_distribution_std_dev(self):
        '''`n88pistoia` returns correct ESS distribution std dev'''
        expected = 3.755E-03

        expression=r'\s+std_dev\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_distribution_minimum(self):
        '''`n88pistoia` returns correct ESS distribution minimum'''
        expected = 7.248E-06

        expression=r'\s+minimum\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_distribution_maximum(self):
        '''`n88pistoia` returns correct ESS distribution maximum'''
        expected = 2.448E-02

        expression=r'\s+maximum\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_distribution_skewness(self):
        '''`n88pistoia` returns correct ESS distribution skewness'''
        expected = 8.473E-01

        expression=r'\s+skewness\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_distribution_kurtosis(self):
        '''`n88pistoia` returns correct ESS distribution kurtosis'''
        expected = 4.639E-01

        expression=r'\s+kurtosis\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)

    def test_pistoia_ess_distribution_median(self):
        '''`n88pistoia` returns correct ESS distribution median'''
        expected = 4.760E-03

        expression=r'\s+median\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        result = float(result[0][0])
        self.assertAlmostEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
