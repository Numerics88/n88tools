from __future__ import division
import os
import unittest
import numpy as np
from .config_regression import cfg
import shutil, tempfile
import subprocess


class TestExtractFieldsRegression(unittest.TestCase):
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

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_extractfields(self):
        '''Can run `n88extractfields` on a file'''
        command = ['n88extractfields', 'Displacement', self.model_filename]
        output = subprocess.check_output(command).decode("utf-8")
        self.assertNotEqual(output, '')

    def test_extractfields_displacement(self):
        '''Can run `n88extractfields Displacement` on a file'''
        command = ['n88extractfields', 'Displacement', self.model_filename]
        output = subprocess.check_output(command).decode("utf-8").split('\n')

        first_line = np.array([float(x) for x in output[0].split('\t')])
        M = np.array([-0.0012419662390147116, -0.000836902353299627, 0.0])

        self.assertTrue(np.allclose(first_line, M), '{}{}{}'.format(first_line, os.linesep, M))

    def test_extractfields_displacement_nodenumber(self):
        '''Can run `n88extractfields NodeNumber,Displacement` on a file'''
        command = ['n88extractfields', 'NodeNumber,Displacement', self.model_filename]
        output = subprocess.check_output(command).decode("utf-8").split('\n')

        first_line = np.array([float(x) for x in output[7162].split('\t')])
        M = np.array([7163, -8.00662123123e-05, -0.000160804893964, -0.0046864000638 ])

        self.assertTrue(np.allclose(first_line, M), '{}{}{}'.format(first_line, os.linesep, M))


if __name__ == '__main__':
    unittest.main()
