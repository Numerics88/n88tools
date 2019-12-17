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
from netCDF4 import Dataset


class TestInterpolateSolution(unittest.TestCase):
    filenames = [
         'test25a_uniaxial_solved.n88model'
        ,'test25a_uniaxial.n88model'
        ,'test25a_uniaxial_coarse.n88model'
    ]

    def setUp(self):
        # Create temporary directory to work in
        self.test_dir = tempfile.mkdtemp()

        for filename in self.filenames:
            # Fetch the data
            input_uncompress_model = cfg['DOWNLOAD_TESTING_DATA'](filename)
            self.assertNotEqual(input_uncompress_model, '', 'Unable to download file ' + filename)

            # Copy to temporary directory
            shutil.copy(input_uncompress_model, self.test_dir)
            model_filename = os.path.join(self.test_dir, filename)
            self.assertTrue(os.path.isfile(model_filename))

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_interpolatesolution(self):
        '''Can run `n88interpolatesolution` on a file'''
        coarse_model = os.path.join(self.test_dir, 'test25a_uniaxial_coarse.n88model')
        unsolved_model = os.path.join(self.test_dir, 'test25a_uniaxial.n88model')

        command = ['n88interpolatesolution', unsolved_model, coarse_model]
        output = subprocess.check_output(command)
        self.assertNotEqual(output, '')

        # Check if readable
        errorObserver = cfg['ERROR_OBSERVER']()
        reader = vtkbone.vtkboneN88ModelReader()
        reader.AddObserver ("ErrorEvent", errorObserver)
        reader.SetFileName(unsolved_model)
        reader.Update()
        self.assertFalse(errorObserver.ErrorOccurred(), errorObserver.ErrorMessage())


if __name__ == '__main__':
    unittest.main()
