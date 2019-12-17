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


class TestExtractSets(unittest.TestCase):
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

        # Run
        command = ['n88extractsets', self.model_filename]
        self.output = subprocess.check_output(command)

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_extractsets(self):
        '''Can run `n88extractsets` on a file'''
        self.assertNotEqual(self.output, '')

    def test_extractsets_outputfiles(self):
        '''`n88extractsets` produces files'''
        output_extensions = [
             '_constraint_bottom_fixed.vtp'
            ,'_constraint_top_displacement.vtp'
            ,'_node_set_face_z0.vtp'
            ,'_node_set_face_z1.vtp'
            ,'_node_set_face_x0.vtp'
            ,'_node_set_face_x1.vtp'
            ,'_node_set_face_y0.vtp'
            ,'_node_set_face_y1.vtp'
            ,'_element_set_face_z0.vtp'
            ,'_element_set_face_z1.vtp'
            ,'_element_set_face_x0.vtp'
            ,'_element_set_face_x1.vtp'
            ,'_element_set_face_y0.vtp'
            ,'_element_set_face_y1.vtp'
        ]

        for extension in output_extensions:
            filename = self.model_filename.replace('.n88model', extension)
            self.assertTrue(os.path.isfile(filename), 'Did not produce file ' + filename)


if __name__ == '__main__':
    unittest.main()
