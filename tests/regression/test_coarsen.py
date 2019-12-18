from __future__ import division
import os
import unittest
from config_regression import cfg
import shutil, tempfile
import vtkbone


class TestCoarsen(unittest.TestCase):
    filename = 'test25a_uniaxial.n88model'

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

    def test_coarsen(self):
        '''Can run `n88coarsen` on a file'''
        output_filename = os.path.join(self.test_dir, 'coarsen_model.n88model')

        # Run command
        command = ['n88coarsen', self.model_filename, output_filename]
        self.assertTrue(
            cfg['RUN_CALL'](command),
            'Cannot call \"{}\"'.format(' '.join(command))
        )

        # File should exist
        self.assertTrue(os.path.isfile(output_filename))

        # Check if readable
        errorObserver = cfg['ERROR_OBSERVER']()
        reader = vtkbone.vtkboneN88ModelReader()
        reader.AddObserver ("ErrorEvent", errorObserver)
        reader.SetFileName(output_filename)
        reader.Update()
        self.assertFalse(errorObserver.ErrorOccurred(), errorObserver.ErrorMessage())


if __name__ == '__main__':
    unittest.main()
