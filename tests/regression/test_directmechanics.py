from __future__ import division
import os
import unittest
from config_regression import cfg
import shutil, tempfile
import vtkbone


class TestDirectMechanics(unittest.TestCase):
    aim_filename = 'test25a.aim'
    n88_filenames = [
        'test25a_strain_xx_solved.n88model', 'test25a_strain_xy_solved.n88model',
        'test25a_strain_yy_solved.n88model', 'test25a_strain_yz_solved.n88model', 
        'test25a_strain_zx_solved.n88model', 'test25a_strain_zz_solved.n88model'
    ]

    def setUp(self):
        # Create temporary directory to work in
        self.test_dir = tempfile.mkdtemp()

        for filename in self.n88_filenames + [self.aim_filename]:
            # Fetch the data
            download_location = cfg['DOWNLOAD_TESTING_DATA'](filename)
            self.assertNotEqual(download_location, '', 'Unable to download file ' + filename)

            # Copy to temporary directory
            shutil.copy(download_location, self.test_dir)
            self.assertTrue(os.path.isfile(os.path.join(self.test_dir, filename)))

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_directmechanics_generate(self):
        '''Can run `n88directmechanics --generate` on a file'''
        # Run command
        command = ['n88directmechanics', '--generate', os.path.join(self.test_dir, self.aim_filename)]
        self.assertTrue(
            cfg['RUN_CALL'](command),
            'Cannot call \"{}\"'.format(' '.join(command))
        )

        # Should create many output files
        output_extensions = [
            '_strain_xx.n88model', '_strain_yy.n88model', '_strain_zz.n88model',
            '_strain_yz.n88model', '_strain_zx.n88model', '_strain_xy.n88model'
        ]

        # File should exist and be readable
        errorObserver = cfg['ERROR_OBSERVER']()
        reader = vtkbone.vtkboneN88ModelReader()
        reader.AddObserver ("ErrorEvent", errorObserver)
        for extension in output_extensions:
            output_filename = os.path.join(self.test_dir, self.aim_filename.replace('.aim', extension))
            self.assertTrue(os.path.isfile(output_filename))

            reader.SetFileName(output_filename)
            reader.Update()
            self.assertFalse(errorObserver.ErrorOccurred(), errorObserver.ErrorMessage())


if __name__ == '__main__':
    unittest.main()
