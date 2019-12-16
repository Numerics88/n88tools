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


class TestDirectMechanicsRegression(unittest.TestCase):
    aim_filename = 'test25a.aim'
    n88_filenames = [
        'test25a_strain_xx_solved.n88model', 'test25a_strain_xy_solved.n88model',
        'test25a_strain_yy_solved.n88model', 'test25a_strain_yz_solved.n88model', 
        'test25a_strain_zx_solved.n88model', 'test25a_strain_zz_solved.n88model'
    ]

    def setUp(self):
        '''Note that we avoid name conflicts with other tests by postpending
            'solved' to the filename. Here, we overwrite the expected file
            with the *solved file so we can run analyze
        '''
        # Create temporary directory to work in
        self.test_dir = tempfile.mkdtemp()

        for filename in self.n88_filenames + [self.aim_filename]:
            # Fetch the data
            download_location = cfg['DOWNLOAD_TESTING_DATA'](filename)
            self.assertNotEqual(download_location, '', 'Unable to download file ' + filename)

            # Copy to temporary directory
            shutil.copy(download_location, self.test_dir)
            self.assertTrue(os.path.isfile(os.path.join(self.test_dir, filename)))

        # Replace solved files
        for filename in self.n88_filenames:
            input_filename = os.path.join(self.test_dir, filename)
            output_filename = os.path.join(self.test_dir, filename.replace('_solved', ''))
            shutil.copyfile(input_filename, output_filename)

        # Run command
        command = ['n88directmechanics', '--analyze', os.path.join(self.test_dir, self.aim_filename)]
        self.output = subprocess.check_output(command)

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_directmechanics_analyze_ran(self):
        '''Can run `n88directmechanics --analyze` on a file'''
        self.assertNotEqual(self.output, '')

    def test_Exx(self):
        '''`n88directmechanics` gives correct Exx in specimen and best orthotropic coordinate system'''
        expression=r'\s+Exx\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 1318.57, 'Exx incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 1983.33, 'Exx incorrect in best orthotropic coordinate system')

    def test_Eyy(self):
        '''`n88directmechanics` gives correct Eyy in specimen and best orthotropic coordinate system'''
        expression=r'\s+Eyy\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 1783.00, 'Eyy incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 1632.04, 'Eyy incorrect in best orthotropic coordinate system')

    def test_Ezz(self):
        '''`n88directmechanics` gives correct Ezz in specimen and best orthotropic coordinate system'''
        expression=r'\s+Ezz\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 1588.87, 'Ezz incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 1243.74, 'Ezz incorrect in best orthotropic coordinate system')

    def test_Gyz(self):
        '''`n88directmechanics` gives correct Gyz in specimen and best orthotropic coordinate system'''
        expression=r'\s+Gyz\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 722.08, 'Gyz incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 653.93, 'Gyz incorrect in best orthotropic coordinate system')

    def test_Gzx(self):
        '''`n88directmechanics` gives correct Gzx in specimen and best orthotropic coordinate system'''
        expression=r'\s+Gzx\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 616.37, 'Gzx incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 680.20, 'Gzx incorrect in best orthotropic coordinate system')

    def test_Gxy(self):
        '''`n88directmechanics` gives correct Gxy in specimen and best orthotropic coordinate system'''
        expression=r'\s+Gxy\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 735.92, 'Gxy incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 703.18, 'Gxy incorrect in best orthotropic coordinate system')

    def test_nu_yx(self):
        '''`n88directmechanics` gives correct nu_yx in specimen and best orthotropic coordinate system'''
        expression=r'\s+nu_yx\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.28394, 'nu_yx incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 0.13908, 'nu_yx incorrect in best orthotropic coordinate system')
    
    def test_nu_zx(self):
        '''`n88directmechanics` gives correct nu_zx in specimen and best orthotropic coordinate system'''
        expression=r'\s+nu_zx\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.27060, 'nu_zx incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 0.17306, 'nu_zx incorrect in best orthotropic coordinate system')

    def test_nu_xy(self):
        '''`n88directmechanics` gives correct nu_xy in specimen and best orthotropic coordinate system'''
        expression=r'\s+nu_xy\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.20998, 'nu_xy incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 0.16902, 'nu_xy incorrect in best orthotropic coordinate system')

    def test_nu_zy(self):
        '''`n88directmechanics` gives correct nu_zy in specimen and best orthotropic coordinate system'''
        expression=r'\s+nu_zy\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.15914, 'nu_zy incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 0.23819, 'nu_zy incorrect in best orthotropic coordinate system')

    def test_nu_xz(self):
        '''`n88directmechanics` gives correct nu_xz in specimen and best orthotropic coordinate system'''
        expression=r'\s+nu_xz\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.22456, 'nu_xz incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 0.27597, 'nu_xz incorrect in best orthotropic coordinate system')

    def test_nu_yz(self):
        '''`n88directmechanics` gives correct nu_yz in specimen and best orthotropic coordinate system'''
        expression=r'\s+nu_yz\s*=\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)'
        result = re.findall(expression, self.output)
        self.assertAlmostEqual(float(result[0][0]), 0.17858, 'nu_yz incorrect in specimen coordinate system')
        self.assertAlmostEqual(float(result[1][0]), 0.31256, 'nu_yz incorrect in best orthotropic coordinate system')

    def test_apparent_stiffness_matrix_in_specimen_coordinate_system(self):
        '''`n88directmechanics` gives correct apparent stiffness matrix in specimen coordinate system'''
        expression=r'(\[\[[0-9\s\.+-\[\]e]*\])'
        result = re.findall(expression, self.output)[0].replace('[','').replace(']', '')
        M = np.fromstring(result, dtype=float, count=6*6, sep=' ').reshape((6,6))

        D = np.array([
            [1571.65,   540.022,  513.776,   -5.1,   -121.256,  -68.126],
            [ 540.022, 2029.057,  470.04,    74.406,  -53.845,  -43.017],
            [ 513.776,  470.04,  1803.981,   10.656,  -57.049,  -20.134],
            [  -5.1,     74.406,   10.656,  730.111,  -32.407,  -49.044],
            [-121.256,  -53.845,  -57.049,  -32.407,  627.399,    5.524],
            [ -68.126,  -43.017,  -20.134,  -49.044,    5.524,  742.336]
        ])

        self.assertTrue(np.allclose(M, D, atol=1e-2), '{} {}'.format(M, D))

    def test_apparent_compliance_matrix_in_specimen_coordinate_system(self):
        '''`n88directmechanics` gives correct apparent compliance matrix in specimen coordinate system'''
        expression=r'(\[\[[0-9\s\.+-\[\]e]*\])'
        result = re.findall(expression, self.output)[1].replace('[','').replace(']', '')
        M = np.fromstring(result, dtype=float, count=6*6, sep=' ').reshape((6,6))

        D = np.array([
            [ 7.584e-04, -1.592e-04, -1.703e-04,  3.311e-05,  1.186e-04,  5.706e-05],
            [-1.592e-04,  5.609e-04, -1.002e-04, -5.580e-05,  5.266e-06,  1.144e-05],
            [-1.703e-04, -1.002e-04,  6.294e-04,  2.312e-07,  1.577e-05, -4.465e-06],
            [ 3.311e-05, -5.580e-05,  2.312e-07,  1.385e-03,  7.237e-05,  9.077e-05],
            [ 1.186e-04,  5.266e-06,  1.577e-05,  7.237e-05,  1.622e-03,  4.328e-06],
            [ 5.706e-05,  1.144e-05, -4.465e-06,  9.077e-05,  4.328e-06,  1.359e-03]
        ])

        self.assertTrue(np.allclose(M, D, atol=1e-3), '{} {}'.format(M, D))

    def test_optimal_rotation_matrix(self):
        '''`n88directmechanics` gives optimal rotation matrix'''
        expression=r'(\[\[[0-9\s\.+-\[\]e]*\])'
        result = re.findall(expression, self.output)[2].replace('[','').replace(']', '')
        M = np.fromstring(result, dtype=float, count=3*3, sep=' ').reshape((3,3))

        D = np.array([
            [-0.369,    0.82505,  0.42795],
            [-0.27759, -0.53726,  0.79643],
            [ 0.88701,  0.17509,  0.42727]
        ])

        self.assertTrue(np.allclose(M, D, atol=1e-2), '{} {}'.format(M, D))

    def test_apparent_stiffness_matrix_in_best_orthotropic_coordinate_system(self):
        '''`n88directmechanics` gives correct apparent stiffness matrix in best orthotropic coordinate system'''
        expression=r'(\[\[[0-9\s\.+-\[\]e]*\])'
        result = re.findall(expression, self.output)[3].replace('[','').replace(']', '')
        M = np.fromstring(result, dtype=float, count=6*6, sep=' ').reshape((6,6))

        D = np.array([
            [2198.279,  460.302,  487.901,  -26.489,  -29.311,    4.707],
            [ 460.302, 1861.202,  521.389,  -36.282,  -10.552,   -2.734],
            [ 487.901,  521.389, 1453.699,    4.375,   24.376,    3.996],
            [ -26.489,  -36.282,    4.375,  657.735,  -38.298,   18.582],
            [ -29.311,  -10.552,   24.376,  -38.298,  683.941,  -12.3  ],
            [   4.707,   -2.734,    3.996,   18.582,  -12.3,    703.922]
        ])

        self.assertTrue(np.allclose(M, D, atol=1e-2), '{} {}'.format(M, D))

    def test_apparent_compliance_matrix_in_best_orthotropic_coordinate_system(self):
        '''`n88directmechanics` gives correct apparent compliance matrix in best orthotropic coordinate system'''
        expression=r'(\[\[[0-9\s\.+-\[\]e]*\])'
        result = re.findall(expression, self.output)[4].replace('[','').replace(']', '')
        M = np.fromstring(result, dtype=float, count=6*6, sep=' ').reshape((6,6))

        D = np.array([
            [ 5.042e-04, -8.522e-05, -1.391e-04,  1.814e-05,  2.622e-05, -2.934e-06],
            [-8.522e-05,  6.127e-04, -1.915e-04,  3.239e-05,  1.450e-05,  3.435e-06],
            [-1.391e-04, -1.915e-04,  8.040e-04, -2.366e-05, -3.898e-05, -4.435e-06],
            [ 1.814e-05,  3.239e-05, -2.366e-05,  1.529e-03,  8.705e-05, -3.871e-05],
            [ 2.622e-05,  1.450e-05, -3.898e-05,  8.705e-05,  1.470e-03,  2.349e-05],
            [-2.934e-06,  3.435e-06, -4.435e-06, -3.871e-05,  2.349e-05,  1.422e-03]
        ])

        self.assertTrue(np.allclose(M, D, atol=1e-3), '{} {}'.format(M, D))


if __name__ == '__main__':
    unittest.main()
