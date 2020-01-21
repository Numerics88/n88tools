from __future__ import division
import os
import unittest
from .config_regression import cfg
import shutil, tempfile


class TestTabulateRegression(unittest.TestCase):
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

        # Run postfaim
        self.output_file = os.path.join(self.test_dir, 'postfaim_output.txt')
        command = ['n88postfaim', '--output_file', self.output_file, os.path.join(self.test_dir, self.filename)]
        self.assertTrue(cfg['RUN_CALL'](command))
        self.assertTrue(os.path.isfile(self.output_file))

        # Run tabulate
        self.csv_file = os.path.join(self.test_dir, 'tabulate.csv')
        command = ['n88tabulate', '-H', '--output_file', self.csv_file, self.output_file]
        self.assertTrue(cfg['RUN_CALL'](command))
        self.assertTrue(os.path.isfile(self.csv_file))

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_tabulate_(self):
        '''`n88tabulate` returns correct csv file'''
        expected='filename\tnum_els\tnum_nodes\tnum_mats\tnum_pp_sets\tsed_avg\tsed_min\tsed_max\tsed_stddev\tsed_skew\tsed_kurt\tsed_median\tsvm_avg\tsvm_min\tsvm_max\tsvm_stddev\tsvm_skew\tsvm_kurt\tsvm_median\tdx_avg_ns1\tdx_avg_ns2\tdx_stddev_ns1\tdx_stddev_ns2\tdx_min_ns1\tdx_min_ns2\tdx_max_ns1\tdx_max_ns2\tdx_median_ns1\tdx_median_ns2\tdy_avg_ns1\tdy_avg_ns2\tdy_stddev_ns1\tdy_stddev_ns2\tdy_min_ns1\tdy_min_ns2\tdy_max_ns1\tdy_max_ns2\tdy_median_ns1\tdy_median_ns2\tdz_avg_ns1\tdz_avg_ns2\tdz_stddev_ns1\tdz_stddev_ns2\tdz_min_ns1\tdz_min_ns2\tdz_max_ns1\tdz_max_ns2\tdz_median_ns1\tdz_median_ns2\tfx_ns1\tfx_ns2\tfx_stddev_ns1\tfx_stddev_ns2\tfx_min_ns1\tfx_min_ns2\tfx_max_ns1\tfx_max_ns2\tfx_median_ns1\tfx_median_ns2\tfy_ns1\tfy_ns2\tfy_stddev_ns1\tfy_stddev_ns2\tfy_min_ns1\tfy_min_ns2\tfy_max_ns1\tfy_max_ns2\tfy_median_ns1\tfy_median_ns2\tfz_ns1\tfz_ns2\tfz_stddev_ns1\tfz_stddev_ns2\tfz_min_ns1\tfz_min_ns2\tfz_max_ns1\tfz_max_ns2\tfz_median_ns1\tfz_median_ns2\tfx_avg_ns1\tfx_avg_ns2\tfy_avg_ns1\tfy_avg_ns2\tfz_avg_ns1\tfz_avg_ns2\ntest25a_uniaxial_solved.n88model\t7087\t9938\t1\t2\t1.494E-01\t1.794E-07\t2.046E+00\t1.913E-01\t2.560E+00\t1.070E+01\t7.737E-02\t3.724E+01\t5.133E-02\t1.651E+02\t2.527E+01\t7.878E-01\t2.688E-01\t3.287E+01\t-6.533E-04\t-3.620E-05\t5.212E-03\t1.418E-03\t-8.831E-03\t-3.627E-03\t1.222E-02\t2.145E-03\t-9.100E-04\t1.639E-04\t2.015E-04\t-1.379E-04\t2.200E-03\t4.155E-04\t-3.293E-03\t-9.885E-04\t5.999E-03\t1.444E-03\t-2.456E-04\t-2.212E-04\t-8.500E-03\t0.000E+00\t0.000E+00\t0.000E+00\t-8.500E-03\t0.000E+00\t-8.500E-03\t0.000E+00\t-8.500E-03\t0.000E+00\t-2.689E-04\t-2.041E-04\t9.938E-06\t6.390E-06\t-2.986E-05\t-1.944E-05\t3.551E-05\t2.142E-05\t-8.065E-07\t-5.811E-07\t1.411E-04\t-4.824E-04\t9.801E-06\t6.175E-06\t-2.690E-05\t-2.133E-05\t3.223E-05\t1.524E-05\t-2.757E-07\t-7.346E-07\t-1.019E+01\t1.019E+01\t4.090E-02\t2.111E-02\t-1.450E-01\t-7.386E-03\t1.856E-02\t8.560E-02\t-2.530E-02\t2.330E-02\t-9.671E-07\t-5.076E-07\t5.076E-07\t-1.200E-06\t-3.665E-02\t2.535E-02\n'

        with open(self.csv_file, 'r') as fp:
            result = fp.read()

        self.assertEqual(expected, result)
        

if __name__ == '__main__':
    unittest.main()
