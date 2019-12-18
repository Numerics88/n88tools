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
        expected='''filename	dx_els	dy_els	dz_els	num_els	num_nodes	num_nodes_per_els	dim_of_problem	num_mats	num_pp_sets	eps_ave_xx	eps_stddev_xx	eps_min_xx	eps_max_xx	eps_skew_xx	eps_kurt_xx	eps_median_xx	eps_perc05_xx	eps_perc25_xx	eps_perc75_xx	eps_perc95_xx	eps_ave_yy	eps_stddev_yy	eps_min_yy	eps_max_yy	eps_skew_yy	eps_kurt_yy	eps_median_yy	eps_perc05_yy	eps_perc25_yy	eps_perc75_yy	eps_perc95_yy	eps_ave_zz	eps_stddev_zz	eps_min_zz	eps_max_zz	eps_skew_zz	eps_kurt_zz	eps_median_zz	eps_perc05_zz	eps_perc25_zz	eps_perc75_zz	eps_perc95_zz	gam_ave_yz	gam_stddev_yz	gam_min_yz	gam_max_yz	gam_skew_yz	gam_kurt_yz	gam_median_yz	gam_perc05_yz	gam_perc25_yz	gam_perc75_yz	gam_perc95_yz	gam_ave_zx	gam_stddev_zx	gam_min_zx	gam_max_zx	gam_skew_zx	gam_kurt_zx	gam_median_zx	gam_perc05_zx	gam_perc25_zx	gam_perc75_zx	gam_perc95_zx	gam_ave_xy	gam_stddev_xy	gam_min_xy	gam_max_xy	gam_skew_xy	gam_kurt_xy	gam_median_xy	gam_perc05_xy	gam_perc25_xy	gam_perc75_xy	gam_perc95_xy	sig_ave_xx	sig_stddev_xx	sig_min_xx	sig_max_xx	sig_skew_xx	sig_kurt_xx	sig_median_xx	sig_perc05_xx	sig_perc25_xx	sig_perc75_xx	sig_perc95_xx	sig_ave_yy	sig_stddev_yy	sig_min_yy	sig_max_yy	sig_skew_yy	sig_kurt_yy	sig_median_yy	sig_perc05_yy	sig_perc25_yy	sig_perc75_yy	sig_perc95_yy	sig_ave_zz	sig_stddev_zz	sig_min_zz	sig_max_zz	sig_skew_zz	sig_kurt_zz	sig_median_zz	sig_perc05_zz	sig_perc25_zz	sig_perc75_zz	sig_perc95_zz	sig_ave_yz	sig_stddev_yz	sig_min_yz	sig_max_yz	sig_skew_yz	sig_kurt_yz	sig_median_yz	sig_perc05_yz	sig_perc25_yz	sig_perc75_yz	sig_perc95_yz	sig_ave_zx	sig_stddev_zx	sig_min_zx	sig_max_zx	sig_skew_zx	sig_kurt_zx	sig_median_zx	sig_perc05_zx	sig_perc25_zx	sig_perc75_zx	sig_perc95_zx	sig_ave_xy	sig_stddev_xy	sig_min_xy	sig_max_xy	sig_skew_xy	sig_kurt_xy	sig_median_xy	sig_perc05_xy	sig_perc25_xy	sig_perc75_xy	sig_perc95_xy	sed_avg	sed_min	sed_max	sed_stddev	sed_skew	sed_kurt	sed_median	sed_perc05	sed_perc25	sed_perc75	sed_perc95	sed_perc05_mat%m	sed_perc25_mat%m	sed_perc75_mat%m	sed_perc95_mat%m	svm_avg	svm_min	svm_max	svm_stddev	svm_skew	svm_kurt	svm_median	svm_perc05	svm_perc25	svm_perc75	svm_perc95	svm_perc05_mat%m	svm_perc25_mat%m	svm_perc75_mat%m	svm_perc95_mat%m	dx_avg_ns1	dx_avg_ns2	dx_stddev_ns1	dx_stddev_ns2	dx_min_ns1	dx_min_ns2	dx_max_ns1	dx_max_ns2	dx_median_ns1	dx_median_ns2	dy_avg_ns1	dy_avg_ns2	dy_stddev_ns1	dy_stddev_ns2	dy_min_ns1	dy_min_ns2	dy_max_ns1	dy_max_ns2	dy_median_ns1	dy_median_ns2	dz_avg_ns1	dz_avg_ns2	dz_stddev_ns1	dz_stddev_ns2	dz_min_ns1	dz_min_ns2	dz_max_ns1	dz_max_ns2	dz_median_ns1	dz_median_ns2	fx_ns1	fx_ns2	fx_avg_ns1	fx_avg_ns2	fx_stddev_ns1	fx_stddev_ns2	fx_min_ns1	fx_min_ns2	fx_max_ns1	fx_max_ns2	fx_median_ns1	fx_median_ns2	fy_ns1	fy_ns2	fy_avg_ns1	fy_avg_ns2	fy_stddev_ns1	fy_stddev_ns2	fy_min_ns1	fy_min_ns2	fy_max_ns1	fy_max_ns2	fy_median_ns1	fy_median_ns2	fz_ns1	fz_ns2	fz_avg_ns1	fz_avg_ns2	fz_stddev_ns1	fz_stddev_ns2	fz_min_ns1	fz_min_ns2	fz_max_ns1	fz_max_ns2	fz_median_ns1	fz_median_ns2
test25a_uniaxial_solved.n88model	0.034	0.034	0.034	7087	9938	8	3	1	2	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	1.494E-01	1.794E-07	2.046E+00	1.913E-01	2.560E+00	1.070E+01	7.737E-02	1.827E-03	1.973E-02	2.092E-01	5.287E-01	-	-	-	-	3.724E+01	5.133E-02	1.651E+02	2.527E+01	7.878E-01	2.688E-01	3.287E+01	5.064E+00	1.654E+01	5.384E+01	8.434E+01	-	-	-	-	-6.533E-04	-3.620E-05	5.212E-03	1.418E-03	-8.831E-03	-3.627E-03	1.222E-02	2.145E-03	-9.100E-04	1.639E-04	2.015E-04	-1.379E-04	2.200E-03	4.155E-04	-3.293E-03	-9.885E-04	5.999E-03	1.444E-03	-2.456E-04	-2.212E-04	-8.500E-03	0.000E+00	0.000E+00	0.000E+00	-8.500E-03	0.000E+00	-8.500E-03	0.000E+00	-8.500E-03	0.000E+00	-2.689E-04	-2.041E-04	-9.671E-07	-5.076E-07	9.938E-06	6.390E-06	-2.986E-05	-1.944E-05	3.551E-05	2.142E-05	-8.065E-07	-5.811E-07	1.411E-04	-4.824E-04	5.076E-07	-1.200E-06	9.801E-06	6.175E-06	-2.690E-05	-2.133E-05	3.223E-05	1.524E-05	-2.757E-07	-7.346E-07	-1.019E+01	1.019E+01	-3.665E-02	2.535E-02	4.090E-02	2.111E-02	-1.450E-01	-7.386E-03	1.856E-02	8.560E-02	-2.530E-02	2.330E-02
'''

        with open(self.csv_file, 'r') as fp:
            result = fp.read()

        self.assertEqual(expected, result)
        

if __name__ == '__main__':
    unittest.main()
