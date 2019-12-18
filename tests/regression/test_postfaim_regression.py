from __future__ import division
import os
import unittest
from n88tools.tables import get_tables, lookups
from config_regression import cfg
import shutil, tempfile
import subprocess


class TestPostFaimRegression(unittest.TestCase):
    '''This test isn't perfect because it doesn't validate the
    entries in the table'''
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
        self.output = subprocess.check_output(command)
        self.assertTrue(os.path.isfile(self.output_file))

        # Split tables
        self.tables = get_tables(self.output_file)

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_postfaim_model_input(self):
        '''`n88postfaim` returns correct model input table'''
        self.assertEqual(lookups['filename'](),             'test25a_uniaxial_solved.n88model')
        self.assertEqual(lookups['dx_els'](),               '0.034')
        self.assertEqual(lookups['dy_els'](),               '0.034')
        self.assertEqual(lookups['dz_els'](),               '0.034')
        self.assertEqual(lookups['num_els'](),              '7087')
        self.assertEqual(lookups['num_nodes'](),            '9938')
        self.assertEqual(lookups['num_nodes_per_els'](),    '8')
        self.assertEqual(lookups['dim_of_problem'](),       '3')

    def test_postfaim_materials(self):
        '''`n88postfaim` returns correct materials table'''
        self.assertEqual(lookups['num_mats'](),             '1')
        self.assertEqual(lookups["id_mat%m"](1),            '127')
        self.assertEqual(lookups["count_mat%m"](1),         '7087')
        
    def test_postfaim_post_processing_sets(self):
        '''`n88postfaim` returns correct post processing sets'''
        self.assertEqual(lookups['num_pp_sets'](),          '2')

    def test_postfaim_strain(self):
        '''`n88postfaim` returns correct strain'''
        self.assertEqual(lookups['eps_ave_xx']('ALL'),      '1.366E-03')
        self.assertEqual(lookups['eps_stddev_xx']('ALL'),   '1.469E-03')
        self.assertEqual(lookups['eps_min_xx']('ALL'),      '-3.082E-03')
        self.assertEqual(lookups['eps_max_xx']('ALL'),      '8.352E-03')
        self.assertEqual(lookups['eps_skew_xx']('ALL'),     '5.199E-01')
        self.assertEqual(lookups['eps_kurt_xx']('ALL'),     '2.617E-01')
        self.assertEqual(lookups['eps_median_xx']('ALL'),   '1.190E-03')
        self.assertEqual(lookups['eps_perc05_xx']('ALL'),   '-6.877E-04')
        self.assertEqual(lookups['eps_perc25_xx']('ALL'),   '2.411E-04')
        self.assertEqual(lookups['eps_perc75_xx']('ALL'),   '2.344E-03')
        self.assertEqual(lookups['eps_perc95_xx']('ALL'),   '3.960E-03')

        self.assertEqual(lookups['eps_ave_yy']('ALL'),      '1.365E-03')
        self.assertEqual(lookups['eps_stddev_yy']('ALL'),   '1.233E-03')
        self.assertEqual(lookups['eps_min_yy']('ALL'),      '-2.343E-03')
        self.assertEqual(lookups['eps_max_yy']('ALL'),      '7.286E-03')
        self.assertEqual(lookups['eps_skew_yy']('ALL'),     '6.196E-01')
        self.assertEqual(lookups['eps_kurt_yy']('ALL'),     '1.411E-01')
        self.assertEqual(lookups['eps_median_yy']('ALL'),   '1.174E-03')
        self.assertEqual(lookups['eps_perc05_yy']('ALL'),   '-2.790E-04')
        self.assertEqual(lookups['eps_perc25_yy']('ALL'),   '4.239E-04')
        self.assertEqual(lookups['eps_perc75_yy']('ALL'),   '2.201E-03')
        self.assertEqual(lookups['eps_perc95_yy']('ALL'),   '3.604E-03')

        self.assertEqual(lookups['eps_ave_zz']('ALL'),      '-4.553E-03')
        self.assertEqual(lookups['eps_stddev_zz']('ALL'),   '3.959E-03')
        self.assertEqual(lookups['eps_min_zz']('ALL'),      '-2.435E-02')
        self.assertEqual(lookups['eps_max_zz']('ALL'),      '5.128E-03')
        self.assertEqual(lookups['eps_skew_zz']('ALL'),     '-7.861E-01')
        self.assertEqual(lookups['eps_kurt_zz']('ALL'),     '2.917E-01')
        self.assertEqual(lookups['eps_median_zz']('ALL'),   '-3.827E-03')
        self.assertEqual(lookups['eps_perc05_zz']('ALL'),   '-1.188E-02')
        self.assertEqual(lookups['eps_perc25_zz']('ALL'),   '-7.131E-03')
        self.assertEqual(lookups['eps_perc75_zz']('ALL'),   '-1.227E-03')
        self.assertEqual(lookups['eps_perc95_zz']('ALL'),   '1.790E-04')

        self.assertEqual(lookups['gam_ave_yz']('ALL'),      '1.619E-06')
        self.assertEqual(lookups['gam_stddev_yz']('ALL'),   '1.967E-03')
        self.assertEqual(lookups['gam_min_yz']('ALL'),      '-1.340E-02')
        self.assertEqual(lookups['gam_max_yz']('ALL'),      '1.403E-02')
        self.assertEqual(lookups['gam_skew_yz']('ALL'),     '-1.347E-01')
        self.assertEqual(lookups['gam_kurt_yz']('ALL'),     '4.668E+00')
        self.assertEqual(lookups['gam_median_yz']('ALL'),   '-7.380E-06')
        self.assertEqual(lookups['gam_perc05_yz']('ALL'),   '-3.305E-03')
        self.assertEqual(lookups['gam_perc25_yz']('ALL'),   '-8.448E-04')
        self.assertEqual(lookups['gam_perc75_yz']('ALL'),   '9.502E-04')
        self.assertEqual(lookups['gam_perc95_yz']('ALL'),   '2.969E-03')

        self.assertEqual(lookups['gam_ave_zx']('ALL'),      '-2.300E-07')
        self.assertEqual(lookups['gam_stddev_zx']('ALL'),   '3.072E-03')
        self.assertEqual(lookups['gam_min_zx']('ALL'),      '-1.665E-02')
        self.assertEqual(lookups['gam_max_zx']('ALL'),      '1.641E-02')
        self.assertEqual(lookups['gam_skew_zx']('ALL'),     '-3.833E-01')
        self.assertEqual(lookups['gam_kurt_zx']('ALL'),     '1.650E+00')
        self.assertEqual(lookups['gam_median_zx']('ALL'),   '1.241E-04')
        self.assertEqual(lookups['gam_perc05_zx']('ALL'),   '-5.541E-03')
        self.assertEqual(lookups['gam_perc25_zx']('ALL'),   '-1.543E-03')
        self.assertEqual(lookups['gam_perc75_zx']('ALL'),   '1.867E-03')
        self.assertEqual(lookups['gam_perc95_zx']('ALL'),   '4.600E-03')

        self.assertEqual(lookups['gam_ave_xy']('ALL'),      '1.691E-06')
        self.assertEqual(lookups['gam_stddev_xy']('ALL'),   '1.349E-03')
        self.assertEqual(lookups['gam_min_xy']('ALL'),      '-8.331E-03')
        self.assertEqual(lookups['gam_max_xy']('ALL'),      '9.056E-03')
        self.assertEqual(lookups['gam_skew_xy']('ALL'),     '1.932E-01')
        self.assertEqual(lookups['gam_kurt_xy']('ALL'),     '2.494E+00')
        self.assertEqual(lookups['gam_median_xy']('ALL'),   '-4.707E-05')
        self.assertEqual(lookups['gam_perc05_xy']('ALL'),   '-2.154E-03')
        self.assertEqual(lookups['gam_perc25_xy']('ALL'),   '-7.068E-04')
        self.assertEqual(lookups['gam_perc75_xy']('ALL'),   '6.621E-04')
        self.assertEqual(lookups['gam_perc95_xy']('ALL'),   '2.289E-03')

    def test_postfaim_stess(self):
        '''`n88postfaim` returns correct stress'''
        self.assertEqual(lookups['sig_ave_xx']('ALL'),      '1.471E-03')
        self.assertEqual(lookups['sig_stddev_xx']('ALL'),   '8.142E+00')
        self.assertEqual(lookups['sig_min_xx']('ALL'),      '-3.315E+01')
        self.assertEqual(lookups['sig_max_xx']('ALL'),      '5.362E+01')
        self.assertEqual(lookups['sig_skew_xx']('ALL'),     '8.357E-01')
        self.assertEqual(lookups['sig_kurt_xx']('ALL'),     '4.373E+00')
        self.assertEqual(lookups['sig_median_xx']('ALL'),   '-2.199E-01')
        self.assertEqual(lookups['sig_perc05_xx']('ALL'),   '-1.274E+01')
        self.assertEqual(lookups['sig_perc25_xx']('ALL'),   '-3.804E+00')
        self.assertEqual(lookups['sig_perc75_xx']('ALL'),   '3.066E+00')
        self.assertEqual(lookups['sig_perc95_xx']('ALL'),   '1.374E+01')

        self.assertEqual(lookups['sig_ave_yy']('ALL'),      '-5.360E-03')
        self.assertEqual(lookups['sig_stddev_yy']('ALL'),   '5.331E+00')
        self.assertEqual(lookups['sig_min_yy']('ALL'),      '-3.269E+01')
        self.assertEqual(lookups['sig_max_yy']('ALL'),      '3.300E+01')
        self.assertEqual(lookups['sig_skew_yy']('ALL'),     '4.997E-02')
        self.assertEqual(lookups['sig_kurt_yy']('ALL'),     '6.425E+00')
        self.assertEqual(lookups['sig_median_yy']('ALL'),   '1.078E-03')
        self.assertEqual(lookups['sig_perc05_yy']('ALL'),   '-7.964E+00')
        self.assertEqual(lookups['sig_perc25_yy']('ALL'),   '-2.122E+00')
        self.assertEqual(lookups['sig_perc75_yy']('ALL'),   '2.016E+00')
        self.assertEqual(lookups['sig_perc95_yy']('ALL'),   '8.143E+00')

        self.assertEqual(lookups['sig_ave_zz']('ALL'),      '-3.110E+01')
        self.assertEqual(lookups['sig_stddev_zz']('ALL'),   '2.815E+01')
        self.assertEqual(lookups['sig_min_zz']('ALL'),      '-1.729E+02')
        self.assertEqual(lookups['sig_max_zz']('ALL'),      '3.485E+01')
        self.assertEqual(lookups['sig_skew_zz']('ALL'),     '-8.286E-01')
        self.assertEqual(lookups['sig_kurt_zz']('ALL'),     '4.442E-01')
        self.assertEqual(lookups['sig_median_zz']('ALL'),   '-2.598E+01')
        self.assertEqual(lookups['sig_perc05_zz']('ALL'),   '-8.272E+01')
        self.assertEqual(lookups['sig_perc25_zz']('ALL'),   '-4.929E+01')
        self.assertEqual(lookups['sig_perc75_zz']('ALL'),   '-7.101E+00')
        self.assertEqual(lookups['sig_perc95_zz']('ALL'),   '1.659E+00')

        self.assertEqual(lookups['sig_ave_yz']('ALL'),      '4.251E-03')
        self.assertEqual(lookups['sig_stddev_yz']('ALL'),   '5.166E+00')
        self.assertEqual(lookups['sig_min_yz']('ALL'),      '-3.519E+01')
        self.assertEqual(lookups['sig_max_yz']('ALL'),      '3.685E+01')
        self.assertEqual(lookups['sig_skew_yz']('ALL'),     '-1.347E-01')
        self.assertEqual(lookups['sig_kurt_yz']('ALL'),     '4.668E+00')
        self.assertEqual(lookups['sig_median_yz']('ALL'),   '-1.938E-02')
        self.assertEqual(lookups['sig_perc05_yz']('ALL'),   '-8.681E+00')
        self.assertEqual(lookups['sig_perc25_yz']('ALL'),   '-2.219E+00')
        self.assertEqual(lookups['sig_perc75_yz']('ALL'),   '2.496E+00')
        self.assertEqual(lookups['sig_perc95_yz']('ALL'),   '7.798E+00')

        self.assertEqual(lookups['sig_ave_zx']('ALL'),      '-6.041E-04')
        self.assertEqual(lookups['sig_stddev_zx']('ALL'),   '8.068E+00')
        self.assertEqual(lookups['sig_min_zx']('ALL'),      '-4.373E+01')
        self.assertEqual(lookups['sig_max_zx']('ALL'),      '4.310E+01')
        self.assertEqual(lookups['sig_skew_zx']('ALL'),     '-3.833E-01')
        self.assertEqual(lookups['sig_kurt_zx']('ALL'),     '1.650E+00')
        self.assertEqual(lookups['sig_median_zx']('ALL'),   '3.260E-01')
        self.assertEqual(lookups['sig_perc05_zx']('ALL'),   '-1.455E+01')
        self.assertEqual(lookups['sig_perc25_zx']('ALL'),   '-4.053E+00')
        self.assertEqual(lookups['sig_perc75_zx']('ALL'),   '4.903E+00')
        self.assertEqual(lookups['sig_perc95_zx']('ALL'),   '1.208E+01')

        self.assertEqual(lookups['sig_ave_xy']('ALL'),      '4.442E-03')
        self.assertEqual(lookups['sig_stddev_xy']('ALL'),   '3.543E+00')
        self.assertEqual(lookups['sig_min_xy']('ALL'),      '-2.188E+01')
        self.assertEqual(lookups['sig_max_xy']('ALL'),      '2.378E+01')
        self.assertEqual(lookups['sig_skew_xy']('ALL'),     '1.932E-01')
        self.assertEqual(lookups['sig_kurt_xy']('ALL'),     '2.494E+00')
        self.assertEqual(lookups['sig_median_xy']('ALL'),   '-1.236E-01')
        self.assertEqual(lookups['sig_perc05_xy']('ALL'),   '-5.659E+00')
        self.assertEqual(lookups['sig_perc25_xy']('ALL'),   '-1.856E+00')
        self.assertEqual(lookups['sig_perc75_xy']('ALL'),   '1.739E+00')
        self.assertEqual(lookups['sig_perc95_xy']('ALL'),   '6.013E+00')

    def test_postfaim_strain_energy_density(self):
        '''`n88postfaim` returns correct strain energy density'''
        self.assertEqual(lookups['sed_avg'](),              '1.494E-01')
        self.assertEqual(lookups['sed_stddev'](),           '1.913E-01')
        self.assertEqual(lookups['sed_min'](),              '1.794E-07')
        self.assertEqual(lookups['sed_max'](),              '2.046E+00')
        self.assertEqual(lookups['sed_skew'](),             '2.560E+00')
        self.assertEqual(lookups['sed_kurt'](),             '1.070E+01')
        self.assertEqual(lookups['sed_median'](),           '7.737E-02')
        self.assertEqual(lookups['sed_perc05'](),           '1.827E-03')
        self.assertEqual(lookups['sed_perc25'](),           '1.973E-02')
        self.assertEqual(lookups['sed_perc75'](),           '2.092E-01')
        self.assertEqual(lookups['sed_perc95'](),           '5.287E-01')

    def test_postfaim_von_mises_stress(self):
        '''`n88postfaim` returns correct von mises stress'''
        self.assertEqual(lookups['svm_avg'](),              '3.724E+01')
        self.assertEqual(lookups['svm_stddev'](),           '2.527E+01')
        self.assertEqual(lookups['svm_min'](),              '5.133E-02')
        self.assertEqual(lookups['svm_max'](),              '1.651E+02')
        self.assertEqual(lookups['svm_skew'](),             '7.878E-01')
        self.assertEqual(lookups['svm_kurt'](),             '2.688E-01')
        self.assertEqual(lookups['svm_median'](),           '3.287E+01')
        self.assertEqual(lookups['svm_perc05'](),           '5.064E+00')
        self.assertEqual(lookups['svm_perc25'](),           '1.654E+01')
        self.assertEqual(lookups['svm_perc75'](),           '5.384E+01')
        self.assertEqual(lookups['svm_perc95'](),           '8.434E+01')

    def test_postfaim_nodal_displacements(self):
        '''`n88postfaim` returns correct nodal displacements'''
        n=1
        self.assertEqual(lookups['dx_avg_ns%n'](n),         '-6.533E-04')
        self.assertEqual(lookups['dy_avg_ns%n'](n),         '2.015E-04')
        self.assertEqual(lookups['dz_avg_ns%n'](n),         '-8.500E-03')
        self.assertEqual(lookups['dx_stddev_ns%n'](n),      '5.212E-03')
        self.assertEqual(lookups['dy_stddev_ns%n'](n),      '2.200E-03')
        self.assertEqual(lookups['dz_stddev_ns%n'](n),      '0.000E+00')
        self.assertEqual(lookups['dx_min_ns%n'](n),         '-8.831E-03')
        self.assertEqual(lookups['dy_min_ns%n'](n),         '-3.293E-03')
        self.assertEqual(lookups['dz_min_ns%n'](n),         '-8.500E-03')
        self.assertEqual(lookups['dx_max_ns%n'](n),         '1.222E-02')
        self.assertEqual(lookups['dy_max_ns%n'](n),         '5.999E-03')
        self.assertEqual(lookups['dz_max_ns%n'](n),         '-8.500E-03')
        self.assertEqual(lookups['dx_median_ns%n'](n),      '-9.100E-04')
        self.assertEqual(lookups['dy_median_ns%n'](n),      '-2.456E-04')
        self.assertEqual(lookups['dz_median_ns%n'](n),      '-8.500E-03')

        n=2
        self.assertEqual(lookups['dx_avg_ns%n'](n),         '-3.620E-05')
        self.assertEqual(lookups['dy_avg_ns%n'](n),         '-1.379E-04')
        self.assertEqual(lookups['dz_avg_ns%n'](n),         '0.000E+00')
        self.assertEqual(lookups['dx_stddev_ns%n'](n),      '1.418E-03')
        self.assertEqual(lookups['dy_stddev_ns%n'](n),      '4.155E-04')
        self.assertEqual(lookups['dz_stddev_ns%n'](n),      '0.000E+00')
        self.assertEqual(lookups['dx_min_ns%n'](n),         '-3.627E-03')
        self.assertEqual(lookups['dy_min_ns%n'](n),         '-9.885E-04')
        self.assertEqual(lookups['dz_min_ns%n'](n),         '0.000E+00')
        self.assertEqual(lookups['dx_max_ns%n'](n),         '2.145E-03')
        self.assertEqual(lookups['dy_max_ns%n'](n),         '1.444E-03')
        self.assertEqual(lookups['dz_max_ns%n'](n),         '0.000E+00')
        self.assertEqual(lookups['dx_median_ns%n'](n),      '1.639E-04')
        self.assertEqual(lookups['dy_median_ns%n'](n),      '-2.212E-04')
        self.assertEqual(lookups['dz_median_ns%n'](n),      '0.000E+00')

    def test_postfaim_nodal_forces(self):
        '''`n88postfaim` returns correct nodal forces'''
        n=1
        self.assertEqual(lookups['fx_ns%n'](n),             '-2.689E-04')
        self.assertEqual(lookups['fy_ns%n'](n),             '1.411E-04')
        self.assertEqual(lookups['fz_ns%n'](n),             '-1.019E+01')
        self.assertEqual(lookups['fx_avg_ns%n'](n),         '-9.671E-07')
        self.assertEqual(lookups['fy_avg_ns%n'](n),         '5.076E-07')
        self.assertEqual(lookups['fz_avg_ns%n'](n),         '-3.665E-02')
        self.assertEqual(lookups['fx_stddev_ns%n'](n),      '9.938E-06')
        self.assertEqual(lookups['fy_stddev_ns%n'](n),      '9.801E-06')
        self.assertEqual(lookups['fz_stddev_ns%n'](n),      '4.090E-02')
        self.assertEqual(lookups['fx_min_ns%n'](n),         '-2.986E-05')
        self.assertEqual(lookups['fy_min_ns%n'](n),         '-2.690E-05')
        self.assertEqual(lookups['fz_min_ns%n'](n),         '-1.450E-01')
        self.assertEqual(lookups['fx_max_ns%n'](n),         '3.551E-05')
        self.assertEqual(lookups['fy_max_ns%n'](n),         '3.223E-05')
        self.assertEqual(lookups['fz_max_ns%n'](n),         '1.856E-02')
        self.assertEqual(lookups['fx_median_ns%n'](n),      '-8.065E-07')
        self.assertEqual(lookups['fy_median_ns%n'](n),      '-2.757E-07')
        self.assertEqual(lookups['fz_median_ns%n'](n),      '-2.530E-02')

        n=2
        self.assertEqual(lookups['fx_ns%n'](n),             '-2.041E-04')
        self.assertEqual(lookups['fy_ns%n'](n),             '-4.824E-04')
        self.assertEqual(lookups['fz_ns%n'](n),             '1.019E+01')
        self.assertEqual(lookups['fx_avg_ns%n'](n),         '-5.076E-07')
        self.assertEqual(lookups['fy_avg_ns%n'](n),         '-1.200E-06')
        self.assertEqual(lookups['fz_avg_ns%n'](n),         '2.535E-02')
        self.assertEqual(lookups['fx_stddev_ns%n'](n),      '6.390E-06')
        self.assertEqual(lookups['fy_stddev_ns%n'](n),      '6.175E-06')
        self.assertEqual(lookups['fz_stddev_ns%n'](n),      '2.111E-02')
        self.assertEqual(lookups['fx_min_ns%n'](n),         '-1.944E-05')
        self.assertEqual(lookups['fy_min_ns%n'](n),         '-2.133E-05')
        self.assertEqual(lookups['fz_min_ns%n'](n),         '-7.386E-03')
        self.assertEqual(lookups['fx_max_ns%n'](n),         '2.142E-05')
        self.assertEqual(lookups['fy_max_ns%n'](n),         '1.524E-05')
        self.assertEqual(lookups['fz_max_ns%n'](n),         '8.560E-02')
        self.assertEqual(lookups['fx_median_ns%n'](n),      '-5.811E-07')
        self.assertEqual(lookups['fy_median_ns%n'](n),      '-7.346E-07')
        self.assertEqual(lookups['fz_median_ns%n'](n),      '2.330E-02')

    def test_postfaim_load_sharing(self):
        '''`n88postfaim` returns correct load sharing'''
        m=1

        n=1
        self.assertEqual(lookups['fx_ns%n_mat%m'](n,m),     '-2.6886E-04')
        self.assertEqual(lookups['fy_ns%n_mat%m'](n,m),     '1.4111E-04')
        self.assertEqual(lookups['fz_ns%n_mat%m'](n,m),     '-1.0190E+01')

        n=2
        self.assertEqual(lookups['fx_ns%n_mat%m'](n,m),     '-2.0406E-04')
        self.assertEqual(lookups['fy_ns%n_mat%m'](n,m),     '-4.8235E-04')
        self.assertEqual(lookups['fz_ns%n_mat%m'](n,m),     '1.0191E+01')
        

if __name__ == '__main__':
    unittest.main()
