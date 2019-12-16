from __future__ import division
import os
import unittest
from config_regression import cfg

class TestRegressionData(unittest.TestCase):

    def test_configuration(self):
        '''The configuration file should define the correct variables'''
        variables = [
             'REGRESSION_DATA_DIRECTORY'
            ,'REGRESSION_FILES'
            ,'REGRESSION_DATA_URL'
            ,'DOWNLOAD_TESTING_DATA'
            ,'ERROR_OBSERVER'
            ,'RUN_CALL'
        ]
        for variable in variables:
            self.assertTrue(variable in cfg, 'Could not find \"' + variable + '\" in configuration file')

    def test_data_exists(self):
        '''All regression data expected was downloaded'''
        for filename in cfg['REGRESSION_FILES']:
            # Download the data
            result = cfg['DOWNLOAD_TESTING_DATA'](filename)
            self.assertNotEqual(result, '', 'Unable to download file ' + filename)

            # Be sure it actually downloaded
            self.assertTrue(os.path.isfile(os.path.join(cfg['REGRESSION_DATA_DIRECTORY'], filename)), 'Cannot locate file ' + filename)


if __name__ == '__main__':
    unittest.main()
