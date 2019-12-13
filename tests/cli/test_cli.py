from __future__ import division
import subprocess
import unittest

class TestCLI(unittest.TestCase):

    def callable(self, command):
        try:
            subprocess.check_output(command)
        except subprocess.CalledProcessError as e:
            return False
        except OSError as e:
            return False
        return True

    def test_cli(self):
        '''Can call help of all functions'''
        commands = [
             'faim'
            ,'n88coarsen'
            ,'n88compress'
            ,'n88copymodel'
            ,'n88directmechanics'
            ,'n88evaluate'
            ,'n88extractfields'
            ,'n88extractsets'
            ,'n88interpolatesolution'
            ,'n88modelinfo'
            ,'n88pistoia'
            ,'n88postfaim'
            ,'n88tabulate'
        ]
        arguments = [
             ['--help']
            ,['-h']
        ]

        # Going to make calls of the form:
        #   n88coarsen --help
        for command in commands:
            for argument in arguments:
                call = [command] + argument
                self.assertTrue(self.callable(call),
                    'Cannot call \"{}\" from command line'.format(
                        ' '.join(call)
                    ))


if __name__ == '__main__':
    unittest.main()
