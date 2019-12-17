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


class TestModelInfoRegression(unittest.TestCase):
    filename = 'test25a_uniaxial_solved.n88model'

    def setUp(self):
        # Create temporary directory to work in
        self.test_dir = tempfile.mkdtemp()

        # Fetch the data
        input_model = cfg['DOWNLOAD_TESTING_DATA'](self.filename)
        self.assertNotEqual(input_model, '', 'Unable to download file ' + self.filename)

        # Copy to temporary directory
        shutil.copy(input_model, self.test_dir)
        self.model_filename = os.path.join(self.test_dir, self.filename)
        self.assertTrue(os.path.isfile(self.model_filename))

        # Run command
        command = ['n88modelinfo', self.model_filename]
        self.output = subprocess.check_output(command)

        # Split tables
        pattern_block = '----------------------------------------------------------------------'
        split = self.output.split(pattern_block)
        self.assertTrue(len(split)%2 == 1)
        self.tables = {}
        for i in range(len(split)//2):
            key = split[2*i].strip().replace('\n', '').replace('\r', '')
            value = split[2*i+1].replace('\r', '')
            self.tables[key] = value

    def tearDown(self):
        # Remove temporary directory and all files
        shutil.rmtree(self.test_dir)

    def test_modelinfo(self):
        '''Ran `n88modelinfo` successfully'''
        self.assertNotEqual(self.output, '')

    def test_modelinfo_history(self):
        '''History output correct'''
        key = 'History:'
        expected = '''
2019-Dec-12 14:05:47 Model created by n88modelgenerator version 8.0-alpha4
2019-Dec-12 14:06:51 Solved by n88solver_slt 8.0-alpha4
2019-Dec-17 15:04:06 Processed by n88derivedfields 8.0-alpha4
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_log(self):
        '''Log output correct'''
        key = 'Log:'
        expected = '''
2019-Dec-12 14:05:47
n88modelgenerator Version 8.0-alpha4
Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
Licensed to Bryce Besler, laptop; lic. no. 167

input_file   = test25a.aim
output_file  = test25a_uniaxial.n88model
connectivity_filter           = on
test                          = uniaxial
test_axis                     = z
normal_strain                 = -0.01
pin                           = off
top_surface                   = intersection
bottom_surface                = intersection
material_table                = homogeneous
youngs_modulus                = 6829
poissons_ratio                = 0.3

   0.00 Reading input data.
   0.00 Read 17576 points.
   0.00 Image bounds:
        6.6300 7.4800 7.2080 8.0580 1.7000 2.5500
   0.00 Applying connectivity filter.
   0.00 Masked out 0 unconnected voxels.
   0.00 Converting to hexahedral mesh.
   0.01 Generated 7087 hexahedrons.
   0.01 Constructing material table.
   0.01 Material table has 1 entry.
   0.01 Constructing finite element model.
   0.03 Model bounds:
        6.6300 7.4800 7.2080 8.0580 1.7000 2.5500
   0.03 Generated the following constraints:
          bottom_fixed : 402 nodes
          top_displacement : 278 nodes

2019-Dec-12 14:06:51
n88solver_slt version 8.0-alpha4
Copyright (c) 2010-2015, Numerics88 Solutions Ltd.
Licensed to Bryce Besler, laptop; lic. no. 167
Problem:
  active solution       = (none)
  active problem        = Problem1
  number of elements    = 7087
  number of nodes       = 9938
Solver engine:
  engine                = mt
  precision             = double
  threads               = 4
Convergence parameters:
  convergence measure   = set
  convergence tolerance = 1e-06
  convergence window    = 16
  maximum iterations    = 30000
Convergence measure tolerance reached.
Number of linear iterations = 382
Peak data allocation of supervisory thread: 1.08 MiB
Peak data allocation, sum over worker threads : 1.99 MiB

2019-Dec-17 15:04:06
n88derivedfields version 8.0-alpha4
Copyright (c) 2010-2015, Numerics88 Solutions Ltd.
Licensed to Bryce Besler, laptop; lic. no. 167
Problem:
  active solution       = Solution1
  active problem        = Problem1
  number of elements    = 7087
  number of nodes       = 9938
Solver engine:
  precision             = single

'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_active_settings(self):
        '''Active settings output correct'''
        key = 'Active Settings:'
        expected = '''
  Active Solution : Solution1
  Active Problem : Problem1
  Active Part : Part1
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_materials(self):
        '''Materials output correct'''
        key = 'Materials:'
        expected = '''

  Name : NewMaterial1
  Type : LinearIsotropic
  E : 6829.0
  nu : 0.3
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_parts(self):
        '''Parts settings output correct'''
        key = 'Parts:'
        expected = '''

  Name : Part1
  NumberOfNodes : 9938
  Hexahedrons :
    NumberOfNodesPerElement : 8
    NumberOfElements : 7087
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_constraints(self):
        '''Constraints output correct'''
        key = 'Constraints:'
        expected = '''

  Name : bottom_fixed
  Part : Part1
  Type : NodeAxisDisplacement
  NumberOfValues : 402

  Name : top_displacement
  Part : Part1
  Type : NodeAxisDisplacement
  NumberOfValues : 278

  Name : convergence_set
  Part : Part1
  Type : NodeAxisForce
  NumberOfValues : 278
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_node_sets(self):
        '''NodeSets output correct'''
        key = 'NodeSets:'
        expected = '''

  Name : face_z0
  Part : Part1
  NumberOfNodes : 402

  Name : face_z1
  Part : Part1
  NumberOfNodes : 278

  Name : face_x0
  Part : Part1
  NumberOfNodes : 333

  Name : face_x1
  Part : Part1
  NumberOfNodes : 312

  Name : face_y0
  Part : Part1
  NumberOfNodes : 409

  Name : face_y1
  Part : Part1
  NumberOfNodes : 401
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_element_sets(self):
        '''ElementSets output correct'''
        key = 'ElementSets:'
        expected = '''

  Name : face_z0
  Part : Part1
  NumberOfElements : 305

  Name : face_z1
  Part : Part1
  NumberOfElements : 207

  Name : face_x0
  Part : Part1
  NumberOfElements : 281

  Name : face_x1
  Part : Part1
  NumberOfElements : 239

  Name : face_y0
  Part : Part1
  NumberOfElements : 328

  Name : face_y1
  Part : Part1
  NumberOfElements : 317
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_problems(self):
        '''Problems output correct'''
        key = 'Problems:'
        expected = '''

  Name : Problem1
  Part : Part1
  Constraints : bottom_fixed,top_displacement
  ConvergenceSet : convergence_set
  PostProcessingNodeSets : face_z1,face_z0
  PostProcessingElementSets : face_z1,face_z0
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

    def test_modelinfo_solutions(self):
        '''Solutions output correct'''
        key = 'Solutions:'
        expected = '''

  Name : Solution1
  Problem : Problem1
  Variables defined on nodes:
    Displacement
    ReactionForce
  Variables defined on elements:
    Strain
    Stress
    StrainEnergyDensity
    VonMisesStress
'''
        
        self.assertTrue(key in self.tables)
        history = self.tables[key]
        self.assertEqual(history, expected)

if __name__ == '__main__':
    unittest.main()
