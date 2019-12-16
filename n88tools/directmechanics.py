"""
directmechanics.py

A tool to perform direct mechanics calculations using n88solver to evaluate
mechanics deformation.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from N88ReportedError import N88ReportedError
import os
import argparse
import StringIO
import ConfigParser
from math import *
import numpy
from numpy.core import *
from numpy import linalg
import scipy
from scipy.optimize import fmin_powell as minimize
import transformations
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import vtkbone
import time
import subprocess

try:
    import pkg_resources
    n88tools_version = pkg_resources.require("n88tools")[0].version
except:
    n88tools_version = "unknown"


# -------------------------------------------------------------------------
#  Utility functions

start_time = time.time()
log_string = ""

def log (msg, append=False, timestamp=True):
    """Print message with time stamp.

    Only the first line is printed with a time stamp; subsequent lines
    are indented but not timestamped. If append=True, the first line is
    also printed indented without a time stamp. If timestamp=False, neither
    a time stamp nor indentation is added to any line.
    
    If msg does not end with a line return, one is added.
    """
    global log_string
    if not timestamp:
       lines = msg.splitlines()
       result = ""
       for l in lines:
           result += l + "\n"
    else:
       lines = msg.splitlines()
       result = ""
       if not append:
           result += "%8.2f %s\n" % (time.time()-start_time, lines[0])
       else:
           result += " " * 9 + lines[0] + "\n"
       for l in lines[1:]:
           result +=  " " * 9 + l + "\n"
    print result, 
    log_string += result


class ErrorObserver:

   def __init__(self):
       self.__ErrorOccurred = False
       self.__ErrorMessage = None
       self.CallDataType = 'string0'

   def __call__(self, obj, event, message):
       self.__ErrorOccurred = True
       self.__ErrorMessage = message

   def ErrorOccurred(self):
       occ = self.__ErrorOccurred
       self.__ErrorOccurred = False
       return occ

   def ErrorMessage(self):
       return self.__ErrorMessage


def directmechanics():

    # -------------------------------------------------------------------------
    # Argument parsing and configuration
    #
    # This is complicated by the fact that we want to support both command-line
    # arguments, and a configuration file, with the former taking precedence.
    # To so this, we:
    #   1. Create an instance of argparse that only checks for the --conf option
    #   2. Use ConfigParser to read the conf file if necessary.
    #   3. Create another instance of argparse to actually parse the command line,
    #      populating its default values with those read from ConfigParser.

    conf_parser = argparse.ArgumentParser(
        # Turn off help, so we print all options in response to -h
        add_help=False
        )
    conf_parser.add_argument("-c", "--config",
        help="""Specify a configuration file. The configuration file may specify
    any arguments that take a value, one per line, in the format name=value (leave
    the double dash "--" off of the argument name).""",
        metavar="FILE")
    args, remaining_args = conf_parser.parse_known_args()

    defaults = {
        "material_table": "homogeneous",
        "youngs_modulus": 6829.0,
        "poissons_ratio": 0.3,
        "homminga_maximum_material_id": 127,
        "homminga_modulus_exponent": 1.7,
        "connectivity_filter": "on"
        }

    # Read config file if requested
    if args.config:
        # Some slightly fancy foot-work to get ConfigParser to read a file
        # without sections: we add a header line to create section [root]
        config_str = '[root]\n' + open(args.config, 'r').read()
        config_fp = StringIO.StringIO(config_str)
        config = ConfigParser.RawConfigParser()
        config.readfp(config_fp)
        for key,value in dict(config.items("root")).items():
            defaults[key] = value

    parser = argparse.ArgumentParser(
        prog="n88directmechanics",
        parents=[conf_parser],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description= \
    """Perform Direct Mechanics calculations. If no actions are specified, the
entire generate, solve and analyze sequence will be performed.

Supported input formats:
  DICOM (a directory)
  Scanco AIM (.aim)
  MetaImage (.mha or .mhd)
  VTK XML ImageData (.vti)
""",
        epilog= \
"""Tip: There are no options passed to the solver. If this is required, run
n88directmechanics with the --generate action, then manually run n88solver_slt
on each of the resulting n88model files, using the desired solver options,
then run n88directmechanics again with the --analyze action."""
        )
    parser.set_defaults(**defaults)

    # input arguments
    inputs_group = parser.add_argument_group('input arguments')
    inputs_group.add_argument ("input_file",
        help="""An image file with segmented data.  Note that
    the file name of the original image file should be used even when the only
    action is --solve and/or --analyze; the n88model file names will
    be derived from the image file name.""")


    # Actions
    action_group = parser.add_argument_group('action arguments')
    action_group.add_argument ("--generate", "-g", action='store_true',
        help="Generate models. (Generates 6 n88model files.) May be used together with --solver and/or --analyze.")

    action_group.add_argument ("--solve", "-s", action='store_true',
        help="Solve models. (n88model files are updated with solutions.) May be used together with --generate and/or --analyze.")

    action_group.add_argument ("--analyze", "-a", action='store_true',
        help="Perform direct mechanics analysis on solved files. May be used together with --generate and/or --solve.")

    # Material specification group
    matspec_group = parser.add_argument_group('material specification arguments')
    matspec_group.add_argument ("--material_table",
        choices=["homogeneous", "homminga"],
        help="""Specify method for generating the material table.
    Values:
      homogeneous: Homogeneous material for
        all non-zero segmentation indices.
      homminga: Homminga et al. (2001) J 
        Biomech 34(4):513-517.
            E = E_max * (i/i_max)^exponent
        For the Homminga model, the input data must be segmented such
        that the full-density
        bone has material ID i_max, with other IDs assigned linearly according
        to CT density. i_max
        is set with the option homminga_maximum_material_id,
        and exponent is set with the option homminga_modulus_exponent. E_max
        depends on the material defined (for an isotropic material it is set
        with the option youngs_modulus). It is possible to
        specify orthotropic material properties, in which case the
        3 shear modulii as well as the 3 Young's moduli will be scaled.""")

    matspec_group.add_argument ("--youngs_modulus", type=float,
        help="""Specify isotropic Young's modulus. 
    Units are MPa (for length units of mm).
    (default: %(default)s)""")

    matspec_group.add_argument ("--poissons_ratio", type=float,
        help="Specify isotropic Poisson's ratio. (default: %(default)s)")

    matspec_group.add_argument ("--orthotropic_parameters",
        help="""Specify orthotropic parameters as 
    Ex,Ey,Ez,nu_yz,nu_zx,nu_xy,G_yz,G_zx,G_xy.""")

    matspec_group.add_argument ("--homminga_maximum_material_id", type=int,
        help="""Specify the maximum material ID used in
    Homminga density to modulus conversion model.
    (default: %(default)s)""")

    matspec_group.add_argument ("--homminga_modulus_exponent", type=float,
       help="""Specify exponent used in Homminga 
    density to modulus conversion model.
    (default: %(default)s)""")


    # Material specification group
    additional_group = parser.add_argument_group('additional arguments')
    additional_group.add_argument ("--connectivity_filter",
        choices=["on", "off"],
        help="""Enable/disable connectivity filtering 
    of the generated mesh which extracts 
    only the largest connected object in 
    the input image.
    (default: %(default)s)""")
    additional_group.add_argument ("--spacing",
        help="""Force the spacing of the input image. This is sometimes required for
    DICOM files, for which the z-spacing sometimes cannot be determined.
    Values should be given as a comma-delimited list of x,y,z values.
    e.g. 0.5,0.5,0.5 .""")

    # Now parse the command line (except for config file name)
    args = parser.parse_args(remaining_args)

    # If no actions specified, select them all
    all_actions = ["generate", "solve", "analyze"]
    if not any(map(lambda a: args.__dict__[a], all_actions)):
        for a in all_actions:
            args.__dict__[a] = True


    # -------------------------------------------------------------------------
    # Some useful values and constants

    basename, extension = os.path.splitext(args.input_file)
    extension = extension.lower()
    strain = 1.00

    normal_strain_files = [ "%s_strain_xx.n88model" % basename,
                            "%s_strain_yy.n88model" % basename,
                            "%s_strain_zz.n88model" % basename ]
    shear_strain_files  = [ "%s_strain_yz.n88model" % basename,
                            "%s_strain_zx.n88model" % basename,
                            "%s_strain_xy.n88model" % basename ]

    indent = " "*4  

    # Create an instance
    errorObserver = ErrorObserver()


    # -------------------------------------------------------------------------
    # Generate the six required models: 3 compression tests and 3 shear tests.

    if args.generate:
        
        # -------------------------------------------------------------------------
        # Summarize configuration for reporting.
        # Report default values where used, but only report relevant values.
        
        log ("""n88directmechanics version %s

    """ % n88tools_version, timestamp=False)

        def print_setting (name, value):
            log ("%-28s = %s\n" % (name, value), timestamp=False)

        print_setting ("input_file", args.input_file)
        print_setting ("material_table", args.material_table)
        if args.orthotropic_parameters != None:
            print_setting ("orthotropic_parameters", args.orthotropic_parameters)
        else:
            print_setting ("youngs_modulus", args.youngs_modulus)
            print_setting ("poissons_ratio", args.poissons_ratio)
        if args.material_table == "homminga":
            print_setting ("homminga_maximum_material_id", args.homminga_maximum_material_id)
            print_setting ("homminga_modulus_exponent", args.homminga_modulus_exponent)
        print_setting ("connectivity_filter", args.connectivity_filter)

        log ("\n", timestamp=False)

        # -------------------------------------------------------------------------
        # Read input data

        if os.path.isdir(args.input_file):
            reader = vtk.vtkDICOMImageReader()        
        elif extension == ".aim":
            reader = vtkbone.vtkboneAIMReader()
            reader.DataOnCellsOn()
        elif extension == ".mha" or extension == ".mhd":
            reader = vtk.vtkMetaImageReader()
        elif extension == ".vti":
            reader = vtk.vtkXMLImageDataReader()
        else:
            raise N88ReportedError ("Unknown input file type.")

        log ("Reading image file " + args.input_file)
        reader.AddObserver ("ErrorEvent", errorObserver)
        reader.SetFileName(args.input_file)
        reader.Update()
        if errorObserver.ErrorOccurred():
            log ("ERROR reading file: %s" % errorObserver.ErrorMessage())
            sys.exit (-1)
        image = reader.GetOutput()
        log ("Read %d points from image file." % image.GetNumberOfPoints())

        if args.spacing != None:
            spacing = args.spacing.split()
            if len(spacing) != 3:
                raise N88ReportedError ("ERROR:spacing must have 3 values")
            image.SetSpacing(float(spacing[0]),float(spacing[1]),float(spacing[2]))

        # -------------------------------------------------------------------------
        # Apply connectivity filter

        if args.connectivity_filter == "on":
            log ("Applying connectivity filter.")
            in_data_vtk = image.GetCellData().GetScalars()
            if in_data_vtk == None:
                in_data_vtk = image.GetPointData().GetScalars()
            in_data = vtk_to_numpy (in_data_vtk)
            connectivity_filter = vtkbone.vtkboneImageConnectivityFilter()
            connectivity_filter.SetExtractionModeToLargestRegion()
            connectivity_filter.SetInputData (image)
            connectivity_filter.Update()
            image = connectivity_filter.GetOutput()
            out_data_vtk = image.GetCellData().GetScalars()
            if out_data_vtk == None:
                out_data_vtk = image.GetPointData().GetScalars()
            out_data = vtk_to_numpy (out_data_vtk)
            maskedVoxelCount = sum ((in_data != 0) * (out_data == 0))
            log ("Masked out %d unconnected voxels." % maskedVoxelCount)

        # Image data can be either on cells or points
        image_data_vtk = image.GetCellData().GetScalars()
        if image_data_vtk is None:
            image_data_vtk = image.GetPointData().GetScalars()

        # -------------------------------------------------------------------------
        # Convert the 3D image data from the input file to hexahedral cells

        log("Converting to hexahedral cells.")
        geometry_generator = vtkbone.vtkboneImageToMesh()
        geometry_generator.SetInputData (image)
        geometry_generator.Update()
        geometry = geometry_generator.GetOutput()
        log("Generated %d hexahedrons" % geometry_generator.GetOutput().GetNumberOfCells())

        # --------------------------------------------------------------------------
        # Create a material table.

        log ("Creating material table.")

        # Create a material
        if args.orthotropic_parameters != None:
            material = vtkbone.vtkboneLinearOrthotropicMaterial()
            orthotropic_parameters = map(float, args.orthotropic_parameters.split(","))
            material.SetYoungsModulusX (orthotropic_parameters[0])
            material.SetYoungsModulusY (orthotropic_parameters[1])
            material.SetYoungsModulusZ (orthotropic_parameters[2])
            material.SetPoissonsRatioYZ(orthotropic_parameters[3])
            material.SetPoissonsRatioZX(orthotropic_parameters[4])
            material.SetPoissonsRatioXY(orthotropic_parameters[5])
            material.SetShearModulusYZ (orthotropic_parameters[6])
            material.SetShearModulusZX (orthotropic_parameters[7])
            material.SetShearModulusXY (orthotropic_parameters[8])
            material.SetName ("OrthotropicMaterial")
        else:
            # Default is a linear isotropic material
            material = vtkbone.vtkboneLinearIsotropicMaterial()
            material.SetYoungsModulus(args.youngs_modulus)
            material.SetPoissonsRatio(args.poissons_ratio)
            material.SetName ("IsotropicMaterial")
            
        # Create a material table
        if args.material_table == "homminga":
            mt_generator = vtkbone.vtkboneGenerateHommingaMaterialTable()
            mt_generator.SetFullScaleMaterial (material)
            mt_generator.SetLastIndex (args.homminga_maximum_material_id)
            mt_generator.SetExponent (args.homminga_modulus_exponent)
            mt_generator.Update()
            material_table = mt_generator.GetOutput()
        else:
            # Default is a homogeneous material.
            mt_generator = vtkbone.vtkboneGenerateHomogeneousMaterialTable()
            mt_generator.SetMaterial (material)
            mt_generator.SetMaterialIdList (image_data_vtk)
            mt_generator.Update()
            material_table = mt_generator.GetOutput()


        # --------------------------------------------------------------------------
        log ("Creating three uniaxial models.")

        def create_uniaxial(test_axis):

            generator = vtkbone.vtkboneApplyCompressionTest()
            generator.SetInputData(0, geometry)
            generator.SetInputData(1, material_table)
            generator.SetAppliedStrain(strain)
            generator.SetTestAxis(test_axis)
            generator.TopSurfaceContactFrictionOff()
            generator.BottomSurfaceContactFrictionOff()
            generator.ConfineSidesOn()
            generator.Update()
            model = generator.GetOutput()

            log( "Writing n88 model file: %s" % normal_strain_files[test_axis])
            model.SetHistory ("")  # Clear History set by vtkboneApplyCompressionTest
            model.AppendHistory ("Created by n88directmechanics version %s" % n88tools_version)
            model.SetLog ("")
            model.AppendLog (log_string)

            writer = vtkbone.vtkboneN88ModelWriter()
            writer.SetInputData(model)
            writer.SetFileName(normal_strain_files[test_axis])
            writer.Update()

        create_uniaxial (0)
        create_uniaxial (1)
        create_uniaxial (2)

        # --------------------------------------------------------------------------
        log ("Creating three symshear models.")

        def create_shear(test_axis):

            generator = vtkbone.vtkboneApplySymmetricShearTest()
            generator.SetInputData(0, geometry)
            generator.SetInputData(1, material_table)
            generator.SetShearStrain(strain)
            generator.SetTestAxis(test_axis)
            generator.ConfineTopAndBottomVerticallyOn()
            generator.ConfineSidesVerticallyOff()
            generator.Update()
            model = generator.GetOutput()

            log ("Writing n88 model file: %s" % shear_strain_files[test_axis])
            model.SetHistory ("")  # Clear History set by vtkboneApplySymmetricShearTest
            model.AppendHistory ("Created by n88directmechanics version %s" % n88tools_version)
            model.SetLog ("")
            model.AppendLog (log_string)

            writer = vtkbone.vtkboneN88ModelWriter()
            writer.SetInputData(model)
            writer.SetFileName(shear_strain_files[test_axis])
            writer.Update()

        create_shear (0)
        create_shear (1)
        create_shear (2)


    # --------------------------------------------------------------------------
    # Solve the 6 models

    if args.solve:

        # --------------------------------------------------------------------------
        # Call n88solver.

        def call_solver(mfile):
            log ("Calling solver on %s" % mfile)
            solver_args = ["n88solver_slt", mfile]
            p = subprocess.call( solver_args, stdout=subprocess.PIPE )
            if p != 0:
                raise N88ReportedError ("Solver returned error code %d" % p,
                                        value=p)
            p = subprocess.call( ["n88derivedfields", mfile], stdout=subprocess.PIPE )
            if p != 0:
                raise N88ReportedError ("n88derivedfields returned error code %d" % p,
                                        value=p)

        for f in normal_strain_files + shear_strain_files:
            call_solver(f)


    # --------------------------------------------------------------------------
    # Utility functions that will be used for post-processing

    def data_frame (sense, test_axis):
        """ Give a sense in the Test Frame, returns the equivalent sense in the Data
        Frame, given a test_axis setting.

        Works on both scalars and arrays.
        """
        return mod(sense+test_axis-2, 3)

    def test_frame (sense, test_axis):
        """ Give a sense in the Data Frame, returns the equivalent sense in the Test
        Frame, given a test_axis setting.

        Works on both scalars and arrays.
        """
        return mod(sense+2-test_axis, 3)


    # --------------------------------------------------------------------------
    # Analyze the 6 solved models

    if args.analyze:
        
        # --------------------------------------------------------------------------
        # Average stresses (across all elements)

        as_matrix = zeros((6,6), float)

        def calc_average_stress_components (test_index, fname):

            log("Reading n88 model file: %s" % fname)
            reader = vtkbone.vtkboneN88ModelReader()
            reader.SetFileName(fname)
            reader.Update()
            model = reader.GetOutput()

            points = vtk_to_numpy (model.GetPoints().GetData())

            element_points = vtk_to_numpy (model.GetAllCellPoints())
            element_points.shape = (model.GetNumberOfCells(), 8)
            number_elements = model.GetNumberOfCells()

            element_stresses = vtk_to_numpy (model.GetCellData().GetArray("Stress"))

            for i in arange(6):
                 as_matrix[i,test_index] = sum(element_stresses[:,i]) / number_elements

        def calc_volume_fraction (fname):
            reader = vtkbone.vtkboneN88ModelReader()
            reader.SetFileName(fname)
            reader.Update()
            model = reader.GetOutput()

            points = vtk_to_numpy (model.GetPoints().GetData())

            element_points = vtk_to_numpy (model.GetAllCellPoints())
            element_points.shape = (model.GetNumberOfCells(), 8)

            # Calculate partial volume
            bounds = array(model.GetBounds())
            dimensions = bounds[1::2] - bounds[::2]
            volume = product(dimensions)
            element_dimensions = points[element_points[0,7]] - points[element_points[0,0]]
            number_elements = model.GetNumberOfCells()

            return number_elements * product(element_dimensions) / volume

        for i in arange(6):
             calc_average_stress_components (i, (normal_strain_files + shear_strain_files)[i])

        volume_fraction = calc_volume_fraction(normal_strain_files[0])
        log ("\n" + indent + "Volume fraction = %.5f" % volume_fraction, timestamp = False)

        # Apparent stiffness matrix, E, is determined by mult of average stress and vol fraction
        Eapp_matrix = as_matrix * volume_fraction

        # Ensure symmetry
        Eapp_matrix_sym = ( Eapp_matrix + Eapp_matrix.transpose() ) / 2

        numpy.set_printoptions(precision=3, suppress=True)
        log ("\n" + indent +
             "Apparent stiffness matrix in specimen coordinate system\n\n" +
             Eapp_matrix_sym.__str__(),
             timestamp=False)

        # We calculate the C_apparent by inverting E_apparent
        C_app = linalg.inv(Eapp_matrix_sym)
        numpy.set_printoptions(precision=3, suppress=False)
        log ("\n" + indent +
             "Apparent compliance matrix in specimen coordinate system\n\n" +
             C_app.__str__(),
             timestamp=False)

        # Engineering constants
        
        def print_engineering_constants(C):
            log (
                indent + "-------------------------------------------\n" +
                indent + "%s = %8.2f\n" % ("Exx",(1/C[0,0])) +
                indent + "%s = %8.2f\n" % ("Eyy",(1/C[1,1])) +
                indent + "%s = %8.2f\n" % ("Ezz",(1/C[2,2])) +
                indent + "%s = %8.2f\n" % ("Gyz",(1/C[3,3])) +
                indent + "%s = %8.2f\n" % ("Gzx",(1/C[4,4])) +
                indent + "%s = %8.2f\n" % ("Gxy",(1/C[5,5])) +
                indent + "-------------------------------------------\n" +
                indent + "%s = %8.5f\n" % ("nu_yx",(1/C[1,1] * -C[0,1])) +
                indent + "%s = %8.5f\n" % ("nu_zx",(1/C[2,2] * -C[0,2])) +
                indent + "%s = %8.5f\n" % ("nu_xy",(1/C[0,0] * -C[1,0])) +
                indent + "%s = %8.5f\n" % ("nu_zy",(1/C[2,2] * -C[1,2])) +
                indent + "%s = %8.5f\n" % ("nu_xz",(1/C[0,0] * -C[2,0])) +
                indent + "%s = %8.5f\n" % ("nu_yz",(1/C[1,1] * -C[2,1])) +
                indent + "-------------------------------------------\n" ,
                timestamp=False)

        log ("\n" + indent + "Material parameters in specimen coordinate system\n",
             timestamp=False)
        print_engineering_constants(C_app)
        
        # --------------------------------------------------------------------------
        # Numerically diagonalize w.r.t. a rotation

        def to_tensor(M):
            T = zeros((3,3,3,3),float)
            i,j = numpy.mgrid[0:3,0:3]
            T[ i,      i,      j,      j     ] = M[i  ,j  ]
            T[ i,      i,     (j+1)%3,(j+2)%3] = M[i  ,j+3]
            T[ i,      i,     (j+2)%3,(j+1)%3] = M[i  ,j+3]
            T[(i+1)%3,(i+2)%3, j,      j     ] = M[i+3,j  ]
            T[(i+2)%3,(i+1)%3, j,      j     ] = M[i+3,j  ]
            T[(i+1)%3,(i+2)%3,(j+1)%3,(j+2)%3] = M[i+3,j+3]
            T[(i+1)%3,(i+2)%3,(j+2)%3,(j+1)%3] = M[i+3,j+3]
            T[(i+2)%3,(i+1)%3,(j+1)%3,(j+2)%3] = M[i+3,j+3]
            T[(i+2)%3,(i+1)%3,(j+2)%3,(j+1)%3] = M[i+3,j+3]
            return T

        def from_tensor(T):
            M = zeros((6,6),float)
            i,j = numpy.mgrid[0:3,0:3]
            M[i  ,j  ] = T[ i,      i,      j,      j     ]
            M[i  ,j+3] = T[ i,      i,     (j+1)%3,(j+2)%3]
            M[i+3,j  ] = T[(i+1)%3,(i+2)%3, j,      j     ]
            M[i+3,j+3] = T[(i+1)%3,(i+2)%3,(j+1)%3,(j+2)%3]
            return M

        def rotate_M (M, R):
            T = to_tensor(M)
            T = numpy.einsum('pi,qj,rk,sl,ijkl', R, R, R, R, T)
            return from_tensor(T)

        def rotation_matrix (angles):
            Rx = transformations.rotationX(angles[0])
            Ry = transformations.rotationY(angles[1])
            Rz = transformations.rotationZ(angles[2])
            return dot(dot(Rz,Ry),Rx)

        def obj_func_R (R, M):
            i = arange(6)
            Mrot = rotate_M (M, R)
            e_sum = sum(Mrot[:3,:3].flatten()**2) \
                  + sum(Mrot[arange(3,6),arange(3,6)]**2)
            others_sum = sum(Mrot.ravel()**2) - e_sum
            return (others_sum / e_sum)

        def obj_func (angles, M):
            return obj_func_R (rotation_matrix(angles), M)

        # Optimization
        angles = array([0, 0, 0])  # Initial value
        angles = minimize(obj_func, angles, args=(Eapp_matrix_sym,), xtol=1e-10, ftol=1E-6, disp=False)
        R = rotation_matrix(angles)
        Eapp_matrix_sym_rot = rotate_M (Eapp_matrix_sym, R)
        
        # Now we're going to permute the axes so that Exx >= Eyy >= Ezz
        key = argsort(abs(Eapp_matrix_sym_rot[arange(3),arange(3)]))[::-1]
        P = zeros((3,3))
        P[key,arange(3)] = 1
        # Ensure we have an even orthogonal permutation (i.e. a rotation):
        # an arbitrary choice is to multiply the last row by the determinant
        P[2,:] *= linalg.det(P)
        R = dot(P,R)
        Eapp_matrix_sym_rot = rotate_M (Eapp_matrix_sym, R)

        numpy.set_printoptions(precision=5, suppress=True)
        log ("\n" + indent +
             "Optimum rotation matrix R:\n\n" +
             R.__str__(),
             timestamp=False)
        Eapp_matrix_sym_rot = rotate_M (Eapp_matrix_sym, R)
        numpy.set_printoptions(precision=3)
        log ("\n" + indent +
             "Apparent stiffness matrix in best orthotropic coordinate system\n\n" +
             Eapp_matrix_sym_rot.__str__(),
             timestamp=False)

        C_app_rot = linalg.inv(Eapp_matrix_sym_rot)
        numpy.set_printoptions(precision=3, suppress=False)
        log ("\n" + indent +
             "Apparent compliance matrix in best orthotropic coordinate system\n\n" +
             C_app_rot.__str__(),
             timestamp=False)

        log ("\n" + indent +
             "Material parameters in best orthotropic coordinate system\n\n",
             timestamp=False)
        print_engineering_constants(C_app_rot)
        log ("\n", timestamp=False)


def main():
    try:
        directmechanics()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
