"""
coarsen.py

A tool to reduce the size of finite element models so that an
approximate solution can be obtained quickly and with less memory.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from N88ReportedError import N88ReportedError

try:
    import pkg_resources
    n88tools_version = pkg_resources.require("n88tools")[0].version
except:
    n88tools_version = "unknown"


def coarsen():

    import argparse
    import os
    import numpy
    import vtk
    import vtkbone

    # ------------------------------------------------------------------------
    # Parse options

    parser = argparse.ArgumentParser (
        prog="n88coarsen",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=
"""This tool is used to reduce the size of finite element models so that an
approximate solution can be obtained quickly and with less memory. It
operates either directly on FE models, or on segmented image data.

When operating on n88model files, it will increase the size of elements by a
factor of 2 in each dimension (thus a factor of 8 in volume), resulting in
fewer larger elements. The output model is coarser: the coarsening is always
done by adding volume to make larger elements, never by removing material, so
that as compared with the input model, the output model always has greater or
equal volume. Material properties are averaged over all the input elements
corresponding to each output element. Empty space is treated as having
identically zero stiffness, wherever empty space in the input is encompassed
within an output element. All essential features of the model, including
boundary conditions, applied forces, and post-processing sets, are translated
to the new coarser mesh.

When the input is an image, the resolution is reduced by exactly 1/2 in each
linear dimension. Each 2x2x2 cube in the input becomes a single voxel in the
output, with value equal to the maximum in the corresponding input 2x2x2.
This operation does not average input values, since it does not make any
sense to average segmentation values, which are just labels. Instead, it
takes the maximum value over all the input voxels corresponding to an output
voxel. This is a somewhat arbitrary choice. For this reason it is preferable
to use n88coarsen directly on n88model files, where material averaging is
possible.

Supported image input formats:
    DICOM (a directory)
    Scanco AIM (.aim)
    MetaImage (.mha or .mhd)
    VTK XML ImageData (.vti)

Supported image output formats:
    VTK XML ImageData (.vti)

""")

    parser.add_argument ("--material_averaging", choices=["linear", "homminga_density"],
        default="homminga_density",
        help="""Determine how the material averaging is done. If 'linear', then stress-strain 
matrices will be linearly averaged. If 'homminga_density', then the stress-strain matrices are
first scaled to a density using the Homminga formula (i.e. raised to the power
1/1.7), then averaged, then converted back the stiffness by raising to the power 1.7.
'homminga_density' nearly always gives more accurate approximations than 'linear'.""")

    parser.add_argument ("input_file", default=None)
    parser.add_argument ("output_file", default=None)

    args = parser.parse_args()


    # -------------------------------------------------------------------------
    # Some useful values and constants

    in_basename, in_extension = os.path.splitext(args.input_file)
    in_extension = in_extension.lower()
    in_basename = os.path.split(in_basename)[1]

    out_basename, out_extension = os.path.splitext(args.output_file)
    out_extension = out_extension.lower()

    # -------------------------------------------------------------------------
    # Error handling

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

    # Create an instance
    errorObserver = ErrorObserver()

    # -------------------------------------------------------------------------
    # Determine input and output types

    input_is_model = False

    if in_extension == ".n88model":
        reader = vtkbone.vtkboneN88ModelReader()
        input_is_model = True
        writer = vtkbone.vtkboneN88ModelWriter()
    else:
        if os.path.isdir(args.input_file):
            reader = vtk.vtkDICOMImageReader()        
        elif in_extension == ".aim":
            reader = vtkbone.vtkboneAIMReader()
            reader.DataOnCellsOn()
        elif in_extension == ".mha" or in_extension == ".mhd":
            reader = vtk.vtkMetaImageReader()
        elif in_extension == ".vti":
            reader = vtk.vtkXMLImageDataReader()
        else:
            raise N88ReportedError ("Unknown input file type.")
        writer = vtk.vtkXMLImageDataWriter()

    if input_is_model:

        if out_extension != ".n88model":
            raise N88ReportedError ("ERROR: Output must be .n88model file if input is n88model.")

        # -------------------------------------------------------------------------
        # Read input data

        print ("Reading model file " + args.input_file)
        reader.AddObserver ("ErrorEvent", errorObserver)
        reader.SetFileName(args.input_file)
        reader.Update()
        if errorObserver.ErrorOccurred():
            raise N88ReportedError ("ERROR reading file: %s" % errorObserver.ErrorMessage())
        input_model = reader.GetOutput()

        nElements = numpy.array(input_model.GetNumberOfCells())
        print ("Input number of elements: %s" % (nElements,))

        # -------------------------------------------------------------------------
        # Apply filter

        print ("Generating coarsened model")
        filter = vtkbone.vtkboneCoarsenModel()
        filter.AddObserver ("ErrorEvent", errorObserver)
        filter.SetInputData (input_model)
        if args.material_averaging == "homminga_density":
            filter.SetMaterialAveragingMethod (vtkbone.vtkboneCoarsenModel.HOMMINGA_DENSITY)
        elif args.material_averaging == "linear":
            filter.SetMaterialAveragingMethod (vtkbone.vtkboneCoarsenModel.LINEAR)
        else:
            N88ReportedError ("Unexpected value for material_averaging.")
        filter.Update()
        if errorObserver.ErrorOccurred():
            raise N88ReportedError ("ERROR: %s" % errorObserver.ErrorMessage())
        output_model = filter.GetOutput()

        nElements = numpy.array(output_model.GetNumberOfCells())
        print ("Output number of elements: %s" % (nElements,))

        # -------------------------------------------------------------------------
        # Set history and log

        output_model.SetHistory ("")  # Clear History set by vtkboneCoarsenModel
        output_model.AppendHistory ("Created by n88coarsen version %s" % n88tools_version)
        log_string = "n88coarsen settings:\n"
        log_string += "  input_file = " + args.input_file + "\n"
        log_string += "  material_averaging = " + args.material_averaging + "\n"
        output_model.SetLog ("")
        output_model.AppendLog (log_string)

        # -------------------------------------------------------------------------
        # Write data

        print ("Writing model file " + args.output_file)
        writer.AddObserver ("ErrorEvent", errorObserver)
        writer.SetInputData (output_model)
        writer.SetFileName(args.output_file)
        writer.Update()
        if errorObserver.ErrorOccurred():
            raise N88ReportedError ("ERROR writing file: %s" % errorObserver.ErrorMessage())

    else:

        if out_extension != ".vti":
            raise N88ReportedError ("ERROR: Output must be .vti file if input is image.")

        # -------------------------------------------------------------------------
        # Read input data

        print ("Reading image file " + args.input_file)
        reader.AddObserver ("ErrorEvent", errorObserver)
        reader.SetFileName(args.input_file)
        reader.Update()
        if errorObserver.ErrorOccurred():
            raise N88ReportedError ("ERROR reading file: %s" % errorObserver.ErrorMessage())
        input_image = reader.GetOutput()

        dims = numpy.array(input_image.GetDimensions())
        if input_image.GetCellData().GetScalars() is None:
            print "Data on points"
        else:
            print "Data on cells"
            dims -= 1
        print ("Input data dimensions: %s" % (dims,))

        # -------------------------------------------------------------------------
        # Apply filter

        print ("Generating reduced resolution image")
        filter = vtkbone.vtkboneDecimateImage()
        filter.AddObserver ("ErrorEvent", errorObserver)
        filter.SetInputData (input_image)
        filter.Update()
        if errorObserver.ErrorOccurred():
            raise N88ReportedError ("ERROR: %s" % errorObserver.ErrorMessage())
        output_image = filter.GetOutput()

        dims = numpy.array(output_image.GetDimensions())
        if not (output_image.GetCellData().GetScalars() is None):
            dims -= 1
        print ("Output data dimensions: %s" % (dims,))

        # -------------------------------------------------------------------------
        # Write data

        print ("Writing image file " + args.output_file)
        writer.AddObserver ("ErrorEvent", errorObserver)
        writer.SetInputData (output_image)
        writer.SetFileName(args.output_file)
        writer.Update()
        if errorObserver.ErrorOccurred():
            raise N88ReportedError ("ERROR writing file: %s" % errorObserver.ErrorMessage())


def main():
    try:
        coarsen()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
