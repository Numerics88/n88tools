"""
interpolatesolution.py

A tool to reduce the size of finite element models so that an
approximate solution can be obtained quickly and with less memory.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""


from __future__ import division
import sys
from .N88ReportedError import N88ReportedError


def interpolatesolution():

    import argparse
    import numpy
    import vtk
    import vtkbone
    import os


    # ------------------------------------------------------------------------
    # Parse options

    parser = argparse.ArgumentParser (
        prog="n88interpolatesolution",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=
    """Copies and interpolates a solution from a solved reduced-resolution n88model
file, produced using n88coarsen, into the original n88model file. This may
help to obtain faster solutions.

WARNING: For linear models only.
""")

    parser.add_argument ("full_model",
        help="An unsolved n88model file.",
        default=None)
    parser.add_argument ("reduced_model",
        help="A solved n88model file, which was generated from the full_model file with n88coarsen.",
        default=None)

    args = parser.parse_args()

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
    # Read input data

    reader1 = vtkbone.vtkboneN88ModelReader()        
    print("Reading full resolution model file " + args.full_model)
    reader1.AddObserver ("ErrorEvent", errorObserver)
    reader1.SetFileName(args.full_model)
    reader1.Update()
    if errorObserver.ErrorOccurred():
        raise N88ReportedError ("ERROR reading file: %s" % errorObserver.ErrorMessage())
    full_model = reader1.GetOutput()

    reader2 = vtkbone.vtkboneN88ModelReader()        
    print("Reading reduced resolution model file " + args.reduced_model)
    reader2.AddObserver ("ErrorEvent", errorObserver)
    reader2.SetFileName(args.reduced_model)
    reader2.Update()
    if errorObserver.ErrorOccurred():
        raise N88ReportedError ("ERROR reading file: %s" % errorObserver.ErrorMessage())
    reduced_model = reader2.GetOutput()


    # -------------------------------------------------------------------------
    # Apply filter

    print("Interpolated reduced resolution solution.")
    filter = vtkbone.vtkboneInterpolateCoarseSolution()
    filter.AddObserver ("ErrorEvent", errorObserver)
    filter.SetInputData (0, full_model)
    filter.SetInputData (1, reduced_model)
    filter.Update()
    if errorObserver.ErrorOccurred():
        raise N88ReportedError ("ERROR: %s" % errorObserver.ErrorMessage())
    data = filter.GetSolutionArray()

    full_model.GetPointData().RemoveArray("Displacement")
    full_model.GetPointData().AddArray(data)


    # -------------------------------------------------------------------------
    # Write data

    writer = vtkbone.vtkboneN88ModelWriter()
    print("Re-writing full model.")
    writer.AddObserver ("ErrorEvent", errorObserver)
    writer.SetInputData (full_model)
    writer.SetFileName(args.full_model)
    writer.Update()
    if errorObserver.ErrorOccurred():
        raise N88ReportedError ("ERROR writing file: %s" % errorObserver.ErrorMessage())


def main():
    try:
        interpolatesolution()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
