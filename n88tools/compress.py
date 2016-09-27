"""
compress.py

A tool to compress n88model files.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from N88ReportedError import N88ReportedError

def compress():

    import os
    import argparse
    import vtk
    import vtkbone

    parser = argparse.ArgumentParser (
        prog="n88compress",
        description="Compress a finite element model file.")

    parser.add_argument ("input")

    args = parser.parse_args()

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

    errorObserver = ErrorObserver()

    reader = vtkbone.vtkboneN88ModelReader()
    writer = vtkbone.vtkboneN88ModelWriter()

    print "Reading n88model file : %s" % args.input
    reader.AddObserver ("ErrorEvent", errorObserver)
    reader.SetFileName (args.input)
    reader.Update()
    if errorObserver.ErrorOccurred():
        raise N88ReportedError ("ERROR reading file: " +
             errorObserver.ErrorMessage())
    model = reader.GetOutput()

    output = args.input
    print "Re-writing n88model file : %s" % output
    writer.SetInputData (model)
    writer.SetFileName (output)
    writer.CompressionOn()
    writer.Update()


def main():
    try:
        compress()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
