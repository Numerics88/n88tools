"""
copymodel.py

A tool to copy and convert finite element model formats. Also supports
compression.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from .N88ReportedError import N88ReportedError

def copymodel():

    import os
    import argparse
    import vtk
    import vtkbone

    supported_input_list = \
"""    .n88model - Numerics88 model file
    .inp - Abaqus input file
    .inp - Faim version 5 input file
    .dat - Faim version 5 output file"""

    supported_output_list = \
"""    .n88model - Numerics88 model file
    .inp - Abaqus input file
    .vtu - VTK XML Unstructured Grid file"""

    parser = argparse.ArgumentParser (
        prog="n88copymodel",
        description="""Copy, convert or compress a finite element model file.""",
        epilog=("supported input formats:\n" + supported_input_list +
                "\n\nsupported output formats:\n" + supported_output_list +
"""

File formats are identified by their extensions. Ambiguous input file
extensions are resolved by examining the file.

Two input files are allowed if the first is a FAIM version 5 input file,
and the second is a FAIM version 5 output file. This will create
a complete solved N88 model file. Either can also be converted individually,
but in the case of .dat files, the resulting N88 model file will be incomplete,
as the .dat file alone does not contain either material definitions or
constraints; it can be rendered and processed by n88postfaim, but not re-solved.
"""),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument (
        "--compress", "-c",
        action="store_true",
        help="""Use compression when writing the output file if the file format supports it.
Compressed n88model files do not need
to be uncompressed in order to use them: they can be
used in all cases exactly like uncompressed files. However,
they may be slower to read, and are particularly slower to write."""
        )
    parser.add_argument ("input")
    parser.add_argument ("input2", nargs='?')
    parser.add_argument ("output")

    args = parser.parse_args()

    reader2 = None
    extension = os.path.splitext(args.input)[1].lower()
    if extension == ".n88model":
        reader_format = "N88 Model"
        reader = vtkbone.vtkboneN88ModelReader()
    elif extension == ".inp":
        # Have to resolve the two possibilities: read first line of file
        f = open(args.input)
        line = f.readline().strip()
        f.close()
        faim_inp_identifier = "# Faim Version"
        if line[:len(faim_inp_identifier)] == faim_inp_identifier:
            reader_format = "Faim version 5 input format"
            reader = vtkbone.vtkboneFaimVersion5InputReader()
            if args.input2:
                if os.path.splitext(args.input2)[1].lower():
                    reader2_format = "Faim version 5 output format"
                    reader2 = vtkbone.vtkboneFaimVersion5OutputReader()
                else:
                    raise N88ReportedError (
                        "The second input with a .inp file must be a .dat file")
        else:
            reader_format = "Abaqus Input"
            reader = vtkbone.vtkboneAbaqusInputReader()
    elif extension == ".dat":
        reader_format = "Faim version 5 output format"
        reader = vtkbone.vtkboneFaimVersion5OutputReader()
        print("Warning: An incomplete model will be written, as Faim version 5 .dat")
        print("         files do not contain material data or constraints.")
        print("         The resulting model can neither be re-solved nor processed")
        print("         with n88postfaim.")
        print(" Tip:    If the corresponding .inp file is available, specify both the")
        print("         the .inp and .dat files as input.  This will create a complete")
        print("         solved model.")
    else:
        raise N88ReportedError (
"""Unsupported input file format.
Currently supported formats are:
""" + supported_input_list)

    extension = os.path.splitext(args.output)[1].lower()
    if extension == ".n88model":
        writer_format = "N88 Model"
        writer = vtkbone.vtkboneN88ModelWriter()
        if args.compress:
            writer.CompressionOn()
        writer_supports_mat_props = True
    elif extension == ".inp":
        writer_format = "Abaqus Input"
        writer = vtkbone.vtkboneAbaqusInputWriter()
        writer_supports_mat_props = True
    elif extension == ".vtu":
        writer_format = "VTK XML Unstructured Grid"
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer_supports_mat_props = False
    else:
        raise N88ReportedError (
"""Unsupported output file format.
Currently supported formats are:
""" + supported_output_list)


    # Don't waste time reading material properties if
    # the writer discards them.
    if (isinstance(reader, vtkbone.vtkboneN88ModelReader) and
            not writer_supports_mat_props):
        reader.ReadMaterialsOff()


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

    print("Reading %s file : %s" % (reader_format, args.input))
    reader.AddObserver ("ErrorEvent", errorObserver)
    reader.SetFileName (args.input)
    reader.Update()
    if errorObserver.ErrorOccurred():
        raise N88ReportedError ("ERROR reading file: " +
             errorObserver.ErrorMessage())
    model = reader.GetOutput()

    if reader2:
        print("Reading file : %s" % args.input2)
        reader2.SetFileName (args.input2)
        reader2.Update()
        solved_model = reader2.GetOutput()
        # Copy all the Element and Point data arrays from the solved model to the
        # input model.
        for i in range(solved_model.GetPointData().GetNumberOfArrays()):
            data = solved_model.GetPointData().GetArray(i)
            if data.GetName():
                model.GetPointData().AddArray(data)
        for i in range(solved_model.GetCellData().GetNumberOfArrays()):
            data = solved_model.GetCellData().GetArray(i)
            if data.GetName():
                model.GetCellData().AddArray(data)
        model.AppendLog("Solution values obtained from %s" % args.input2)

    print("Writing %s file : %s" % (writer_format, args.output))
    writer.SetInputData (model)
    writer.SetFileName (args.output)
    writer.Update()


def main():
    try:
        copymodel()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
