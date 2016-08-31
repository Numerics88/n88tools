"""
extractsets.py

A tool to extract node sets, element sets, and constraints from n88model files.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from N88ReportedError import N88ReportedError


def extractsets():

    import os
    import argparse
    import vtk
    import vtkbone
    from vtk.util.numpy_support import vtk_to_numpy


    parser = argparse.ArgumentParser (
        prog="n88extractsets",
        description="""Extract node sets, element sets and constraint sets
         and write them to VTK
    PolyData (.vtp) files which can be opened with ParaView. This is mostly useful for
    visualizing sets, boundary conditions and applied forces. For nodes sets and
    constraints, these files will consist of a collection of vertices.""",
        epilog=""" If no arguments are specified, all
    types of sets will be extracted.

    Warning: This utility will generate file names based on the
    names of the sets, and will overwrite existing files with the same names
    without warning.""",
        )

    parser.add_argument ("--constraints", "-C", action='store_true',
        help="Extract constraints.")
    parser.add_argument ("--node_sets", "-N", action='store_true',
        help="Extract node sets.")
    parser.add_argument ("--element_sets", "-E", action='store_true',
        help="Extract element sets.")

    parser.add_argument ("input")

    args = parser.parse_args()

    # If no arguments specified, select them all
    all_args = ["constraints","node_sets","element_sets"]
    if not any(map(lambda a: args.__dict__[a], all_args)):
        for a in all_args:
            args.__dict__[a] = True

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

    print "Reading %s file : %s" % ("N88 Model", args.input)
    reader = vtkbone.vtkboneN88ModelReader()
    reader.AddObserver ("ErrorEvent", errorObserver)
    reader.SetFileName (args.input)
    reader.Update()
    if errorObserver.ErrorOccurred():
        raise N88ReportedError ("ERROR reading file: " + errorObserver.ErrorMessage())
    model = reader.GetOutput()

    fileroot = os.path.splitext(os.path.basename(args.input))[0]

    # We will need this to be able to track back original Ids in the written-out data.
    vtkbone.vtkboneSelectionUtilities.AddPointPedigreeIdsArray(model)

    if args.constraints:
        constraints = model.GetConstraints()
        constraints.InitTraversal()
        constraint = constraints.GetNextConstraint()
        while constraint:
            name = constraint.GetName()
            output_file = "%s_constraint_%s.vtp" % (fileroot, name)
            print "Writing constraint :", output_file
            data = model.DataSetFromConstraint (constraint.GetName())
            reduceToPolyData = vtk.vtkGeometryFilter()
            reduceToPolyData.SetInputData (data)
            reduceToPolyData.MergingOff()
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName (output_file)
            writer.SetInputConnection (reduceToPolyData.GetOutputPort())
            writer.Update()
            constraint = constraints.GetNextConstraint()

    if args.node_sets:
        node_sets = model.GetNodeSets()
        node_sets.InitTraversal()
        node_set = node_sets.GetNextItem()
        while node_set:
            name = node_set.GetName()
            output_file = "%s_node_set_%s.vtp" % (fileroot, name)
            print "Writing node set :", output_file
            data = model.DataSetFromNodeSet (node_set.GetName())
            reduceToPolyData = vtk.vtkGeometryFilter()
            reduceToPolyData.SetInputData (data)
            reduceToPolyData.MergingOff()
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName (output_file)
            writer.SetInputConnection (reduceToPolyData.GetOutputPort())
            writer.Update()
            node_set = node_sets.GetNextItem()

    if args.element_sets:
        element_sets = model.GetElementSets()
        element_sets.InitTraversal()
        element_set = element_sets.GetNextItem()
        while element_set:
            name = element_set.GetName()
            output_file = "%s_element_set_%s.vtp" % (fileroot, name)
            print "Writing element set :", output_file
            data = model.DataSetFromElementSet (element_set.GetName())
            reduceToPolyData = vtk.vtkGeometryFilter()
            reduceToPolyData.SetInputData (data)
            reduceToPolyData.MergingOff()
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName (output_file)
            writer.SetInputConnection (reduceToPolyData.GetOutputPort())
            writer.Update()
            element_set = element_sets.GetNextItem()


def main():
    try:
        extractsets()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
