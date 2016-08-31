"""
extractfields.py

A tool to extract field data from solutions in n88model files and write
the data in text format.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from N88ReportedError import N88ReportedError
import numpy
from numpy.core import *


def extractfields():

    import argparse
    from netCDF4 import Dataset

    parser = argparse.ArgumentParser (
        prog="n88extractfields",
        description="""Extracts specified solution fields as tabular text data.""",
        epilog="""examples:

  This will extract the displacement values as an Nx3 array, where N is
  the number of nodes:

      n88extractfields Displacement mymodel.n88model

  This example will extract both element number and corresponding stress.
  As the stress is 6-valued, this makes a total of 7 values per row of
  output:
        
      n88extractfields ElementNumber,Stress mymodel.n88model

tip:
  To determine which solution fields are available to be extracted, run
  "n88modelinfo --solutions" on your model file.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    parser.add_argument ("--output_file", "-o",
        help="Output file. If not specified, output will go to STDOUT.")


    parser.add_argument ("fields",
        help="""A comma-delimited list of names of solution fields to extract.
          May be the name of
    either a solution node
    value (e.g. "Displacement") or a solution element value (e.g. "Stress").
    Additionally, "NodeNumber", "ElementNumber", "MaterialID", "NodeCoordinates", 
    "ElementCoordinates" or "Topology" can be specified.
    Element coordinates are the centers of the elements.
    Topology gives for each element the 8 node numbers of the nodes
    constituting the element. The order is as in the n88model file; refer
    to the manual.
    All fields must be the same length, which practically means that they
    must be either all node values or all element values.    
    """)
    parser.add_argument ("input_file", help="An n88 model file.")

    args = parser.parse_args()

    if args.output_file == None:
        out = sys.stdout
    else:
        out = open (args.output_file, "wt")

    rootGroup = Dataset (args.input_file, 'r')

    try:
        activeSolutionGroup = rootGroup.groups['Solutions'].groups[rootGroup.ActiveSolution]
    except:
        activeSolutionGroup = None
    if activeSolutionGroup != None:
        activeProblemGroup = rootGroup.groups['Problems'].groups[activeSolutionGroup.Problem]
    else:
        activeProblemGroup = rootGroup.groups['Problems'].groups[rootGroup.ActiveProblem]
    activePartGroup = rootGroup.groups['Parts'].groups[activeProblemGroup.Part]
        
    try:
        activeProblemGroup = rootGroup.groups['Problems'].groups[rootGroup.ActiveProblem]
        activePartGroup = rootGroup.groups['Parts'].groups[activeProblemGroup.Part]
    except (AttributeError, KeyError):
        raise N88ReportedError ("No active part in file.")


    def get_field (fieldName):

        if fieldName == "NodeNumber":
            try:
                fieldVariable = activePartGroup.variables["NodeCoordinates"]
            except KeyError:
                raise N88ReportedError ("NodeCoordinates of part not found.")
            fieldVariable = arange(len(fieldVariable)) + 1

        elif fieldName == "ElementNumber":
            try:
                fieldVariable = activePartGroup.groups["Elements"].groups["Hexahedrons"].variables["ElementNumber"]
            except (AttributeError, KeyError):
                raise N88ReportedError ("ElementNumber not found.")
        
        elif fieldName == "Topology":
            try:
                fieldVariable = activePartGroup.groups["Elements"].groups["Hexahedrons"].variables["NodeNumbers"]
            except (AttributeError, KeyError):
                raise N88ReportedError ("NodeNumbers not found.")

        elif fieldName == "MaterialID":
            try:
                fieldVariable = activePartGroup.groups["Elements"].groups["Hexahedrons"].variables["MaterialID"]
            except (AttributeError, KeyError):
                raise N88ReportedError ("MaterialID not found.")

        elif (fieldName == "NodeCoordinates") or (fieldName == "ElementCoordinates"):
            try:
                fieldVariable = activePartGroup.variables["NodeCoordinates"]
            except KeyError:
                raise N88ReportedError ("NodeCoordinates of part not found.")

            if fieldName == "ElementCoordinates":
                # Actually have to calculate these
                nodeCoord = fieldVariable[:]
                try:
                    nodeNumbers = activePartGroup.groups["Elements"].groups["Hexahedrons"].variables["NodeNumbers"]
                except (AttributeError, KeyError):
                    raise N88ReportedError ("NodeNumbers not found.")
                fieldVariable = 0.5*(nodeCoord[nodeNumbers[:,0]-1,:] + nodeCoord[nodeNumbers[:,7]-1,:])

        else:
            if activeSolutionGroup == None:
                raise N88ReportedError ("No active solution in file.")
            try:
                fieldVariable = activeSolutionGroup.groups["NodeValues"].variables[fieldName]
            except KeyError:
                try:
                    fieldVariable = activeSolutionGroup.groups["ElementValues"].variables[fieldName]
                except KeyError:
                    raise N88ReportedError ("Solution field \"%s\" not found." % fieldName)
        
        return fieldVariable


    variables = []
    for name in args.fields.split(","):
        variables.append (get_field(name))

    for v in variables[1:]:
        if len(v) != len(variables[0]):
            raise N88ReportedError ("Fields of different lengths.  Did you mix node and elements values?")

    for i in xrange(len(variables[0])):
        values = []
        for v in variables:
            values += map(str, numpy.atleast_1d(v[i]))
        out.write ("\t".join(values) + "\n")


def main():
    try:
        extractfields()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
