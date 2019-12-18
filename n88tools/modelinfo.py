"""
modelinfo.py

A tool to print summary data about an n88model file.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from .N88ReportedError import N88ReportedError

def modelinfo():

    import os
    import argparse
    from netCDF4 import Dataset

    parser = argparse.ArgumentParser (
        prog="n88modelinfo",
        description="Print summary information for a Numerics88 finite element model file.",
        epilog="Multiple action arguments may be specified.  If no action argument"
               " is specified, it is equivalent to specifying them all." 
        )

    parser.add_argument ("input_file", help="The .n88model file to read.")

    # Actions
    action_group = parser.add_argument_group('action arguments')
    action_group.add_argument ("--active", action='store_true', help="list active solution, problem and part")
    action_group.add_argument ("--history", action='store_true', help="show history")
    action_group.add_argument ("--log", action='store_true', help="show log")
    action_group.add_argument ("--materials", action='store_true', help="list defined materials")
    action_group.add_argument ("--parts", action='store_true', help="show parts")
    action_group.add_argument ("--node_sets", action='store_true', help="show node sets")
    action_group.add_argument ("--element_sets", action='store_true', help="show element sets")
    action_group.add_argument ("--sets", action='store_true',
        help="show both node and element sets (equivalent to node_sets and element_sets)")
    action_group.add_argument ("--constraints", action='store_true', help="show constraints")
    action_group.add_argument ("--problems", action='store_true', help="show problems")
    action_group.add_argument ("--solutions", action='store_true', help="show solutions")

    args = parser.parse_args()

    if args.sets:
        args.node_sets = True
        args.element_sets = True

    # If no arguments specified, select them all
    all_args = ["active",
                "history",
                "log",
                "materials",
                "parts",
                "node_sets",
                "element_sets",
                "constraints",
                "problems",
                "solutions",
                ]
    if not any(map(lambda a: args.__dict__[a], all_args)):
        for a in all_args:
            args.__dict__[a] = True


    extension = os.path.splitext(args.input_file)[1].lower()
    if extension != ".n88model":
        raise N88ReportedError ("Unsupported input file format.  Only .n88model files are supported.")

    divider = "-"*70

    rootGroup = Dataset (args.input_file, 'r')

    def enumerate_group (group, list_vars=False, postfunc=None):
        for name,subgroup in group.groups.items():
            print('')
            print("  Name : " + name)
            for attrName,attrValue in subgroup.__dict__.items():
                print("  " + attrName + " : " + str(attrValue))
            for dimName,dimValue in subgroup.dimensions.items():
                print("  " + dimName + " : " + str(len(dimValue)))
            if list_vars:
                variables = []
                for varName in subgroup.variables.keys():
                    variables.append (varName)
                print("  Variables : " + ' '.join(variables))
            if postfunc:
                postfunc(subgroup)

    def enumerate_subgroups (parent_group, group_name, list_vars=False, postfunc=None):
        if group_name in parent_group.groups:
            enumerate_group (parent_group.groups[group_name], list_vars=list_vars, postfunc=postfunc)

    first = True

    if args.history:
        if first:
            first = False
        else:
            print('')
        print("History: ")
        print(divider)
        if 'History' in rootGroup.__dict__:
            print(rootGroup.History)
        print(divider)

    if args.log:
        if first:
            first = False
        else:
            print('')
        print("Log: ")
        print(divider)
        # Support older Comments field as well
        if 'Comments' in rootGroup.__dict__:
            print(rootGroup.Comments)
        if 'Log' in rootGroup.__dict__:
            print(rootGroup.Log)
        print(divider)

    if args.active:
        if first:
            first = False
        else:
            print('')
        print("Active Settings:")
        print(divider)
        if 'ActiveSolution' in rootGroup.__dict__:
            activeSolution = rootGroup.ActiveSolution
            # Get active pproblem from active solution
            activeSolutionGroup = rootGroup.groups['Solutions'].groups[activeSolution]
            activeProblem = activeSolutionGroup.Problem
        else:
            activeSolution = None
            activeProblem = rootGroup.ActiveProblem
        activeProblemGroup = rootGroup.groups['Problems'].groups[activeProblem]
        activePart = activeProblemGroup.Part
        print("  Active Solution : " + str(activeSolution))
        print("  Active Problem : " + activeProblem)
        print("  Active Part : " + activePart)
        print(divider)

    if args.materials:
        if first:
            first = False
        else:
            print('')
        print("Materials:")
        print(divider)
        enumerate_subgroups (rootGroup, 'MaterialDefinitions')
        print(divider)

    def parts_postfunc (group):
        if "Elements" in group.groups:
            subgroup = group.groups["Elements"]
            for name,subsubgroup in subgroup.groups.items():
                print("  " + name + " :")
                for dimName,dimValue in subsubgroup.dimensions.items():
                    print("    " + dimName + " : " + str(len(dimValue)))

    if args.parts:
        if first:
            first = False
        else:
            print('')
        print("Parts:")
        print(divider)
        enumerate_subgroups (rootGroup, 'Parts', postfunc=parts_postfunc)
        print(divider)

    if args.constraints:
        if first:
            first = False
        else:
            print('')
        print("Constraints:")
        print(divider)
        enumerate_subgroups (rootGroup, 'Constraints')
        print(divider)

    if args.node_sets:
        if first:
            first = False
        else:
            print('')
        print("NodeSets:")
        print(divider)
        if 'Sets' in rootGroup.groups:
            enumerate_subgroups (rootGroup.groups['Sets'], 'NodeSets')
        print(divider)

    if args.element_sets:
        if first:
            first = False
        else:
            print('')
        print("ElementSets:")
        print(divider)
        if 'Sets' in rootGroup.groups:
            enumerate_subgroups (rootGroup.groups['Sets'], 'ElementSets')
        print(divider)

    if args.problems:
        if first:
            first = False
        else:
            print('')
        print("Problems:")
        print(divider)
        enumerate_subgroups (rootGroup, 'Problems')
        print(divider)

    def solutions_postfunc (group):
        if "NodeValues" in group.groups:
            subgroup = group.groups["NodeValues"]
            print("  Variables defined on nodes:")
            for varName in subgroup.variables.keys():
                print("    " + varName)
        if "ElementValues" in group.groups:
            subgroup = group.groups["ElementValues"]
            print("  Variables defined on elements:")
            for varName in subgroup.variables.keys():
                print("    " + varName)
            if len(subgroup.groups) > 0:
                print("  Variables defined on gauss points:")
                for subsubgroup in subgroup.groups.values():
                    for varName in subsubgroup.variables.keys():
                        print("    " + varName)

    if args.solutions:
        if first:
            first = False
        else:
            print('')
        print("Solutions:")
        print(divider)
        enumerate_subgroups (rootGroup, 'Solutions', postfunc=solutions_postfunc)
        print(divider)


def main():
    try:
        modelinfo()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
