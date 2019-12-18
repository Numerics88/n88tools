"""
faim.py

A utility to simply running n88solver, n88derivedfields and n88postfaim.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from .N88ReportedError import N88ReportedError

def faim():

    import os
    import argparse
    import subprocess
    from pipes import quote
    from netCDF4 import Dataset
    import copy
    import numpy


    # ------------------------------------------------------------------------
    # Parse options

    parser = argparse.ArgumentParser (
        prog="faim",
        description=
"""A convenience utility that calls a solver, n88derivedfields and
n88postfaim in sequence on the specified model file.

The most appropriate solver for your input model will be used.
  - If your model contains elastoplastic material definitions, n88solver_spt
    will be used.
  - If your model contains exactly one material array with length equal to
    the number of elements, n88solver_sla will be used.
  - In all other cases, n88solver_slt will be used.

The analysis file produced by n88postfaim will be named after the input
file, with "_analysis.txt" replacing the extension of the input file.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
        )

    parser.add_argument("--engine", "-g", choices=["mt", "nv"],
        help="Set the software engine (default: mt).")
    parser.add_argument ("--threads", "-t",
        help="Set the number of threads for the mt engine.")
    parser.add_argument ("--device", "-d",
        help="""Set the NVidia device for the nv engine.
    A comma-separated list may be specified if multiple devices are to be used.""")
    parser.add_argument ("--precision", choices=["single", "mixed", "double"],
        help="""Set the floating point precision used.""")
    parser.add_argument ("--restart", "-r", action='store_true',
        help="Do not use existing solution as initial value.")
    parser.add_argument ("--convergence_measure", "-m",
        help="Set the convergence measure.", choices=["set","dumax","durms","auto"])
    parser.add_argument ("--convergence_tolerance", "-e",
        help="Set the convergence tolerance.")
    parser.add_argument ("--convergence_window", "-w",
        help="Set the linear convergence window.")
    parser.add_argument ("--plastic_convergence_window", "-W",
        help="Set the plastic convergence window.")
    parser.add_argument ("--maximum_iterations", "-n",
        help="Set maximum number of linear iterations.")
    parser.add_argument ("--maximum_plastic_iterations",
        help="Set maximum number of plastic iterations.")
    parser.add_argument ("--iterations_file",
        help="Specify a file to output all iteration data.")

    parser.add_argument ("--use_coarsen", action='store_true',
        help="Generate a coarsened model if possible to calculate an initial estimate "
             "of the solution. For certain models, this may speed-up finding the final "
             "solution.")
    parser.add_argument ("--material_averaging", choices=["linear", "homminga_density"],
        help="Determine how the material averaging is done. Applies only if the flag "
             "--use_coarsen is used.")

    parser.add_argument ("--no_post", action='store_true',
        help="""Do not run the postprocessor.""")

    parser.add_argument ("--node_sets", "-N",
        help="""Specify the node sets to use for analysis.
    The argument should be a list of node set names separated by commas.
    Specifying this option will override any value in the input file.""")
    parser.add_argument ("--element_sets", "-E",
        help="""Specify the node sets to use for analysis.
    The argument should be a list of element set names separated by commas.
    Specifying this option will override any value in the input file.""")
    parser.add_argument ("--sets", "-s",
        help="A convenience option that sets both node_sets and elements_sets.")
    parser.add_argument ("--rotation_center", "-c",
        help="""Specify the spatial center used for calculation of angular 
    quantities in post-processing. The argument must be given as a triplet of coordinates.
    For bending and torsion tests, a default is taken on the twist axis or the bending axis
    respectively.  Other types of test have no default - for these tests no angular quantities
    will be calculated unless this parameter is specified.""")
    parser.add_argument ("--license_check", "-l", action='store_true',
        help="Print licensing information.")
    parser.add_argument ("--quiet", "-q", action='store_true',
        help="Suppress output to terminal (except for error messages).")
    parser.add_argument ("--verbose", "-v", action='store_true',
        help="Output additional information messages to terminal.")

    parser.add_argument ("input_file", nargs='?', default=None)

    input_args = parser.parse_args()


    # ------------------------------------------------------------------------
    # Check what is in file

    has_existing_solution = False
    has_elastoplastic = False
    has_max_length_material_array = False
    if not input_args.input_file is None:
        try:
            rootGroup = Dataset (input_args.input_file, 'r')
        except RuntimeError:
            raise N88ReportedError ("ERROR: Cannot open file as n88model file.")
        assert (rootGroup.Conventions == "Numerics88/Finite_Element_Model-1.0")

        has_existing_solution = "ActiveSolution" in rootGroup.ncattrs()

        materialGroup = rootGroup.groups['MaterialDefinitions']
        materials = materialGroup.groups.values()
        # Examine material types of all materials: look for any elastoplastic
        for entry in materials:
            if (entry.Type == "VonMisesIsotropic" or
                entry.Type == "MohrCoulombIsotropic" or
                entry.Type == "MaximumPrincipalStrainIsotropic"):
                has_elastoplastic = True
                break
            elif not (entry.Type == "LinearIsotropic" or
                      entry.Type == "LinearOrthotropic" or
                      entry.Type == "LinearAnisotropic"):
                raise N88ReportedError (
"""ERROR: Cannot identify material types.
Select and run solver manually.""")
       # Determine number of elements
        if has_existing_solution:
            activeSolution = rootGroup.ActiveSolution
            activeProblem = rootGroup.groups['Solutions'].groups[activeSolution].Problem
        else:
            activeProblem = rootGroup.ActiveProblem
        activePart = rootGroup.groups['Problems'].groups[activeProblem].Part
        activePartGroup = rootGroup.groups['Parts'].groups[activePart]
        elementNumbers = activePartGroup.groups['Elements'].groups['Hexahedrons'].variables['ElementNumber']
        number_elements = elementNumbers.shape[0]
        # Now check if we can use n88solver_sla.
        # Only need to check length of material array is there is exactly 1 material.
        # The test for a material array is if there are any variables defined.
        if (not has_elastoplastic and
            len(materials) == 1 and
            len(materials[0].variables) > 0):
            # Compare length of first variable of first material to number_elements
            if (materials[0].variables.values()[0].shape[0] == number_elements):
                has_max_length_material_array = True
        rootGroup.close()

    if not input_args.quiet:
        sys.stdout.write ("Model has %d elements.\n" % number_elements)
        if has_elastoplastic:
            sys.stdout.write ("Model contains elastoplastic material definitions.\n")
        else:
            sys.stdout.write ("Model contains only linear material definitions.\n")
        if has_max_length_material_array:
            sys.stdout.write ("Model contains a single material array of length equal to number of elements.\n")
        if has_existing_solution:
            sys.stdout.write ("Model contains existing solution.\n")

    # ------------------------------------------------------------------------
    # Run n88coarsen (and solve) to generate approximate starting solution

    use_coarsen = input_args.use_coarsen
    if has_elastoplastic:
        use_coarsen = False
    if has_existing_solution and not input_args.restart:
        use_coarsen = False
    if input_args.license_check:
        use_coarsen = False
    if input_args.input_file is None:
        use_coarsen = False
    if number_elements < 1000:
        use_coarsen = False
    if use_coarsen:
        if not input_args.quiet:
            sys.stdout.write ("\nGenerating approximate initial solution for starting value.\n")
        coarsened_file = os.path.splitext(input_args.input_file)[0] + "_coarse.n88model"
        call_args = ["n88coarsen"]
        # arguments with a parameter
        arg_list = ["material_averaging"]
        for a in arg_list:
            if not (input_args.__dict__[a] is None):
                call_args += ["--%s=%s" % (a, input_args.__dict__[a])]
        call_args += [input_args.input_file, coarsened_file]
        if not input_args.quiet:
            sys.stdout.write ("Running: %s\n" %(" ".join(map(quote,call_args))))
        sys.stdout.flush()
        p = subprocess.call(call_args)
        if p != 0:
            raise N88ReportedError ("n88coarsen returned error code %d" % p,
                                    value=p)
        solver = "n88solver_slt"
        call_args = [solver]
        # arguments with a parameter
        arg_list = ["engine", "threads", "device", "precision",
                    "convergence_measure",
                    "maximum_iterations", 
                    "convergence_tolerance",
                    "convergence_window"]
        for a in arg_list:
            if not (input_args.__dict__[a] is None):
                call_args += ["--%s=%s" % (a, input_args.__dict__[a])]
        # arguments without a parameter
        arg_list = ["quiet"]
        for a in arg_list:
            if input_args.__dict__[a]:
                call_args += ["--%s" % a]
        call_args += [coarsened_file]
        if not input_args.quiet:
            sys.stdout.write ("\nRunning solver: %s\n" %(" ".join(map(quote,call_args))))
        sys.stdout.flush()
        p = subprocess.call(call_args)
        if p != 0:
            raise N88ReportedError ("Solver returned error code %d" % p,
                                    value=p)
        call_args = ["n88interpolatesolution", input_args.input_file, coarsened_file]
        if not input_args.quiet:
            sys.stdout.write ("\nCopying approximate solution to original model:\n%s\n" %
                (" ".join(map(quote,call_args))))
        sys.stdout.flush()
        p = subprocess.call(call_args)
        if p != 0:
            raise N88ReportedError ("n88interpolatesolution returned error code %d" % p,
                                    value=p)

    # ------------------------------------------------------------------------
    # Run n88solver

    if has_elastoplastic:
        solver = "n88solver_spt"
    else:
        if has_max_length_material_array:
            solver = "n88solver_sla"
        else:
            solver = "n88solver_slt"

    if not input_args.quiet:
        sys.stdout.write ("\nChoosing %s.\n" % solver)
    call_args = [solver]
    # arguments with a parameter
    arg_list = ["engine", "threads", "device", "precision",
                "convergence_measure",
                "maximum_iterations", "maximum_plastic_iterations",
                "convergence_tolerance",
                "convergence_window", "plastic_convergence_window",
                "iterations_file"]
    for a in arg_list:
        if not (input_args.__dict__[a] is None):
            call_args += ["--%s=%s" % (a, input_args.__dict__[a])]
    # arguments without a parameter
    arg_list = ["license_check", "quiet"]
    if not use_coarsen:
        arg_list += ["restart"]
    for a in arg_list:
        if input_args.__dict__[a]:
            call_args += ["--%s" % a]
    if not input_args.input_file is None:
        call_args += [input_args.input_file]
    if not input_args.quiet:
        sys.stdout.write ("Running solver: %s\n" %(" ".join(map(quote,call_args))))
    sys.stdout.flush()
    p = subprocess.call(call_args)
    if p != 0:
        raise N88ReportedError ("Solver returned error code %d" % p,
                                value=p)

    if input_args.license_check:
        return


    # ------------------------------------------------------------------------
    # Run n88derivedfields

    if not input_args.quiet:
        sys.stdout.write("\n")

    call_args = ["n88derivedfields"]
    # arguments with a parameter
    arg_list = []
    for a in arg_list:
        if not (input_args.__dict__[a] is None):
            call_args += ["--%s" % a, input_args.__dict__[a]]
    # arguments without a parameter
    arg_list = ["quiet"]
    for a in arg_list:
        if input_args.__dict__[a]:
            call_args += ["--%s" % a]
    call_args += [input_args.input_file]
    if not input_args.quiet:
        sys.stdout.write ("Running field calculator: %s\n" %(" ".join(map(quote,call_args))))
    sys.stdout.flush()
    p = subprocess.call(call_args)
    if p != 0:
        raise N88ReportedError ("Field calculator returned error code %d" % p,
                                value=p)


    # ------------------------------------------------------------------------
    # Run n88postfaim

    if input_args.no_post:
        return

    if not input_args.quiet:
        sys.stdout.write ("\n")

    call_args = ["n88postfaim"]
    # arguments with a parameter
    arg_list = ["node_sets", "element_sets", "sets", "rotation_center"]
    for a in arg_list:
        if not (input_args.__dict__[a] is None):
            call_args += ["--%s" % a, input_args.__dict__[a]]
    analysis_file = os.path.splitext(input_args.input_file)[0] + "_analysis.txt"
    call_args += ["--output_file", analysis_file, input_args.input_file]
    if not input_args.quiet:
        sys.stdout.write ("Running analysis tool: %s\n" %(" ".join(map(quote,call_args))))
    sys.stdout.flush()
    if sys.platform == "win32":
        p = subprocess.call(call_args, shell=True)
    else:
        p = subprocess.call(call_args)
    if p != 0:
        raise N88ReportedError ("Analysis tool returned error code %d" % p,
                                value=p)


def main():
    try:
        faim()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
