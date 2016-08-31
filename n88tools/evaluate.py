"""
evaluate.py

A tool to evaluate the quality of solutions.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from N88ReportedError import N88ReportedError
from numpy.core import *
from finiteelement import *

def evaluate():

    import time
    import argparse
    import numpy
    from numpy import column_stack
    from N88ModelReader import N88ModelReader

    # -------------------------------------------------------------------------
    #    Routines for standard tasks

    def log(msg, *additionalLines):
        """Print message with time stamp.

        The first argument is printed with the a time stamp.
        Subsequent arguments are printed one to a line without a timestamp.
        """
        print "%8.2f %s" % (time.time()-start_time, msg)
        for line in additionalLines:
            print " " * 9 + line
        sys.stdout.flush()

    start_time = time.time()


    # ---------------------------------------------------------------------------
    # Parse options and load python modules.

    parser = argparse.ArgumentParser (
        prog="n88evaluate",
        description="""A tool to evaluate the quality of solutions."""
        )

    parser.add_argument ("input")
    parser_group = parser.add_mutually_exclusive_group()
    parser_group.add_argument("-s", "--sparse",
        action="store_true", dest="useSparse", default=True,
        help="Use sparse matrices (default). Requires scipy.")
    parser_group.add_argument("-d", "--dense",
        action="store_false", dest="useSparse", default=True,
        help="Do not use sparse matrices.")

    args = parser.parse_args()

    log ("Using %s matrices." % ["dense", "sparse"][args.useSparse])
    if args.useSparse:
        try:
            import scipy
            from scipy import sparse
            from scipy.sparse import linalg
        except:
            raise N88ReportedError (
"""ERROR: Unable to import scipy.
You can try to rerun with the -d flag to avoid loading scipy, but much more
memory will be required.""")
    else:
        from numpy import linalg


    # -------------------------------------------------------------------------
    # More useful functions

    def rms (x):
        return sqrt(sum(x**2)/len(x))


    # -------------------------------------------------------------------------
    # Main program

    results = ""

    reader_format = "N88 Model"
    log ("Reading %s file: %s" % (reader_format, args.input))
    reader = N88ModelReader()
    reader.read(args.input)

    log ("Reading solution displacements.")
    x = reader.Displacement
    nn = x.shape[0]

    log ("Reading boundary conditions.")
    bc_index = 3*reader.DisplacementConstraints["NodeNumber"] \
               + reader.DisplacementConstraints["Sense"]
    bc_value = reader.DisplacementConstraints["Value"]
    bc_err = x.flatten()[bc_index] - bc_value
    def stats (e):
        return (max(e), rms(e))
    results += """
    Analysis of solution displacements at boundary conditions:
        max err           : %.2E
        rms err           : %.2E
    """ % stats(bc_err)

    log ("Generating local stiffness matrices.")
    mat_D, mat_km = generate_local_stiffness(reader)
    g_mat = reader.ElementMaterialIds
    id_key, mat_D, mat_km, g_mat = renumber_material_ids(mat_D, mat_km, g_mat)

    log ("Reading applied forces.")
    b = generate_force_terms(
            reader.ForceConstraints["NodeNumber"],
            reader.ForceConstraints["Sense"],
            reader.ForceConstraints["Value"],
            nn)

    g_num = reader.ElementNodeNumbers

    nels = len(reader.ElementNodeNumbers)
    if 'PlasticStrain' in reader.__dict__:

        log ("Calculating strain.")
        a = reader.ElementSize
        strain = calculate_strain (x, g_num, a)
        plastic_strain = reader.PlasticStrain
        elastic_strain = strain - plastic_strain
        del strain  # Don't need total strain after this

        log ("Calculating stress.")
        stress = calculate_stress (elastic_strain, mat_D, g_mat)

        log ("Calculating body loads.")
        body_load = calculate_body_load (plastic_strain, mat_D, g_mat, g_num, a, nn)
        b += body_load.flatten()

        log ("Calculating yield function.")
        f = calculate_yield_function (stress, reader.MaterialTable, reader.MaterialDefinitions, id_key, g_mat)

        # Now do some statistics on yield function.
        # max_stress should probably be max principal stress, but this is
        # close enough for comparing the values of F with.
        max_stress = max(abs(stress.flatten()))
        yielded = any(plastic_strain!=0, axis=1)
        # Shouldn't be necessary to compress, but doesn't seem to work otherwise.
        f_unyielded_and_posf = numpy.ma.masked_array(f, mask=logical_or(yielded,f<0)).compressed()
        def stats (e, max_s, no_min=False):
            count = len(e)
            if count > 0:
                min_e = min(e)
                max_e = max(e)
                rms_e = rms(e)
                vals = count, "%.2E" % min_e, "%.2E" % max_e, "%.2E" % rms_e, \
                       "%.2E" % (min_e/max_s), "%.2E" % (max_e/max_s), "%.2E" % (rms_e/max_s)
            else:
                vals = count, '-', '-', '-', '-', '-', '-'
            if no_min:
                return tuple([vals[i] for i in (0,2,3,5,6)])
            else:
                return tuple(vals)
        results += """
    Analysis of yield function F:
        Unyielded elements with F>0:
            count             : %d
            max F             : %s
            rms F             : %s
            max F/max stress  : %s
            rms F/max stress  : %s""" % stats (f_unyielded_and_posf, max_stress, no_min=True)
        f_yielded = numpy.ma.masked_array(f, mask=~yielded).compressed()
        results += """
        Yielded elements:
            count             : %d
            min F             : %s
            max F             : %s
            rms F             : %s
            min F/max stress  : %s
            max F/max stress  : %s
            rms F/max stress  : %s
    """ % stats (f_yielded, max_stress)


    log ("Assembling global stiffness matrix.")
    K = generate_global_stiffness(mat_km, g_mat, g_num, nn, sparse=args.useSparse)

    log ("Calculating forces.")
    b_solution = calculate_forces (K, x, sparse=args.useSparse)

    # We need a mask so that we can avoid evaluating the force errors
    # on degrees of freedom subject to boundary conditions.
    mask = ones(nn*3,int16)
    mask[bc_index] = 0
    b_err = (b_solution - b)*mask

    def stats (e, b):
        max_e = max(e)
        rms_e = rms(e)
        max_b = max(b)
        return max_e, rms_e, max_e/max_b, rms_e/max_b
    results += """
    Analysis of forces (residuals):
        max err           : %.2E
        rms err           : %.2E
        max err/max force : %.2E
        rms err/max force : %.2E
    """ % stats(b_err, b_solution)

    sys.stdout.write(results)


def main():
    try:
        evaluate()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
