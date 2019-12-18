"""
tabulate.py

A tool to extract numerical values from n88postfaim output files and
collate in a table.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from .N88ReportedError import N88ReportedError

def tabulate():
    from .tables import get_tables, lookups
    import platform
    import argparse
    import collections
    import re

    parser = argparse.ArgumentParser (
        prog="n88tabulate",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Extract and tabulate values from n88postfaim output files.
The result is suitable for importing into a spreadsheet.""",
        epilog="""
For any variables not found in the analysis file, a dash ('-') will be
inserted into the output table.

examples:

  In the following example, every possible variable will be extracted from
  every file ending in "_analysis.txt" in the directory. This will result
  in a very large number of values. The output will be written to the file
  summary.txt.

    $ n88tabulate -H -o summary.txt *_analysis.txt

  In the following example, we request the total forces along the z direction
  on the first two node sets for two different analysis files.
 
    $ n88tabulate -H -V "filename,fz_ns1,fz_ns2" test25a_analysis.txt test42a_analysis.txt 
    filename	fz_ns1	fz_ns2
    test25a_uniaxial.n88model	-0.1019E+02	0.1019E+02
    test42a_uniaxial.n88model	-0.4641E+02	0.4641E+02

  The following example is exactly the same, but the output goes to a file.

    $ n88tabulate -H -V "filename,fz_ns1,fz_ns2" -o summary.txt test25a_analysis.txt test42a_analysis.txt 

  Now suppose we want the same selection of variables on a different analysis
  file, we could do the following. Notice that we are getting the list
  of variables from the existing file summary.txt.

    $ n88tabulate -H --from summary.txt test99a_analysis.txt 
    filename	fz_ns1	fz_ns2
    test99a_uniaxial.n88model	-0.6442E+02	0.6442E+02

  In the following example, we are interested in the strain energy density,
  and request some statistical values. We also use a comma as the delimiter.

    $ n88tabulate -H -d ',' -V "sed_avg,sed_stddev,sed_skew,sed_kurt" analysis.txt
    sed_avg,sed_stddev,sed_skew,sed_kurt
    0.2990E+00,0.1127E+00,-0.1512E+01,0.8999E+00

variables:

In the variable list, certain variables can take one or more numeric indices.
This is denoted with %. For example, for id_mat%m, actual variable names
are id_mat1 for the first defined material, id_mat2 for the second defined
material, and so on. Similarly for dx_avg_ns%n, actual variables names are
dx_avg_ns1 for the first node set, dx_avg_ns2 for the second node set, and
so on.

variables corresponding to values in table "Model Input":

    variable name       value in analysis file table
    ------------------------------------------------
    filename            Filename
    dx_els              Element dim X
    dy_els              Element dim Y
    dz_els              Element dim Z
    num_els             Number of elements
    num_nodes           Number of nodes
    num_nodes_per_els   Number of nodes per element
    dim_of_problem      Dimension of problem

variables corresponding to values in table "Materials":

    variable name    value in analysis file table
    ---------------------------------------------
    num_mats         Number of materials
    id_mat%m         The material ID of the nth defined material
    count_mat%m      Number of elements for the nth defined material

variables corresponding to values in table "Post-processing Sets":

    variable name    value in analysis file table
    ---------------------------------------------
    num_pp_sets      The number of sets used in post-processing.

variables corresponding to values in table "Strain":

    Variable names are listed below. See table "Strain Energy Density"
    for additional description. 'eps' standard for epsilon and 'gam'
    for 'gamma'.
    ------------------------------------------------------------------
    eps_ave_xx      eps_ave_yy      eps_ave_zz
    eps_kurt_xx     eps_kurt_yy     eps_kurt_zz
    eps_max_xx      eps_max_yy      eps_max_zz
    eps_median_xx   eps_median_yy   eps_median_zz
    eps_min_xx      eps_min_yy      eps_min_zz
    eps_perc05_xx   eps_perc05_yy   eps_perc05_zz
    eps_perc25_xx   eps_perc25_yy   eps_perc25_zz
    eps_perc75_xx   eps_perc75_yy   eps_perc75_zz
    eps_perc95_xx   eps_perc95_yy   eps_perc95_zz
    eps_skew_xx     eps_skew_yy     eps_skew_zz
    eps_stddev_xx   eps_stddev_yy   eps_stddev_zz
    gam_ave_xy      gam_ave_yz      gam_ave_zx
    gam_kurt_xy     gam_kurt_yz     gam_kurt_zx
    gam_max_xy      gam_max_yz      gam_max_zx
    gam_median_xy   gam_median_yz   gam_median_zx
    gam_min_xy      gam_min_yz      gam_min_zx
    gam_perc05_xy   gam_perc05_yz   gam_perc05_zx
    gam_perc25_xy   gam_perc25_yz   gam_perc25_zx
    gam_perc75_xy   gam_perc75_yz   gam_perc75_zx
    gam_perc95_xy   gam_perc95_yz   gam_perc95_zx
    gam_skew_xy     gam_skew_yz     gam_skew_zx
    gam_stddev_xy   gam_stddev_yz   gam_stddev_zx

variables corresponding to values in table "Stress":

    Variable names are listed below. See table "Strain Energy Density"
    for additional description. 'sig' stands for 'sigma'
    ------------------------------------------------------------------
    sig_ave_xx      sig_ave_xy      sig_ave_yy
    sig_ave_yz      sig_ave_zx      sig_ave_zz
    sig_kurt_xx     sig_kurt_xy     sig_kurt_yy
    sig_kurt_yz     sig_kurt_zx     sig_kurt_zz
    sig_max_xx      sig_max_xy      sig_max_yy
    sig_max_yz      sig_max_zx      sig_max_zz
    sig_median_xx   sig_median_xy   sig_median_yy
    sig_median_yz   sig_median_zx   sig_median_zz
    sig_min_xx      sig_min_xy      sig_min_yy
    sig_min_yz      sig_min_zx      sig_min_zz
    sig_perc05_xx   sig_perc05_xy   sig_perc05_yy
    sig_perc05_yz   sig_perc05_zx   sig_perc05_zz
    sig_perc25_xx   sig_perc25_xy   sig_perc25_yy
    sig_perc25_yz   sig_perc25_zx   sig_perc25_zz
    sig_perc75_xx   sig_perc75_xy   sig_perc75_yy
    sig_perc75_yz   sig_perc75_zx   sig_perc75_zz
    sig_perc95_xx   sig_perc95_xy   sig_perc95_yy
    sig_perc95_yz   sig_perc95_zx   sig_perc95_zz
    sig_skew_xx     sig_skew_xy     sig_skew_yy
    sig_skew_yz     sig_skew_zx     sig_skew_zz
    sig_stddev_xx   sig_stddev_xy   sig_stddev_yy
    sig_stddev_yz   sig_stddev_zx   sig_stddev_zz

variables corresponding to values in table "Strain Energy Density":

    variable name    value in analysis file table
    ---------------------------------------------
    sed_avg          average (all materials)
    sed_stddev       std_dev (all materials)
    sed_skew         skewness (all materials)
    sed_kurt         kurtosis (all materials)
    sed_min          minimum (all materials)
    sed_max          maximum (all materials)
    sed_median       median (all materials)
    sed_perc05       5th percentile (all materials)
    sed_perc25       25th percentile (all materials)
    sed_perc75       75th percentile (all materials)
    sed_perc95       95th percentile (all materials)
    sed_avg_mat%m    average over the nth defined material
    sed_stddev_mat%m sed_dev over the nth defined material
    sed_skew_mat%m   skewness over the nth defined material
    sed_kurt_mat%m   kurtosis over the nth defined material
    sed_min_mat%m    minimum over the nth defined material
    sed_max_mat%m    maximum over the nth defined material
    sed_median_mat%m median over the nth defined material
    sed_perc05_mat%m 5th percentile over the nth defined material
    sed_perc25_mat%m 25th percentile over the nth defined material
    sed_perc75_mat%m 75th percentile over the nth defined material
    sed_perc95_mat%m 95th percentile over the nth defined material

variables corresponding to values in table "Von Mises Stress":

    variable name    value in analysis file table
    ---------------------------------------------
    svm_avg          average (all materials)
    svm_stddev       std_dev (all materials)
    svm_skew         skewness (all materials)
    svm_kurt         kurtosis (all materials)
    svm_min          minimum (all materials)
    svm_max          maximum (all materials)
    svm_median       median (all materials)
    svm_perc05       5th percentile (all materials)
    svm_perc25       25th percentile (all materials)
    svm_perc75       75th percentile (all materials)
    svm_perc95       95th percentile (all materials)
    svm_avg_mat%m    average over the nth defined material
    svm_stddev_mat%m svm_dev over the nth defined material
    svm_skew_mat%m   skewness over the nth defined material
    svm_kurt_mat%m   kurtosis over the nth defined material
    svm_min_mat%m    minimum over the nth defined material
    svm_max_mat%m    maximum over the nth defined material
    svm_median_mat%m median over the nth defined material
    svm_perc05_mat%m 5th percentile over the nth defined material
    svm_perc25_mat%m 25th percentile over the nth defined material
    svm_perc75_mat%m 75th percentile over the nth defined material
    svm_perc95_mat%m 95th percentile over the nth defined material

variables corresponding to values in table "Nodal Displacements":

    variable name    value in analysis file table
    ---------------------------------------------
    dx_avg_ns%n      The average displacement in the x direction over all
                     nodes in the nth node set.
    dx_stddev_ns%n   The standard deviation of the displacement in the x
                     direction over all nodes in the nth node set.
    dx_min_ns%n      The minimum displacement in the x direction over all
                     nodes in the nth node set.
    dx_max_ns%n      The maximum displacement in the x direction over all
                     nodes in the nth node set.
    dx_median_ns%n   The median displacement in the x direction over all
                     nodes in the nth node set.
    dy_avg_ns%n      The average displacement in the y direction over all
                     nodes in the nth node set.
    dy_stddev_ns%n   The standard deviation of the displacement in the y
                     direction over all nodes in the nth node set.
    dy_min_ns%n      The minimum displacement in the y direction over all
                     nodes in the nth node set.
    dy_max_ns%n      The maximum displacement in the y direction over all
                     nodes in the nth node set.
    dy_median_ns%n   The median displacement in the y direction over all
                     nodes in the nth node set.
    dz_avg_ns%n      The average displacement in the z direction over all
                     nodes in the nth node set.
    dz_stddev_ns%n   The standard deviation of the displacement in the z
                     direction over all nodes in the nth node set.
    dz_min_ns%n      The minimum displacement in the z direction over all
                     nodes in the nth node set.
    dz_max_ns%n      The maximum displacement in the z direction over all
                     nodes in the nth node set.
    dz_median_ns%n   The median displacement in the z direction over all
                     nodes in the nth node set.

variables corresponding to values in table "Nodal Forces":

    variable name    value in analysis file table
    ---------------------------------------------
    fx_ns%n          The total force in the x direction over all nodes in
                     the nth node set.
    fx_stddev_ns%n   The standard deviation of the force in the x
                     direction over all nodes in the nth node set.
    fx_min_ns%n      The minimum force in the x direction over all
                     nodes in the nth node set.
    fx_max_ns%n      The maximum force in the x direction over all
                     nodes in the nth node set.
    fx_median_ns%n   The median force in the x direction over all
                     nodes in the nth node set.
    fy_ns%n          The total force in the y direction over all nodes in
                     the nth node set.
    fy_stddev_ns%n   The standard deviation of the force in the y
                     direction over all nodes in the nth node set.
    fy_min_ns%n      The minimum force in the y direction over all
                     nodes in the nth node set.
    fy_max_ns%n      The maximum force in the y direction over all
                     nodes in the nth node set.
    fy_median_ns%n   The median force in the y direction over all
                     nodes in the nth node set.
    fz_ns%n          The total force in the z direction over all nodes in
                     the nth node set.
    fz_stddev_ns%n   The standard deviation of the force in the z
                     direction over all nodes in the nth node set.
    fz_min_ns%n      The minimum force in the z direction over all
                     nodes in the nth node set.
    fz_max_ns%n      The maximum force in the z direction over all
                     nodes in the nth node set.
    fz_median_ns%n   The median force in the z direction over all
                     nodes in the nth node set.

variables corresponding to values in table "Load Sharing":

    variable name    value in analysis file table
    ---------------------------------------------
    fx_ns%n_mat%m    The force in the x direction over nodes in the nth node
                     set, summed over nodes belonging to the mth defined
                     material.
    fy_ns%n_mat%m    The force in the y direction over nodes in the nth node
                     set, summed over nodes belonging to the mth defined
                     material.
    fz_ns%n_mat%m    The force in the z direction over nodes in the nth node
                     set, summed over nodes belonging to the mth defined
                     material.

variables corresponding to values in table "Pistoia Failure Load Estimate":

NOTE: The Pistoia Failure Load Estimate table is no longer generated by
      n88postfaim, but may be generated by the tool n88pistoia if desired.
      The following variables are by default not included, unless a specific
      list which includes them is specified.

    variable name    value in analysis file table
    ---------------------------------------------
    pis_stiffx       The stiffness in the x direction.
    pis_stiffy       The stiffness in the y direction.
    pis_stiffz       The stiffness in the z direction.
    pis_fx_fail      The estimated failure load in the z direction.
    pis_fy_fail      The estimated failure load in the z direction.
    pis_fz_fail      The estimated failure load in the z direction.
""")

    parser.add_argument ("--variables", "-V",
        help="""A list of the variables to extract from the input files. Separate
    variable names with commas. See below for a list of valid variable names.
    If not specified all possible variables will be selected (this makes for a
    very large table).""")

    parser.add_argument ("--from", "-f",
        help="""Obtain the list of variables from the first line of a text file.
    The variable names may be separated by any kind of delimiter (white space
    or commas). Note that an output file (from n88tabulate, generated with the
    --header option) can be used as a --from argument, in which case the same
    selection of variables will be used.""")

    parser.add_argument ("--header", "-H",
        action="store_true",
        help="""Print a header line first. May be used even if no input files are
    specified.""")

    parser.add_argument ("--delimiter", "-d",
        default = "\t",
        help="""Delimiter character to separate columns. Default is a tab
    ('\\t').""")

    parser.add_argument ("--output_file", "-o",
        help="Output file. If not specified, output will go to STDOUT.")

    parser.add_argument ("input_files", nargs='*',
       help="""n88postfaim output files to process. Any number may be specified,
    and wildcard expansion of * and ? is performed on systems where the shell does
    not do this.""")

    args = parser.parse_args()


    if vars(args)["from"] != None:
        keys = re.findall (r'\w+', open(vars(args)["from"], "rt").readline())
    elif args.variables == None:
        # Variables not explicitly selected; so select them all.
        # This will require opening the first input file to check for the
        # number of materials and number of sets.
        tables = get_tables (args.input_files[0])
        num_mats = int(lookups["num_mats"]())
        num_pp_sets = int(lookups["num_pp_sets"]())
        keys = collections.OrderedDict()
        for key,value in lookups.iteritems():
            key_added = False
            if not key_added:
                s = re.search ('\A(\D+)_ns%n_mat%m\Z', key)
                if s:
                    # Don't break down by material if a single material
                    if num_mats > 1:
                        for n in range(num_pp_sets):
                            for m in range(num_mats):
                                keys["%s_ns%d_mat%d" % (s.group(1),n+1,m+1)] = value
                    key_added = True
            if not key_added:
                s = re.search ('\A(\D+)_mat%m\Z', key)
                if s:
                    if num_mats > 1:
                        # Don't break down by material if a single material
                        for m in range(num_mats):
                            keys["%s_mat%d" % (s.group(1),m+1)] = value
                    key_added = True
            if not key_added:
                s = re.search ('\A(\D+)_ns%n\Z', key)
                if s:
                    for n in range(num_pp_sets):
                        keys["%s_ns%d" % (s.group(1),n+1)] = value
                    key_added = True
            if not key_added:
                # From version 7, skip Pistoia unless specifically requested.
                if not (key in ["pis_stiffx", "pis_stiffy", "pis_stiffz", "pis_fx_fail", "pis_fy_fail", "pis_fz_fail"]):
                    keys[key] = value
    else:
        keys = args.variables.split(",")

    if args.output_file == None:
        out = sys.stdout
    else:
        out = open (args.output_file, "wt")

    if args.header:
        out.write (args.delimiter.join(keys) + "\n")

    # On Windows, we have to do wildcard expansion ourselves; on Unix systems
    # the shell has already done this for us.
    if platform.system() == "Windows":
        import glob
        input_files = []
        for input_arg in args.input_files:
            input_files += glob.glob(input_arg)
    else:
        input_files = args.input_files

    for input_file in input_files:
        tables = get_tables (input_file)
        count = 0
        for key in keys:
            value = None
            try:
                m = re.search ('\A(\D+)_ns(\d+)_mat(\d+)\Z', key)
                if m:
                    value = lookups["%s_ns%%n_mat%%m" % m.group(1)](int(m.group(2)),int(m.group(3)))
            except:
                pass
            if value == None:
                try:
                    m = re.search ('\A(\D+)_mat(\d+)\Z', key)
                    if m:
                        value = lookups["%s_mat%%m" % m.group(1)](int(m.group(2)))
                except:
                    pass
            if value == None:
                try:
                    m = re.search ('\A(\D+)_ns(\d+)\Z', key)
                    if m:
                        value = lookups["%s_ns%%n" % m.group(1)](int(m.group(2)))
                except:
                    pass
            if value == None:
                try:
                    value = lookups[key]()
                except:
                    pass
            if value == None:
                value = '-'
            if count:
                out.write (args.delimiter)
            out.write (value)
            count += 1
        out.write ("\n")


def main():
    try:
        tabulate()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
