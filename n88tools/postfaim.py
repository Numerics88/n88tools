"""
postfaim.py

Generate tables of standard post-processing quantities.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

from __future__ import division
import sys
from .N88ReportedError import N88ReportedError
import numpy
from numpy.core import *


def postfaim():

    import os
    import argparse
    from collections import OrderedDict
    import warnings
    import numbers
    from netCDF4 import Dataset

    global table_count

    # ------------------------------------------------------------------------
    # Parse options

    parser = argparse.ArgumentParser (
        description="Generate tables of standard post-processing quantities.",
        )

    # parser.add_argument ("--version", "-v", action="store_true",
    #     help="Print version information.")

    parser.add_argument ("--output_file", "-o",
        help="Output file. If not specified, output will go to STDOUT.")

    parser.add_argument ("--node_sets", "-N",
        help="""Specify the node sets to use for analysis.
    The argument should be a list of node set names,
    separated by commas.""")

    parser.add_argument ("--element_sets", "-E",
        help="""Specify the element sets to use for analysis.
    The argument should be a list of element set names,
    separated by commas. Each node set must be matched
    by a corresponding set of elements.""")

    parser.add_argument ("--sets", "-s",
        help="""A convenience option that sets both node_sets and elements_sets.
    This is only useful if corresponding node and element sets are identically named.""")
                                 
    parser.add_argument ("--rotation_center", "-c",
        type = lambda x: numpy.fromstring(x, dtype=float, sep=","),
        help="""Specify the spatial center used for calculation of angular 
    quantities. The argument must be given as a triplet of coordinates.
    The default is read from the n88model. If not specified, either on the
    command line or in the n88model file, no angular quantities will be
    calculated.""")

    parser.add_argument ("model_file")

    args = parser.parse_args()

    if not (args.rotation_center is None):
        if len(args.rotation_center) != 3:
            raise N88ReportedError ("ERROR: rotation_center must be a triplet of values.")

    if args.node_sets == None:
        args.node_sets = args.sets
    if args.element_sets == None:
        args.element_sets = args.sets

    # For now this is hard-wired.
    args.twist_threshold = 1E-6

    if args.output_file == None:
        out = sys.stdout
    else:
        out = open (args.output_file, "wt")



    # ------------------------------------------------------------------------
    # Get information about the n88model file

    # Get the netCDF group handles from the n88model file.

    root = Dataset (args.model_file, "r")
    try:
        activeSolution = root.groups["Solutions"].groups[root.ActiveSolution]
    except KeyError:
        raise N88ReportedError ("No solution found in model file.")
    activeProblem = root.groups["Problems"].groups[activeSolution.Problem]
    activePart = root.groups["Parts"].groups[activeProblem.Part]
    hexahedrons = activePart.groups["Elements"].groups["Hexahedrons"]
    nodeValues = activeSolution.groups["NodeValues"]
    elementValues = activeSolution.groups["ElementValues"]
    materialTable = activePart.groups["MaterialTable"]
    materialDefinitions = root.groups["MaterialDefinitions"]
    try:
        nodeSetsGroup = root.groups["Sets"].groups["NodeSets"]
    except KeyError:
        nodeSetsGroup = None
    try:
        elementSetsGroup = root.groups["Sets"].groups["ElementSets"]
    except KeyError:
        elementSetsGroup = None

    # Determine some constants to do with array sizes

    num = zeros (8,int)
    num[:] = hexahedrons.variables["NodeNumbers"][0] - 1
    coord_var = activePart.variables["NodeCoordinates"]
    spacing = array ([coord_var[num[1],0] - coord_var[num[0],0],
                     coord_var[num[2],1] - coord_var[num[0],1],
                     coord_var[num[4],2] - coord_var[num[0],2]])
    numberOfNodes = coord_var.shape[0]
    dimensionality = coord_var.shape[1]
    del num, coord_var
    numberOfElements = hexahedrons.variables["NodeNumbers"].shape[0]
    numberOfNodesPerElement = hexahedrons.variables["NodeNumbers"].shape[1]
    sizeOfMaterialTable = materialTable.variables["ID"].shape[0]
    nst = 6

    # Material table data
    materialIds = materialTable.variables["ID"][:]
    materialNames = materialTable.variables['MaterialName'][:]

    # Now we need to determine how many of our materials are material arrays.
    materialArrayLength = zeros (sizeOfMaterialTable, dtype=int)
    for m in range(sizeOfMaterialTable):
        material = materialDefinitions.groups[materialNames[m]]
        if len(material.variables) != 0:
            # Has at least one netcdf variable -> is a material array
            materialArrayLength[m] = material.variables.values()[0].shape[0]

    def SplitOnCommaOrSemicolon (s):
        tokens = s.split(",")
        if len(tokens) == 1:
            tokens = s.split(";")
        return tokens

    if args.node_sets != None:
        ppNodeSets = SplitOnCommaOrSemicolon (args.node_sets)
    else:
        try:
            ppNodeSets = SplitOnCommaOrSemicolon (activeProblem.PostProcessingNodeSets)
        except AttributeError:
            ppNodeSets = []
    if args.element_sets != None:
        ppElementSets = SplitOnCommaOrSemicolon (args.element_sets)
    else:
        try:
            ppElementSets = SplitOnCommaOrSemicolon (activeProblem.PostProcessingElementSets)
        except AttributeError:
            ppElementSets = []
    numPPSets = len(ppNodeSets)

    if args.rotation_center is None:
        try:
            args.rotation_center = activeProblem.RotationCenter
        except:
            args.rotation_center = None

    # Array data that it is convenient to read now and keep around.
    # Note that most array data will be read as needed and subsequently
    # discarded for memory efficiency.

    elementMaterials = hexahedrons.variables["MaterialID"][:]


    # ------------------------------------------------------------------------
    # Some constants used for formatting

    width = 76
    table_delimiter = "="*width + "\n"
    section_delimiter = "-"*width + "\n"
    subtable_delimiter = "."*width + "\n"
    any_entry = "%%-32s%%%d" % (width-32)
    text_entry = any_entry + "s\n"
    integer_entry = any_entry + "d\n"


    # ------------------------------------------------------------------------
    # Calculation functions

    def CalculateStats (data, stats):
        if "median" in stats:
            # If sorting we make a copy so we don't modify the original data
            if isinstance(data, numpy.ma.MaskedArray):
                x = sort(data.compressed())
            else:
                x = sort(data)
        else:
            if isinstance(data, numpy.ma.MaskedArray):
                # If a MaskedArray, get a copy of the compressed array - may be faster
                x = data.compressed()
            else:
                x = data
        nan = float('nan')
        n = len(x)
        if n == 0:
            total = 0
            average = nan
            rms = nan
            std_dev = nan
            minimum = nan
            maximum = nan
            skewness = nan
            kurtosis = nan
            median = nan
            perc05 = nan
            perc25 = nan
            perc75 = nan
            perc95 = nan
        else:
            total = sum(x)
            average = total/n
            if "std_dev" in stats:
                std_dev = sqrt(sum((x-average)**2)/n)
            if "rms" in stats:
                rms = sqrt(sum(x**2)/n)
            if "minimum" in stats:
                minimum = min(x)
            if "maximum" in stats:
                maximum = max(x)
            # Might get warning if stddev is zero: ignore these
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")    
                if "skewness" in stats:
                    skewness = sum(((x-average)/std_dev)**3)/n
                if "kurtosis" in stats:
                    kurtosis = sum(((x-average)/std_dev)**4)/n - 3
            median = x[n//2]
            perc05 = x[n//20]
            perc25 = x[n//4]
            perc75 = x[(3*n)//4]
            perc95 = x[(19*n)//20]
        values = OrderedDict()
        for v in stats:
            values[v] = locals()[v]
        return values


    # ------------------------------------------------------------------------
    # Output functions

    def BeginTable (name):
        global table_count
        table_count += 1
        out.write (
            "\n" + table_delimiter
          + "Table %d: %s\n" % (table_count, name)
          + section_delimiter)
        

    def StatsTable (data, stats, labels=None):
        if len(data.shape) == 1:
            num_col = 1
        else:
            num_col = data.shape[1]
        left_width = width - 11*6
        if labels:
            out.write (" "*left_width + ("%11s"*num_col) % tuple(labels) + "\n")
        int_row_format = "%%-%ds" % left_width + "%11d"*num_col + "\n"
        float_row_format = "%%-%ds" % left_width + "%11.3E"*num_col + "\n"
        if num_col == 1:
            values = CalculateStats(data, stats)
            for stat_name, stat_value in values.items():
                if isinstance(stat_value, numbers.Integral):
                    row_format = int_row_format
                else:
                    row_format = float_row_format
                out.write (row_format % (stat_name, stat_value))
        else:
            values = []
            for i in range(num_col):
                values.append (CalculateStats(data[:,i], stats))
            for v in values[0].keys():
                row_values = []
                for i in range(num_col):
                    row_values.append (values[i][v])
                if isinstance(row_values[0], numbers.Integral):
                    row_format = int_row_format
                else:
                    row_format = float_row_format
                out.write (row_format % (tuple([v] + row_values)))
        # Note that the return value is typically not used.
        return values


    def StatsSubtablesByMaterial (data, stats, labels=None):
        if len(data.shape) == 1:
            num_col = 1
        else:
            num_col = data.shape[1]
        if num_col > 1:
            mask = zeros (data.shape, bool)
        for m in range(sizeOfMaterialTable):
            out.write (section_delimiter)
            out.write (text_entry % ("m:", m+1))
            if materialArrayLength[m] <= 1:
                id_field = materialIds[m]
                mask1 = (elementMaterials != materialIds[m])
            else:
                stop_id = materialIds[m] + materialArrayLength[m] -1
                id_field = "%d-%d" % (materialIds[m], stop_id)
                mask1 = logical_or(elementMaterials < materialIds[m],
                                   elementMaterials > stop_id)
            out.write (text_entry % ("Material ID:", id_field))
            out.write (subtable_delimiter)
            if num_col == 1:
                mask = mask1
            else:
                # There must be a better way to do the following!
                for i in range(num_col):
                    mask[:,i] = mask1[:]
            subset = numpy.ma.MaskedArray (data, mask)
            if len(subset) > 0:
                StatsTable (subset, stats, labels)
            else:
                out.write ("No elements.\n")


    def SubtablesBySet (data, stats, labels=None):
        return_values = []
        for n in range(numPPSets):
            if n>0:
                out.write (section_delimiter)
            out.write (text_entry % ("Node set:", n+1))
            out.write (text_entry % ("Name:", ppNodeSets[n]))
            out.write (subtable_delimiter)
            var = nodeSetsGroup.groups[ppNodeSets[n]].variables["NodeNumber"]
            numSetNodes = var.shape[0]
            if numSetNodes > 0:
                set_indices = zeros (numSetNodes, int)
                set_indices[:] = var[:] - 1
                set_data = data[set_indices]
                # Note that the return values are typically not used
                return_values.append (StatsTable (set_data, stats, labels))
            else:
                out.write ("No nodes.\n")
                return_values.append (None)
        return return_values
        

    # ------------------------------------------------------------------------
    # Functions to generate particular tables

    def ModelInputTable():
        BeginTable ("Model Input")
        out.write (text_entry % ("Filename:", os.path.basename(args.model_file)))
        float_entry = any_entry + ".3f\n"
        out.write (
            float_entry % ("Element dim X:", spacing[0])
          + float_entry % ("Element dim Y:", spacing[1])
          + float_entry % ("Element dim Z:", spacing[2]))
        out.write (
            integer_entry % ("Number of elements:", numberOfElements)
          + integer_entry % ("Number of nodes:", numberOfNodes)
          + integer_entry % ("Number of nodes per element:", numberOfNodesPerElement)
          + integer_entry % ("Dimension of problem:", dimensionality))
        out.write (table_delimiter)


    def MaterialsTable():
        BeginTable ("Materials")
        out.write (integer_entry % ("Number of materials:",sizeOfMaterialTable))
        out.write (section_delimiter)
        types = []
        Emaxes = []
        for m in range(sizeOfMaterialTable):
            material = materialDefinitions.groups[materialNames[m]]
            types.append (material.Type)
            if materialArrayLength[m] == 0:
                if types[-1][-9:] == "Isotropic":
                    Emaxes.append ("%.1f" % material.E)
                elif types[-1] == "LinearOrthotropic":
                    Emaxes.append ("%.1f" % max(material.E))
                else:
                    Emaxes.append ("-")
            else:
                Emaxes.append ("-")
        col_width = (width - 34)/2
        row_format = "%%4s %%6s   %%-%ds %%-%ds %%9s %%10s\n" % (col_width, col_width-2)
        out.write (row_format % ("m", "ID", "Name", "Type", "E_ii_max", "Elements"))
        row_format = "%%4d %%8s %%-%ds %%-%ds %%7s %%10d\n" % (col_width, col_width)
        for m,id,name,t,Emax in zip(range(sizeOfMaterialTable),
                                    materialIds,
                                    materialNames,
                                    types,
                                    Emaxes):
            if materialArrayLength[m] == 0:
                count = sum(elementMaterials == id)
                id_field = id
            else:
                stop_id = id + materialArrayLength[m] -1
                id_field = "%d..%d" % (id,stop_id)
                count = sum(logical_and(elementMaterials >= id,
                                        elementMaterials <= stop_id))
            out.write (row_format % (m+1, id_field, name, t, Emax, count))
        out.write (table_delimiter)


    def SetsTable():
        BeginTable ("Post-processing Sets")
        out.write (integer_entry % ("Number of sets:",numPPSets))
        out.write (section_delimiter)
        counts = range(1,numPPSets+1)
        col_width = (width - 28)/2
        row_format = "%%6s   %%-%ds %%9s   %%-%ds %%9s\n" % (col_width-2,col_width-2)
        out.write (row_format % ("n","Node set name","Nodes","Element set name","Elements"))
        row_format = "%%6d   %%-%ds %%7d   %%-%ds %%7d\n" % (col_width,col_width)
        for n,nodeSet,elementSet in zip(counts,ppNodeSets,ppElementSets):
            numNodes = nodeSetsGroup.groups[nodeSet].variables["NodeNumber"].shape[0]
            numElements = elementSetsGroup.groups[elementSet].variables["ElementNumber"].shape[0]
            out.write (row_format %(n,nodeSet,numNodes,elementSet,numElements))
        out.write (table_delimiter)


    def ResidualsTable():
        try:
            residuals = zeros (nodeValues.variables["Residual"].shape, float64)
            residuals[:] = nodeValues.variables["Residual"][:]
        except KeyError:
            return
        BeginTable ("Residuals")
        residuals.shape = (numberOfNodes*dimensionality,)
        residuals = abs(residuals)
        stats = ["maximum", "rms" , "median",
                 "perc05", "perc25", "perc75", "perc95"]
        StatsTable (residuals, stats)
        out.write (table_delimiter)


    def StressTable():
        try:
            stress = zeros (elementValues.variables["Stress"].shape, float64)
            stress[:] = elementValues.variables["Stress"][:]
        except KeyError:
            return
        BeginTable ("Stress")
        stats = ["average", "std_dev", "minimum", "maximum", "skewness", "kurtosis",
                 "median", "perc05", "perc25", "perc75", "perc95"]
        labels=["sigma_xx", "sigma_yy", "sigma_zz", "sigma_yz", "sigma_zx", "sigma_xy"]
        out.write (text_entry % ("m:", "ALL"))
        out.write (text_entry % ("Material ID:", "ALL"))
        out.write (subtable_delimiter)
        StatsTable (stress, stats, labels)
        if sizeOfMaterialTable > 1:
            stats = ["average", "std_dev", "minimum", "maximum",
                     "skewness", "kurtosis", "median"]
            StatsSubtablesByMaterial (stress, stats, labels)
        out.write (table_delimiter)


    def StrainTable():
        try:
            strain = zeros (elementValues.variables["Strain"].shape, float64)
            strain[:] = elementValues.variables["Strain"][:]
        except KeyError:
            return
        BeginTable ("Strain")
        stats = ["average", "std_dev", "minimum", "maximum", "skewness", "kurtosis",
                 "median", "perc05", "perc25", "perc75", "perc95"]
        labels=["epsilon_xx", "epsilon_yy", "epsilon_zz", "gamma_yz", "gamma_zx", "gamma_xy"]
        out.write (text_entry % ("m:", "ALL"))
        out.write (text_entry % ("Material ID:", "ALL"))
        out.write (subtable_delimiter)
        StatsTable (strain, stats, labels)
        if sizeOfMaterialTable > 1:
            stats = ["average", "std_dev", "minimum", "maximum",
                     "skewness", "kurtosis", "median"]
            StatsSubtablesByMaterial (strain, stats, labels)
        out.write (table_delimiter)


    def PlasticStrainTable():
        try:
            plastic_strain = zeros (elementValues.variables["PlasticStrain"].shape, float64)
            plastic_strain[:] = elementValues.variables["PlasticStrain"][:]
        except KeyError:
            return
        BeginTable ("Plastic Strain")
        stats = ["average", "std_dev", "minimum", "maximum", "skewness", "kurtosis",
                 "median", "perc05", "perc25", "perc75", "perc95"]
        labels=["epsilon_xx", "epsilon_yy", "epsilon_zz", "gamma_yz", "gamma_zx", "gamma_xy"]
        out.write (text_entry % ("m:", "ALL"))
        out.write (text_entry % ("Material ID:", "ALL"))
        out.write (subtable_delimiter)
        StatsTable (plastic_strain, stats, labels)
        if sizeOfMaterialTable > 1:
            stats = ["average", "std_dev", "minimum", "maximum",
                     "skewness", "kurtosis", "median"]
            StatsSubtablesByMaterial (plastic_strain, stats, labels)
        out.write (table_delimiter)


    def StrainEnergyDensityTable():
        try:
            sed = zeros (elementValues.variables["StrainEnergyDensity"].shape, float64)
            sed[:] = elementValues.variables["StrainEnergyDensity"][:]
        except KeyError:
            return
        BeginTable ("Strain Energy Density")
        stats = ["average", "std_dev", "minimum", "maximum", "skewness", "kurtosis",
                 "median", "perc05", "perc25", "perc75", "perc95"]
        out.write (text_entry % ("m:", "ALL"))
        out.write (text_entry % ("Material ID:", "ALL"))
        out.write (subtable_delimiter)
        StatsTable (sed, stats)
        if sizeOfMaterialTable > 1:
            stats = ["average", "std_dev", "minimum", "maximum",
                     "skewness", "kurtosis", "median"]
            StatsSubtablesByMaterial (sed, stats)
        out.write (table_delimiter)


    def VonMisesStressTable():
        try:
            vm_stress = zeros (elementValues.variables["VonMisesStress"].shape, float64)
            vm_stress[:] = elementValues.variables["VonMisesStress"][:]
        except KeyError:
            return
        BeginTable ("Von Mises Stress")
        stats = ["average", "std_dev", "minimum", "maximum", "skewness", "kurtosis",
                 "median", "perc05", "perc25", "perc75", "perc95"]
        out.write (text_entry % ("m:", "ALL"))
        out.write (text_entry % ("Material ID:", "ALL"))
        out.write (subtable_delimiter)
        StatsTable (vm_stress, stats)
        if sizeOfMaterialTable > 1:
            stats = ["average", "std_dev", "minimum", "maximum",
                     "skewness", "kurtosis", "median"]
            StatsSubtablesByMaterial (vm_stress, stats)
        out.write (table_delimiter)


    def NodalDisplacementsTable():
        try:
            displacements = zeros (nodeValues.variables["Displacement"].shape, float64)
            displacements[:] = nodeValues.variables["Displacement"][:]
        except KeyError:
            return
        BeginTable ("Nodal Displacements")
        stats = ["average", "std_dev", "minimum", "maximum", "median"]
        values = SubtablesBySet (displacements, stats, labels=["ux","uy","uz"])
        out.write (table_delimiter)


    def NodalForcesTable():
        try:
            forces = zeros (nodeValues.variables["ReactionForce"].shape, float64)
            forces[:] = nodeValues.variables["ReactionForce"][:]
        except KeyError:
            return
        BeginTable ("Nodal Forces")
        stats = ["total", "average", "std_dev", "minimum", "maximum", "median"]
        values = SubtablesBySet (forces, stats, labels=["Fx","Fy","Fz"])
        out.write (table_delimiter)


    def NodalTwistTable():
        try:
            disp = zeros (nodeValues.variables["Displacement"].shape, float64)
            disp[:] = nodeValues.variables["Displacement"][:]
        except KeyError:
            return
        BeginTable ("Nodal Twist")
        center_text = "%.3f %.3f %.3f" % tuple(args.rotation_center)
        row_format = "%%-%ds %%s\n" % (width-1-len(center_text))
        out.write (row_format % ("Twist measured relative to:", center_text))
        float_entry = any_entry + ".1E\n"
        out.write (float_entry % ("Exclusion radius:", args.twist_threshold))
        out.write (section_delimiter)
        p = zeros (activePart.variables["NodeCoordinates"].shape, float64)
        p[:] = activePart.variables["NodeCoordinates"][:]
        p -= args.rotation_center
        p_dp = p + disp
        twist = numpy.arctan2 (numpy.roll(p_dp,1,axis=1), numpy.roll(p_dp,2,axis=1)) \
              - numpy.arctan2 (numpy.roll(p,1,axis=1), numpy.roll(p,2,axis=1))
        # Correct for going over branch cut
        twist += 2*pi * (twist < pi)
        twist -= 2*pi * (twist > pi)
       # Mask out any points (initial or final) that are too close to the rotation center
       # in the rotation plane.
        mask = (numpy.roll(p_dp,1,axis=1)**2 + numpy.roll(p_dp,2,axis=1)**2
                   > args.twist_threshold**2)[:, numpy.newaxis]
        mask *= (numpy.roll(p,1,axis=1)**2 + numpy.roll(p,2,axis=1)**2
                   > args.twist_threshold**2)[:, numpy.newaxis]
        twist_masked = numpy.ma.MaskedArray (twist, numpy.logical_not(mask))
        stats = ["n", "average", "std_dev", "minimum", "maximum", "median"]
        values = SubtablesBySet (twist_masked, stats, labels=["ROTx","ROTy","ROTz"])
        out.write (table_delimiter)


    def NodalTorquesTable():
        try:
            rf = zeros (nodeValues.variables["ReactionForce"].shape, float64)
            rf[:] = nodeValues.variables["ReactionForce"][:]
        except KeyError:
            return
        BeginTable ("Nodal Torques")
        center_text = "%.3f %.3f %.3f" % tuple(args.rotation_center)
        row_format = "%%-%ds %%s\n" % (width-1-len(center_text))
        out.write (row_format % ("Torque measured relative to:", center_text))
        out.write (section_delimiter)
        p = zeros (activePart.variables["NodeCoordinates"].shape, float64)
        p[:] = activePart.variables["NodeCoordinates"][:]
        p -= args.rotation_center
        torque = numpy.roll(rf,1,axis=1) * numpy.roll(p,2,axis=1) \
               - numpy.roll(rf,2,axis=1) * numpy.roll(p,1,axis=1)
        stats = ["total", "average", "std_dev", "minimum", "maximum", "median"]
        values = SubtablesBySet (torque, stats, labels=["Tx","Ty","Tz"])
        out.write (table_delimiter)


    def LoadSharingTable():
        try:
            g_rf = zeros (nodeValues.variables["ReactionForce"].shape, float64)
            g_rf[:] = nodeValues.variables["ReactionForce"][:]
        except KeyError:
            return
        BeginTable ("Load Sharing")
        if not (args.rotation_center is None):
            center_text = "%.4f %.4f %.4f" % tuple(args.rotation_center)
            row_format = "%%-%ds %%s\n" % (width-1-len(center_text))
            out.write (row_format % ("Torque measured relative to:", center_text))
            out.write (section_delimiter)
            g_p = zeros (activePart.variables["NodeCoordinates"].shape, float64)
            g_p[:] = activePart.variables["NodeCoordinates"][:]
            g_p -= args.rotation_center
        g_num = hexahedrons.variables["NodeNumbers"][:] - 1
        # node_material is an array (nn,8) which gives for each node and
        # local element index (0-8) the material index of that element.
        # Is zero for no element there.
        node_material = zeros ((numberOfNodes,8),int)
        indices = numpy.mgrid[0:numberOfElements,0:8]
        node_material[g_num[indices[0],indices[1]],indices[1]] = elementMaterials[indices[0]]
        del indices
        loads = zeros((numPPSets,sizeOfMaterialTable,3), float)
        if not (args.rotation_center is None):
            torques = zeros((numPPSets,sizeOfMaterialTable,3), float)
        materialMask = zeros (g_num.shape, bool)
        for m in range(sizeOfMaterialTable):
            if materialArrayLength[m] <= 1:
                materialMask[:,:] = (elementMaterials == materialIds[m])[:,numpy.newaxis]
            else:
                stop_id = materialIds[m] + materialArrayLength[m] -1
                materialMask[:,:] = (logical_and(
                    elementMaterials >= materialIds[m],
                    elementMaterials <= stop_id))[:,numpy.newaxis]
            nodeListByMaterial = numpy.unique (g_num[numpy.nonzero(materialMask)])
            node_material_reduced = node_material[nodeListByMaterial,:]
            if materialArrayLength[m] <= 1:
                weightsByMaterial = numpy.sum(node_material_reduced == materialIds[m], axis=1) \
                                   / numpy.sum(node_material_reduced > 0, axis=1)
            else:
                weightsByMaterial = numpy.sum(logical_and(node_material_reduced >= materialIds[m],
                                                          node_material_reduced <= stop_id), axis=1) \
                                   / numpy.sum(node_material_reduced > 0, axis=1)
            for n in range(numPPSets):
                nodesInSet = nodeSetsGroup.groups[ppNodeSets[n]].variables["NodeNumber"][:] - 1
                setIndices = numpy.nonzero (numpy.in1d (nodeListByMaterial, nodesInSet, assume_unique=True))
                nodeListByMaterialAndSet = nodeListByMaterial[setIndices]
                weightsByMaterialAndSet = weightsByMaterial[setIndices]
                rf = g_rf[nodeListByMaterialAndSet,:]
                loads[n,m] = numpy.sum (rf * weightsByMaterialAndSet[:,numpy.newaxis], axis=0)
                if not (args.rotation_center is None):
                    p = g_p[nodeListByMaterialAndSet,:]
                    torques[n,m] = numpy.sum (
                        (numpy.roll(rf,1,axis=1) * numpy.roll(p,2,axis=1)
                         - numpy.roll(rf,2,axis=1) * numpy.roll(p,1,axis=1))\
                        * weightsByMaterialAndSet[:,numpy.newaxis],
                        axis = 0)
        col_width = 12
        left_margin = (width - 4*col_width)//2
        heading_format = " "*left_margin + ("%%%ds" % col_width)*4 + "\n"
        table_line = "-"*10
        row_format = " "*left_margin + "%%%ds" % col_width + ("%%%d.4E" % col_width)*3 + "\n"
        for n in range(numPPSets):
            if n > 0:
                out.write (section_delimiter)
            out.write (text_entry % ("Node set:", n+1))
            out.write (text_entry % ("Name:", ppNodeSets[n]))
            out.write (subtable_delimiter)
            out.write (heading_format % ("material", "Fx", "Fy", "Fz"))
            total = zeros(3,float)
            for m in range(sizeOfMaterialTable):
                if materialArrayLength[m] <= 1:
                    id_field = materialIds[m]
                else:
                    stop_id = materialIds[m] + materialArrayLength[m] -1
                    id_field = "%d-%d" % (materialIds[m],stop_id)
                out.write (row_format % ((id_field,) + tuple(loads[n,m])))
                total += loads[n,m]
            out.write (heading_format % ((table_line,)*4))
            out.write (row_format % (("total",) + tuple(total)))
            if not (args.rotation_center is None):
                out.write ("\n");
                out.write (heading_format % ("material", "Tx", "Ty", "Tz"))
                total = zeros(3,float)
                for m in range(sizeOfMaterialTable):
                    if materialArrayLength[m] <= 1:
                        id_field = materialIds[m]
                    else:
                        stop_id = materialIds[m] + materialArrayLength[m] -1
                        id_field = "%d-%d" % (materialIds[m],stop_id)
                    out.write (row_format % ((id_field,) + tuple(torques[n,m])))
                    total += torques[n,m]
                out.write (heading_format % ((table_line,)*4))
                out.write (row_format % (("total",) + tuple(total)))
        out.write (table_delimiter)


    # ------------------------------------------------------------------------
    # Call the desired table functions

    table_count = 0

    ModelInputTable()
    MaterialsTable()
    SetsTable()
    ResidualsTable()
    StrainTable()
    PlasticStrainTable()
    StressTable()
    StrainEnergyDensityTable()
    VonMisesStressTable()
    NodalDisplacementsTable()
    NodalForcesTable()
    if not (args.rotation_center is None):
        NodalTwistTable()
        NodalTorquesTable()
    LoadSharingTable()


def main():
    try:
        postfaim()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
