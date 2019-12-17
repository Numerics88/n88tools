"""
N88ModelReader

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""

import sys
import numpy
from numpy.core import *
from netCDF4 import Dataset

class N88ModelReader:
    """A class for reading n88model files. It is particularly suited to adapting
    the model data to finite element calculations.
    """

    # Data is stored in n88model files with VTK topology. For FE calculations,
    # we prefer to use the topology of Smith and Griffiths.
    ConvertToSGTopology = True

    def read(self, file_name):
        
        self.file_name = file_name     # Put into class namespace
        rootGroup = Dataset (file_name, 'r')
        rootGroup.set_auto_mask(False)

        assert (rootGroup.Conventions == "Numerics88/Finite_Element_Model-1.0")
        
        self.Dimensionality = len(rootGroup.dimensions['Dimensionality'])
        if 'History' in rootGroup.__dict__:
            self.History = rootGroup.History
        else:
            self.History = None
        if 'Log' in rootGroup.__dict__:        
            self.Log = rootGroup.Log
        else:
            self.Log = None

        if 'ActiveSolution' in rootGroup.__dict__:
            self.ActiveSolution = rootGroup.ActiveSolution
            # Get active problem from active solution
            activeSolutionGroup = rootGroup.groups['Solutions'].groups[self.ActiveSolution]
            self.ActiveProblem = activeSolutionGroup.Problem
        else:
            self.ActiveSolution = None
            self.ActiveProblem = rootGroup.ActiveProblem
        # Get active part from active problem
        activeProblemGroup = rootGroup.groups['Problems'].groups[self.ActiveProblem]
        self.ActivePart = activeProblemGroup.Part

        materialGroup = rootGroup.groups['MaterialDefinitions']
        self.MaterialDefinitions = {}
        for key,entry in materialGroup.groups.items():
            self.MaterialDefinitions[key] = entry.__dict__

        activePartGroup = rootGroup.groups['Parts'].groups[self.ActivePart]

        materialIDs = activePartGroup.groups['MaterialTable'].variables['ID']
        materialNames = activePartGroup.groups['MaterialTable'].variables['MaterialName']
        self.MaterialTable = {}
        for id,name in zip(materialIDs[:], materialNames[:]):
            self.MaterialTable[id] = name

        nodecoord_netCDF = activePartGroup.variables['NodeCoordinates']
        self.NodeCoordinates = nodecoord_netCDF[:]

        element_ids_netCDF = activePartGroup.groups['Elements'].groups['Hexahedrons'].variables['ElementNumber']
        self.ElementIds = element_ids_netCDF[:]
        self.ElementIds -= 1  # Convert to 0-indexed

        node_ids_netCDF = activePartGroup.groups['Elements'].groups['Hexahedrons'].variables['NodeNumbers']
        self.ElementNodeNumbers = node_ids_netCDF[:]
        self.ElementNodeNumbers -= 1  # Convert to 0-indexed
        # Convert to S&G topology
        if self.ConvertToSGTopology:
            key = array((0,4,5,1,2,6,7,3))
            self.ElementNodeNumbers = self.ElementNodeNumbers[:,key]

        materialIDs_netCDF = activePartGroup.groups['Elements'].groups['Hexahedrons'].variables['MaterialID']
        self.ElementMaterialIds = materialIDs_netCDF[:]

        constraintsList = activeProblemGroup.Constraints.split(',')

        constraintsGroup = rootGroup.groups['Constraints']
        self.Constraints = {}
        for key,entry in constraintsGroup.groups.items():
            if key in constraintsList:
                item = {}
                item["Type"] = entry.Type
                item["NodeNumber"] = entry.variables['NodeNumber'][:]
                item["NodeNumber"] -= 1
                item["Sense"] = entry.variables['Sense'][:]
                item["Sense"] -= 1
                item["Value"] = entry.variables['Value'][:]
            self.Constraints[key] = item

        if self.ActiveSolution:
            nodeValuesGroup = activeSolutionGroup.groups['NodeValues']
            self.Displacement = nodeValuesGroup.variables['Displacement'][:]
            if 'ElementValues' in activeSolutionGroup.groups:
                elementValuesGroup = activeSolutionGroup.groups['ElementValues']
                if 'Strain' in elementValuesGroup.variables:
                    self.Strain = elementValuesGroup.variables['Strain'][:]
                if 'PlasticStrain' in elementValuesGroup.variables:
                    self.PlasticStrain = elementValuesGroup.variables['PlasticStrain'][:]
                if 'Stress' in elementValuesGroup.variables:
                    self.Stress = elementValuesGroup.variables['Stress'][:]
                if 'StrainEnergyDensity' in elementValuesGroup.variables:
                    self.StrainEnergyDensity = elementValuesGroup.variables['StrainEnergyDensity'][:]

        rootGroup.close()
        
        # Some derived stuff
        
        if self.ConvertToSGTopology:
            far_corner = 6
        else:
            far_corner = 7
        self.ElementSize = self.NodeCoordinates[self.ElementNodeNumbers[0,far_corner]] \
                         - self.NodeCoordinates[self.ElementNodeNumbers[0,0]]

        # Combine all the Displacement constraints: sort and eliminate duplicates
        self.DisplacementConstraints = {}
        self.DisplacementConstraints["NodeNumber"] = array([], int)
        self.DisplacementConstraints["Sense"] = array([], int)
        self.DisplacementConstraints["Value"] = array([], float64)
        for constraint in self.Constraints.itervalues():
            if constraint["Type"] == "NodeAxisDisplacement":
                self.DisplacementConstraints["NodeNumber"] = numpy.append(self.DisplacementConstraints["NodeNumber"], constraint["NodeNumber"])
                self.DisplacementConstraints["Sense"] = numpy.append(self.DisplacementConstraints["Sense"], constraint["Sense"])
                self.DisplacementConstraints["Value"] = numpy.append(self.DisplacementConstraints["Value"], constraint["Value"])
        if len(self.DisplacementConstraints["NodeNumber"]) > 0:
            key = numpy.lexsort((self.DisplacementConstraints["Sense"],self.DisplacementConstraints["NodeNumber"]))
            mask_unique = numpy.hstack([True,
                                        (numpy.diff(self.DisplacementConstraints["NodeNumber"][key]) != 0) |
                                        (numpy.diff(self.DisplacementConstraints["Sense"][key]) != 0)])
            key = key[mask_unique]
            self.DisplacementConstraints["NodeNumber"] = self.DisplacementConstraints["NodeNumber"][key]
            self.DisplacementConstraints["Sense"] = self.DisplacementConstraints["Sense"][key]
            self.DisplacementConstraints["Value"] = self.DisplacementConstraints["Value"][key]

        # Combine all the Force constraints: sort and eliminate duplicates
        self.ForceConstraints = {}
        self.ForceConstraints["NodeNumber"] = array([], int)
        self.ForceConstraints["Sense"] = array([], int)
        self.ForceConstraints["Value"] = array([], float64)
        for constraint in self.Constraints.itervalues():
            if constraint["Type"] == "NodeAxisForce":
                self.ForceConstraints["NodeNumber"] = numpy.append(self.ForceConstraints["NodeNumber"], constraint["NodeNumber"])
                self.ForceConstraints["Sense"] = numpy.append(self.ForceConstraints["Sense"], constraint["Sense"])
                self.ForceConstraints["Value"] = numpy.append(self.ForceConstraints["Value"], constraint["Value"])
        if len(self.ForceConstraints["NodeNumber"]) > 0:
            key = numpy.lexsort((self.ForceConstraints["Sense"],self.ForceConstraints["NodeNumber"]))
            mask_unique = numpy.hstack([True,
                                        (numpy.diff(self.ForceConstraints["NodeNumber"][key]) != 0) |
                                        (numpy.diff(self.ForceConstraints["Sense"][key]) != 0)])
            key = key[mask_unique]
            self.ForceConstraints["NodeNumber"] = self.ForceConstraints["NodeNumber"][key]
            self.ForceConstraints["Sense"] = self.ForceConstraints["Sense"][key]
            self.ForceConstraints["Value"] = self.ForceConstraints["Value"][key]
