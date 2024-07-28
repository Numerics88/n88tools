"""
transformations.py

Some elementary geometric manipulations.
"""

from math import cos, sin
from numpy import array

def rotationZ (phi):
    return array (((( cos(phi), -sin(phi), 0 )),
                   (( sin(phi),  cos(phi), 0 )),
                   ((        0,         0, 1 ))))


def rotationX (phi):
    return array ((((        1,         0,  0        )),
                   ((        0,  cos(phi), -sin(phi) )),
                   ((        0,  sin(phi),  cos(phi) ))))


def rotationY (phi):
    return array (((( cos(phi),         0,  sin(phi) )),
                   ((        0,         1,  0        )),
                   ((-sin(phi),         0,  cos(phi) ))))
