# n88tools

Command-line utilities implemented in python for manipulating finite element models.

[![Build Status](https://dev.azure.com/babesler/n88/_apis/build/status/Numerics88.n88tools?branchName=master)](https://dev.azure.com/babesler/n88/_build/latest?definitionId=11&branchName=master)
[![Anaconda-Server Badge](https://anaconda.org/numerics88/n88tools/badges/installer/conda.svg)](https://anaconda.org/Numerics88/n88tools/)

n88tools works together with the Faim Finite Element package, but it can also be
used to manipulate certain other kinds of finite element files, as well as certain
micro-CT file formats.

## Documentation

A manual exists for the Faim Finite Element package that includes documentation on n88tools,
and a chapter with tutorials. It can be found at http://numerics88.com/documentation/ .

## Compiling and linking

n88tools requires numpy, scipy, VTK (http://www.vtk.org) and
vtkbone (https://github.com/Numerics88/vtkbone).

It can be built and installed with

```sh
python setup.py install
```

On Mac, an extra step may be required
```sh
export CFLAGS="-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk -mmacosx-version-min=10.9 ${CFLAGS}"
export CXXFLAGS="-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk -mmacosx-version-min=10.9 ${CXXFLAGS}"
pip install --no-deps -e .
```
Note that the SDK version - `MacOSX10.15.sdk` - may be different on your machine.


## Authors and Contributors

n88tools is maintained and supported by Numerics88
Solutions Ltd (http://numerics88.com). It was originally developed
by Eric Nodwell (eric.nodwell@numerics88.com) and Steven K. Boyd
(https://bonelab.ucalgary.ca/).

## Licence

n88tools is licensed under a MIT-style open source license. See the file LICENSE.
