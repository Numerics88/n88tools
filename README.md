# n88tools

Command-line utilities implemented in python for manipulating finite element models.

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

## Authors and Contributors

n88tools is maintained and supported by Numerics88
Solutions Ltd (http://numerics88.com). It was originally developed
by Eric Nodwell (eric.nodwell@numerics88.com) and Steven K. Boyd
(https://bonelab.ucalgary.ca/).

## Licence

n88tools is licensed under a MIT-style open source license. See the file LICENSE.
