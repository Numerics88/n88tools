from setuptools import setup
from setuptools.extension import Extension
import numpy

setup (

    name = "n88tools",
    version = "8.0.0.dev1",
    packages = ['n88tools'],
    ext_modules = [Extension('n88tools.finiteelementfunctions',
                             include_dirs = [numpy.get_include()],
    	                     sources = ['extensions/finiteelementfunctions.cpp'])],

    install_requires = ['numpy',
                        'scipy'],

    entry_points = {
        'console_scripts': [
            'n88coarsen = n88tools.coarsen:main',
            'n88copymodel = n88tools.copymodel:main',
            'n88directmechanics = n88tools.directmechanics:main',
            'n88evaluate = n88tools.evaluate:main',
            'n88extractfields = n88tools.extractfields:main',
            'n88extractsets = n88tools.extractsets:main',
            'faim = n88tools.faim:main',
            'n88interpolatesolution = n88tools.interpolatesolution:main',
            'n88modelinfo = n88tools.modelinfo:main',
            'n88pistoia = n88tools.pistoia:main',
            'n88postfaim = n88tools.postfaim:main',
            'n88tabulate = n88tools.tabulate:main',
        ],
    },

    zip_safe = False,  # force egg to unzip on installation

    author = "Eric Nodwell",
    author_email = "eric.nodwell@numerics88.com",
    description = "Various useful tools for use with n88model finite element models.",
    url = "http://numerics88.com/",

)
