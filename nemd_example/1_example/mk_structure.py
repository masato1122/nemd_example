# -*- coding: utf-8 -*-
import numpy as np
from optparse import OptionParser
import pymatgen as mg
import ase.build.graphene_nanoribbon as build_ribbon

def _make_graphene(n1=3, n2=4, type='armchair'):
    return build_ribbon(n1, n2, type=type, saturated=False, vaccum=0.0)

def main(options):
    
    crystal = _make_graphene(options.n1, options.n2)

    
if __name__ == '__main__':
    parser = OptionParser()
    #parser.add_option("-f", "--filename", dest="filename", type="string",
    #    help="input file name")
    parser.add_option("--n1", dest="n1", type="int",
            default=3, help="number 1")
    parser.add_option("--n2", dest="n2", type="int",
            default=4, help="number 2")

    (options, args) = parser.parse_args()
    #file_check(options.filename)
    main(options)

