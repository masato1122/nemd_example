#'damping_time': 0.05,
#'pressure_damping': 0.5,

nemd_default_parameters = {
        'time_step': 0.0005,
        'damping_time': 100.00,
        'pressure_damping': 1000.00,
        'pressure': 0.0,
        'read_data': None,
        'restart_file': None,
        'potential_file': 'opt.tersoff',
        'dump_minimization': 'mini.dump',
        'dump_file': None,
        'tcold': 290.,
        'thot': 310.,
        'ntotal_thermo': 100,
        'ntotal_dump': 10,
        'time_npt': None,
        'time_increase': 10.0,
        'time_nemd': None
        }

## S. Zhibo, et. al., RCS Adv., 2018, 8, 38706
## Note that this paper does not give the cutoff parameters
#lj_parameters = {
#        "C":{
#            'cutoff': 12.00,       # Angstrom
#            'epsilon': 0.00239,    # eV
#            'sigma': 3.4           # Angstrom
#            },
#        "Fe3+":{
#            'cutoff': 10.0,
#            'epsilon': 0.02166,
#            'sigma': 1.91
#            },
#        "Fe2+":{
#            'cutoff': 10.0,
#            'epsilon': 0.02106,
#            'sigma': 1.91
#            },
#        "Cl-":{
#            'cutoff': 10.0,
#            'epsilon': 0.012889,
#            'sigma': 3.470
#            },
#        "O":{
#            'cutoff': 10.0,
#            'epsilon': 0.006502,
#            'sigma': 3.166 
#            }
#        }

##
## Carbon-carbon
## L.A. Girifalco, M. Hodak, and R.S.Lee, Phys. Rev. B 62, 19, 13104 (2000).
##
## A.K. Rappe, etal., J. Am. Chem. Soc. 114, 25, 10025 (1992).
##
#        "C":{
#            'cutoff' : 10.00,        # Angstrom
#            'epsilon': 0.0046616,    # eV
#            'sigma'  : 3.35          # Angstrom
#            },
lj_parameters_rappe = {
        "C":{
            'cutoff' : 10.00,         # Angstrom
            'epsilon': 0.0023967,     # eV
            'sigma'  : 3.4147762      # Angstrom
            },
        "Cl":{
            'cutoff' : 10.00,
            'epsilon': 0.0098434,
            'sigma'  : 3.51638
            },
        "Fe":{
            'cutoff' : 10.0,
            'epsilon': 0.0005637,
            'sigma'  : 2.59430 
            }
        }


