
nemd_default_parameters = {
        'time_step': 0.0005,
        'damping_time': 0.05,
        'pressure': 1.0,
        'pressure_damping': 0.5,
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
lj_parameters = {
        "C":{
            'cutoff': 12.00,       # Angstrom
            'epsilon': 0.00239,    # eV
            'sigma': 3.4           # Angstrom
            },
        "Fe3+":{
            'cutoff': 10.0,
            'epsilon': 0.02166,
            'sigma': 1.91
            },
        "Fe2+":{
            'cutoff': 10.0,
            'epsilon': 0.02106,
            'sigma': 1.91
            },
        "Cl-":{
            'cutoff': 10.0,
            'epsilon': 0.012889,
            'sigma': 3.470
            },
        "O":{
            'cutoff': 10.0,
            'epsilon': 0.006502,
            'sigma': 3.166 
            }
        }


