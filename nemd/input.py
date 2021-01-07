# -*- coding: utf-8 -*-
import numpy as np
from nemd.default import nemd_default_parameters, lj_parameters_rappe

def write_nemd_inputs(atoms, 
        datafile='data.lammps',
        index_nemd=None, lthermo=None,
        thot=310, tcold=290, 
        time_npt=10., time_increase=None, time_nemd=100., 
        nloop=2, 
        output_minimization=True, 
        charge=False, pppm=1e-5, thermostat='Nose',
        **kwargs):
    
    import copy
    parameters = copy.copy(nemd_default_parameters)
    
    ## modify default values
    for kw in kwargs.keys():
        parameters[kw] = kwargs[kw]
    parameters['time_npt'] = time_npt
    parameters['time_increase'] = time_increase
    parameters['time_nemd'] = time_nemd
    parameters['thot'] = thot
    parameters['tcols'] = tcold
    if output_minimization is None:
        parameters['dump_minimization'] = None
    
    ## Check parameters
    for kw in parameters.keys():
        if kw is None:
            print(" Error: %s is not defined."%(kw))
            exit()
    
    ## NPT
    parameters['read_data'] = datafile
    parameters['restart_file'] = 'restart0.nemd'
    parameters['dump_file'] = 'npt.dump'
    write_nemd_input(atoms, 
            index_nemd=index_nemd, lthermo=lthermo,
            type='npt', output='nemd0.in',
            parameters=parameters, charge=charge, pppm=pppm,
            thermostat=thermostat)
    
    ## NEMD
    for ical in range(nloop):
        num = ical + 1
        parameters['read_data'] = parameters['restart_file']
        parameters['restart_file'] = 'restart%d.nemd'%(num)
        parameters['dump_file'] = 'nemd%d.dump'%(num)
        outfile = "nemd%d.in"%(num)
        write_nemd_input(atoms, 
                index_nemd=index_nemd, lthermo=lthermo,
                type='nemd', output=outfile,
                parameters=parameters, charge=charge, pppm=pppm,
                thermostat=thermostat)

def write_nemd_input(atoms, 
        index_nemd=None, lthermo=None,
        type='npt', output='nemd.in',
        parameters=nemd_default_parameters, 
        charge=False, 
        pppm=1e-5,
        thermostat=None
        ):
    """ Make LAMMPS scripts for NEMD simulations. Heat flows along z-axis.
    atoms : Atoms object
    index_nemd : dictionary of array of integers
        index_nemd[name][2]
    """
    taverage = (parameters['tcold'] + parameters['thot']) * 0.5

    ## get the list of tags which denote the layer index
    tags_list = list(set(atoms.get_tags()))
    tags_list.sort()
    
    ## get species correspongind to the tags
    species_corresp_tags = {}
    for tag in tags_list:
        idx = np.where(atoms.get_tags() == tag)[0][0]
        species_corresp_tags[tag] = atoms[idx].symbol
    
    ## output
    ofs = open(output, 'w')
    ofs.write("clear\n")
    ofs.write("\n")
    ofs.write("# ----- MD parameters -----\n")
    ofs.write("variable tstep   equal %f\n"%(parameters['time_step']))
    ofs.write("variable tdamp   equal %f\n"%(parameters['damping_time']))
    if type == 'npt':
        ofs.write("variable MDTIME  equal %f\n"%(parameters['time_npt']))
    else:
        ofs.write("variable MDTIME  equal %f\n"%(parameters['time_nemd']))
    ofs.write("\n")
    ofs.write("variable NRUN      equal ${MDTIME}/${tstep}\n")
    ofs.write("variable THERSTEP  equal ${NRUN}/%d\n"%(
        int(parameters['ntotal_thermo'])))
    ofs.write("variable DUMPSTEP  equal ${NRUN}/%d\n"%(
        int(parameters['ntotal_dump'])))
    
    ## time for increasing temperature
    if type == 'npt':
        if parameters['time_increase'] is not None:
            ofs.write("variable MDTIME2 equal %f\n"%(
                parameters['time_increase']))
            ofs.write("variable NRUN2     equal ${MDTIME2}/${tstep}\n")
            ofs.write("variable DUMPSTEP2 equal ${NRUN2}/%d\n"%(
                int(parameters['ntotal_dump'])))
    
    ofs.write("variable TEMP      equal %.2f\n"%(taverage))
    ofs.write("\n")
    ofs.write("variable  kb   equal 1.3806488e-23\n")
    ofs.write("variable  ee   equal 1.60217657e-19\n")
    ofs.write("\n")
    ofs.write("# ----- Initialize simulation -----\n")
    ofs.write("units         metal\n")
    if charge == False:
        ofs.write("atom_style    atomic\n")
    else:
        ofs.write("atom_style    charge\n")
    ofs.write("dimension     3\n")
    ofs.write("boundary      p p p\n")
    ofs.write("atom_modify   map array\n")
    if type == 'npt':
        ofs.write("read_data     %s\n"%(parameters['read_data']))
    elif type == 'nemd':
        ofs.write("read_restart  %s\n"%(parameters['read_data']))
    else:
        print(" Error: %s is not defined."%(type))
        exit()
    
    ofs.write("timestep      ${tstep}\n")
    ofs.write("\n")
    
    """ This part should be modified for more complex materials.
    """
    ofs.write("# ----- Define potentials -----\n")
    #nlayers=len(tags_list),
    _write_lj_potential_graphite(ofs, 
            species=species_corresp_tags,
            potfile=parameters['potential_file'],
            charge=charge, pppm=pppm
            )
    
    ofs.write("# ---- Bin ----- \n")
    ofs.write("neighbor      2.0 bin\n")
    ofs.write("neigh_modify  delay 0 every 1 check yes\n")
    ofs.write("\n")
    
    ofs.write("# ----- Groups -----\n")
    if index_nemd is not None:
        _write_groups_with_id(ofs, index_nemd)
    elif lthermo is not None:
        _write_groups_MP_region(ofs, 
                thermo_ratio=lthermo, 
                lsystem=atoms.cell[2,2]
                )
    else:
        print(" Error: index_nemd or lthermo must be given.")
        exit()

    ofs.write("# ----- compute -----\n")
    #ofs.write("compute   displace all displace/atom\n")
    #ofs.write("compute   eng      all pe/atom\n")
    #ofs.write("compute   stress   all stress/atom NULL\n")
    #ofs.write("compute   eatoms   all reduce sum c_eng\n")
    #ofs.write("compute   stressz  all reduce sum c_stress[3]\n")
    ofs.write("compute   1        all temp\n")
    ofs.write("compute   ke       all ke/atom\n")
    ofs.write("\n")
    
    ofs.write("reset_timestep 0\n")
    ofs.write("\n")
    
    if type == 'npt':
        ofs.write("#----- NPT simulation -----\n")
        ofs.write("\n")
        
        ofs.write("#--- 1. minimization\n")
        _write_minimization(ofs, atoms)
        
        ofs.write("#--- 2-1. NPT\n")
        _write_npt_simulation(
                ofs, 
                temperature=taverage, 
                pressure=parameters['pressure'], 
                pdamp=0.05,
                outfile=parameters['dump_file'],
                reset_time=True,
                create_velocity=True,
                drag=2.0
                )
        ofs.write("\n")
        ofs.write("#--- 2-2. NPT\n")
        _write_npt_simulation(
                ofs, 
                temperature=taverage, 
                pressure=parameters['pressure'], 
                pdamp=parameters['pressure_damping'], 
                outfile=parameters['dump_file'],
                reset_time=False,
                create_velocity=False,
                drag=0.0
                )
        
        if parameters['time_increase'] is not None:
            ofs.write("# --- 3. set temperature gradient\n")
            _write_nemd_simulation(ofs, 
                    tinit=taverage,
                    tcold=parameters['tcold'], 
                    thot=parameters['thot'], 
                    outfile='nemd0.dump',
                    names={'dump':'DUMPSTEP2', 'nrun':'NRUN2'},
                    thermostat=thermostat
                    )
        
        ofs.write("#--- prepare fore NEMD simulation\n")
        ofs.write("reset_timestep 0\n")
    
    elif type == 'nemd':
        ofs.write("# ----- NEMD simulation -----\n")
        _write_nemd_simulation(ofs, 
                tcold=parameters['tcold'], 
                thot=parameters['thot'], 
                outfile=parameters['dump_file'],
                thermostat=thermostat
                )
    
    ##
    ofs.write("write_restart  %s\n"%(parameters['restart_file']))
    ofs.write("\n")
    ofs.write("print \"All done!\"\n")
    ofs.close()
    print(" Output", output)

def _write_lj_potential_graphite(ofs, species=None,
        potfile='opt.tersoff', 
        charge=False, pppm=1e-5
        ):
    """ Write LJ potential for graphite
    species : dictionary
        species[tag] : string, element (tag = 1, 2, ...)
    """
    nlayers = len(species)
    
    if charge:
        lj_line = "lj/cut/coul/long"
    else:
        lj_line = "lj/cut"

    ## define potentials
    ofs.write("pair_style hybrid %s 10.0 "%(lj_line))
    for key in species.keys():
        if species[key] == "C":
            ofs.write("tersoff ")
    ofs.write("\n")
    
    ## define interlayer potentials
    ofs.write("\n")
    count = 0
    for i in range(nlayers):
        if species[i+1] != 'C':
            continue
        ofs.write("pair_coeff * * tersoff %d %s"%(count+1, potfile))
        for j in range(nlayers):
            if i == j:
                ofs.write(" C")
            else:
                ofs.write(" NULL")
        count += 1
        ofs.write("\n")
    ofs.write("\n")
    
    ## define intralayer potentials
    for i1 in range(nlayers):
        el1 = species[i1+1]
        for i2 in range(i1, nlayers):
            el2 = species[i2+1]
            if i1 == i2 and el1 == "C":
                continue
            ##
            p1 = _get_lj_parameter_each(el1)
            p2 = _get_lj_parameter_each(el2)
            params = _get_combined_lj_parameter(p1, p2)
            ofs.write("pair_coeff %d %d %s %10.8f %8.5f %5.2f\n"%(
                i1+1, i2+1, lj_line,
                params['epsilon'], 
                params['sigma'], 
                params['cutoff']
                ))
    ofs.write("\n")
    
    if charge:
        ofs.write("kspace_style pppm %.2e\n"%(pppm))
        ofs.write("\n")

def _get_lj_parameter_each(element):
    if element == "Fe":
        return lj_parameters_rappe['Fe']
    elif element == "Cl":
        return lj_parameters_rappe['Cl']
    elif element == "C":
        return lj_parameters_rappe['C']
    else:
        print(" Error", element)
        exit()

def _get_combined_lj_parameter(p1, p2):
    params = {}
    key = 'cutoff'
    params[key] = (p1[key] + p2[key]) * 0.5
    key = 'sigma'
    params[key] = (p1[key] + p2[key]) * 0.5
    key = 'epsilon'
    params[key] = np.sqrt(p1[key] * p2[key])
    return params

#
def _write_groups_with_id(ofs, index, 
        group_list1=['bottom','hot','middle','cold','top'],
        group_list2=['hot','mid1','cold','mid2'],
        ):
    """
    Parameters
    -------------
    index : 
        from index[group_name][0]+1 to index[group_name][1]+1 or
        from index[group_name][ith][0]+1 to index[group_name][ith][1]+1
    """
    for key in index.keys():
        ofs.write('group %8s id '%(key))
        if np.isscalar(index[key][0]):
            ofs.write(' %d:%d'%(index[key][0]+1, index[key][1]+1))
        else:
            for idx_each in index[key]:
                ofs.write(' %d:%d'%(idx_each[0]+1, idx_each[1]+1))
        ofs.write('\n')
    ofs.write("\n")
    
    flag1 = True
    for group in group_list1:
        if group not in index.keys():
            flag1 = False
    
    flag2 = True
    for group in group_list2:
        if group not in index.keys():
            flag2 = False
    
    if flag1 == False and flag2 == False:
        print(" Error in gorup names:", index.keys())
        exit()

    if flag1:
        ofs.write("group const  union bottom top\n")
        ofs.write("group mobile union hot middle cold\n")
        ofs.write("group free   union middle\n")
    else:
        ofs.write("group mobile union mid1 mid2\n")

    ofs.write("\n")

def _write_groups_MP_region(ofs, thermo_ratio=0.1, lsystem=None):
    
    lthermo = thermo_ratio * lsystem
    lmiddle = (0.5 - thermo_ratio) * lsystem
    
    ofs.write("region Rm1 block INF INF INF INF INF %f units box\n"%(
        lmiddle))
    ofs.write("region Rh  block INF INF INF INF %f %f units box\n"%(
        lmiddle, lmiddle+lthermo))
    ofs.write("region Rm2 block INF INF INF INF %f %f units box\n"%(
        lmiddle+lthermo, 2*lmiddle+lthermo))
    ofs.write("region Rc  block INF INF INF INF %f INF units box\n"%(
        2*lmiddle+lthermo))
    ofs.write("\n")
    ofs.write("group hot    region Rh\n")
    ofs.write("group cold   region Rc\n")
    ofs.write("group mid1   region Rm1\n")
    ofs.write("group mid2   region Rm2\n")
    ofs.write("group middle union  mid1 mid2\n")
    ofs.write("\n")

def _write_minimization(ofs, atoms, etol=1e-20, ftol=0.0001, 
        maxiter=50000, maxeval=50000, outfile='mini.dump'):
    """ Minimization is stoped by force tolerance with the default setting.
    """
    size = np.zeros(3)
    for j in range(3):
        size[j] = np.linalg.norm(atoms.cell[j])
    ##ofs.write("fix       1 all box/relax x %.2f y %.2f z %.2f\n"%(
    ##    size[0]*0.1, size[1]*0.1, size[2]*0.1))
    if outfile is not None:
        ofs.write("dump         id1 all custom 100000 mini.dump id type x y z\n")
        ofs.write("dump_modify  id1 sort id\n") 
        ofs.write("dump_modify  id1 format float %13.8f \n") 
    ofs.write("thermo    10000\n")
    ofs.write("minimize  %.5e %.5e %d %d\n"%(etol, ftol, maxiter, maxeval))
    ##ofs.write("unfix     1\n")
    if outfile is not None:
        ofs.write("undump    id1\n")
    ofs.write("\n")
    
def _write_npt_simulation(ofs, outfile=None,
        temperature=300., pressure=1.0, pdamp=0.1, drag=0.0,
        reset_time=True, create_velocity=True):
    """
    parameters
    ---------------
    pdamp : float or list of three values of float
    """
    if reset_time:
        ofs.write("reset_timestep 0\n")
        ofs.write("\n")
    if create_velocity:
        ofs.write("velocity  all create ${TEMP} %d rot yes dist gaussian\n"%(
            np.random.randint(1e7)))
        ofs.write("velocity  all zero linear\n")
        ofs.write("\n")
    
    ofs.write("fix  2 all npt temp %f %f ${tdamp} "%(
        temperature, temperature))
    if type(pdamp) is float:
        ofs.write("aniso %f %f %f drag %f\n"%(pressure, pressure, pdamp, drag))
    elif len(pdamp) == 3:
        #ofs.write("aniso %f %f %f drag %f\n"%(pressure, pressure, pdamp, drag))
        dds = ['x', 'y', 'z']
        ofs.write("&\n")
        for j in range(3):
            ofs.write("     %s %f %f %f "%(
                dds[j], pressure, pressure, pdamp[j]))
            if j != 2:
                ofs.write("&\n")
        ofs.write("\n")
    else:
        print(" Error in pressure damping", pdamp)
        exit()
    ofs.write("\n")
    ofs.write("variable  temp atom c_ke*${ee}/1.5/${kb}\n")
    ofs.write("fix       tempave all ave/atom 1 ${DUMPSTEP} ${DUMPSTEP} v_temp\n")
    ofs.write("\n")
    ofs.write("thermo_style custom step temp press pxx pyy pzz pe ke etotal ")
    ofs.write("evdwl ecoul epair ebond elong etail vol\n")
    ofs.write("thermo       ${THERSTEP}\n")
    if outfile is not None:
        #ofs.write("dump         id2 all custom ${DUMPSTEP} %s "\
        #        "id type x y z f_tempave\n"%(outfile))
        ofs.write("dump         id2 all custom ${DUMPSTEP} %s "\
                "id type x y z f_tempave\n"%(outfile))
        ofs.write("dump_modify  id2 sort id\n") 
        ofs.write("dump_modify  id2 format float %13.8f \n") 
    ofs.write("run       ${NRUN}\n")
    ofs.write("unfix     2\n")
    ofs.write("unfix     tempave\n")
    if outfile is not None:
        ofs.write("undump    id2\n")
    ofs.write("\n")

def _write_nemd_simulation(ofs, thot=310, tcold=290, tinit=None, outfile=None,
        names={'dump':'DUMPSTEP', 'nrun':'NRUN'}, thermostat=None):
    
    ofs.write("compute Thot hot temp\n")
    ofs.write("compute Tcold cold temp\n")
    #ofs.write("fix frozen const setforce 0.0 0.0 0.0\n")
    ofs.write("\n")
    if "nose" in thermostat.lower():
        if tinit is not None:
            ofs.write("fix t1 hot  nvt temp %.2f %.2f ${tdamp}\n"%(tinit, thot))
            ofs.write("fix t2 cold nvt temp %.2f %.2f ${tdamp}\n"%(tinit, tcold))
        else:
            ofs.write("fix t1 hot  nvt temp %.2f %.2f ${tdamp}\n"%(thot, thot))
            ofs.write("fix t2 cold nvt temp %.2f %.2f ${tdamp}\n"%(tcold, tcold))
        ofs.write("fix m1 free nve\n")
    elif 'lan' in thermostat.lower():
        if tinit is not None:
            ofs.write("fix t1 hot  langevin %.2f %.2f ${tdamp} %d tally yes \n"%(
                tinit, thot, np.random.randint(1e7)))
            ofs.write("fix t2 cold langevin %.2f %.2f ${tdamp} %d tally yes \n"%(
                tinit, tcold, np.random.randint(1e7)))
        else:
            ofs.write("fix t1 hot  langevin %.2f %.2f ${tdamp} %d tally yes \n"%(
                thot, thot, np.random.randint(1e7)))
            ofs.write("fix t2 cold langevin %.2f %.2f ${tdamp} %d tally yes \n"%(
                tcold, tcold, np.random.randint(1e7)))
        ##
        #ofs.write("fix m1 mobile nve\n")
        #ofs.write("fix m1 middle nve\n")
        ofs.write("fix a1 all nve\n")
    
    ofs.write("fix_modify t1 temp Thot\n")
    ofs.write("fix_modify t2 temp Tcold\n")
    ofs.write("\n")
    #ofs.write("fix 2 mobile ave/spatial 1 ${%s} ${%s} v_temp &\n"%(
    #    names['dump'], names['dump']))
    #ofs.write("      file tmp.profile units reduced\n")
    
    #ofs.write("fix mom all momentum 1 linear 1 1 1\n")
    ofs.write("variable     temp atom c_ke*${ee}/1.5/${kb}\n")
    ofs.write("variable     diff equal f_t2-f_t1\n")
    ofs.write("fix          tempave all ave/atom "\
            "1 ${%s} ${%s} v_temp\n"%(names['dump'], names['dump']))
    #ofs.write("fix          tempave all ave/time "\
    #        "1 ${%s} ${%s} c_1 file temp_ave.txt\n"%(names['dump'], names['dump']))
    ofs.write("\n")
    ofs.write("thermo_style custom step pe etotal temp c_Thot c_Tcold f_t1 f_t2 v_diff\n")
    ofs.write("thermo       ${THERSTEP}\n")
    if outfile is not None:
        #ofs.write("dump         nemd all custom ${DUMPSTEP} %s "\
        #        "id type x y z vx vy vz f_tempave\n"%(outfile))
        ofs.write("dump         nemd all custom ${%s} %s "\
                "id type x y z f_tempave\n"%(names['dump'], outfile))
        ofs.write("dump_modify  nemd sort id\n") 
        ofs.write("dump_modify  nemd format float %13.8f \n") 
    
    ofs.write("run          ${%s}\n"%(names['nrun']))
    if outfile is not None:
        ofs.write("undump     nemd\n")
    ofs.write("\n")
    
