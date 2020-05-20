# -*- coding: utf-8 -*-
import numpy as np
from nemd import (
        nemd_default_parameters, lj_parameters_C)

def write_nemd_inputs(atoms, index_nemd, datafile='data.lammps',
        thot=310, tcold=290, 
        time_npt=10., time_increase=100., time_nemd=100., 
        nloop=2, 
        output_minimization=True,
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
    write_nemd_input(atoms, index_nemd, type='npt', output='nemd0.in',
            parameters=parameters,
            lj_parameters=lj_parameters_C)

    ## NEMD
    for ical in range(nloop):
        num = ical + 1
        parameters['read_data'] = parameters['restart_file']
        parameters['restart_file'] = 'restart%d.nemd'%(num)
        parameters['dump_file'] = 'nemd%d.dump'%(num)
        outfile = "nemd%d.in"%(num)
        write_nemd_input(atoms, index_nemd, type='nemd', output=outfile,
                parameters=parameters,
                lj_parameters=lj_parameters_C)

def write_nemd_input(atoms, index_nemd, type='npt', output='nemd.in',
        parameters=nemd_default_parameters,
        lj_parameters=lj_parameters_C):
    
    """ Make LAMMPS scripts for NEMD simulations
    atoms : Atoms object
    index_nemd : dictionary of array of integers
        index_nemd[name][2]
    """
    taverage = (parameters['tcold'] + parameters['thot']) * 0.5

    ## get the list of tags which denote the layer index
    tags_list = list(set(atoms.get_tags()))
    tags_list.sort()
    
    ## output
    ofs = open(output, 'w')
    ofs.write("clear\n")
    ofs.write("\n")
    ofs.write("# ----- MD parameters -----\n")
    ofs.write("variable tstep   equal %f\n"%(parameters['time_step']))
    ofs.write("variable tdamp   equal %f\n"%(parameters['damping_time']))
    if type == 'npt':
        md_time = parameters['time_npt']
        md_time2 = parameters['time_increase']
    else:
        md_time = parameters['time_nemd']
    ofs.write("variable MDTIME  equal %f\n"%(md_time))
    ofs.write("\n")
    ofs.write("variable NRUN      equal ${MDTIME}/${tstep}\n")
    ofs.write("variable THERSTEP  equal ${NRUN}/%d\n"%(
        int(parameters['ntotal_thermo'])))
    ofs.write("variable DUMPSTEP  equal ${NRUN}/%d\n"%(
        int(parameters['ntotal_dump'])))
    
    ## time for increasing temperature
    if type == 'npt':
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
    ofs.write("atom_style    atomic\n")
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
    ofs.write("# ----- Define LJ potential -----\n")
    _write_lj_potential_graphite(ofs, lj_parameters, nlayers=len(tags_list),
            potfile=parameters['potential_file'])
    
    ofs.write("# ---- Bin ----- \n")
    ofs.write("neighbor      0.3 bin\n")
    ofs.write("neigh_modify  delay 0\n")
    ofs.write("\n")
    
    ofs.write("# ----- Group -----\n")
    _write_groups(ofs, index_nemd)
    
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
        
        ofs.write("#--- 2. NPT\n")
        _write_npt_simulation(
                ofs, 
                temperature=taverage, 
                pressure=parameters['pressure'], 
                pdamp=parameters['pressure_damping'], 
                outfile=parameters['dump_file']
                )
        
        ofs.write("# --- 3. set temperature gradient\n")
        _write_nemd_simulation(ofs, 
                tinit=taverage,
                tcold=parameters['tcold'], 
                thot=parameters['thot'], 
                outfile='nemd0.dump',
                names={'dump':'DUMPSTEP2', 'nrun':'NRUN2'})
        
        ofs.write("#--- prepare fore NEMD simulation\n")
        ofs.write("reset_timestep 0\n")
    
    elif type == 'nemd':
        ofs.write("# ----- NEMD simulation -----\n")
        _write_nemd_simulation(ofs, 
                tcold=parameters['tcold'], 
                thot=parameters['thot'], 
                outfile=parameters['dump_file'])
    
    ##
    ofs.write("write_restart  %s\n"%(parameters['restart_file']))
    ofs.write("\n")
    ofs.write("print \"All done!\"\n")
    ofs.close()
    print(" Output", output)

def _write_lj_potential_graphite(ofs, lj_params, nlayers=None,
        potfile='opt.tersoff'):
    """ Write LJ potential for graphite
    lj_params : dictionary 
        lj_params[key] : float
        key = 'cutoff', 'epsilon', and 'sigma'
    """
    
    ## output LJ parameters
    for key in lj_params.keys():
        ofs.write("variable %8s equal %.7f\n"%(key, lj_params[key]))
    ofs.write("\n")
    
    ## define potentials
    ofs.write("pair_style hybrid lj/cut ${cutoff}")
    for i in range(nlayers):
        ofs.write(" tersoff")
     
    ## define interlayer potentials
    ofs.write("\n")
    for i in range(nlayers):
        ofs.write("pair_coeff * * tersoff %d %s"%(i+1, potfile))
        for j in range(nlayers):
            if i == j:
                ofs.write(" C")
            else:
                ofs.write(" NULL")
        ofs.write("\n")
    
    ## define intralayer potentials
    for i1 in range(nlayers-1):
        for i2 in range(i1+1, nlayers):
            ofs.write("pair_coeff %d %d lj/cut ${epsilon} ${sigma}\n"%(
                i1+1, i2+1))
    ofs.write("\n")
    
#
def _write_groups(
        ofs, index, 
        group_list=['bottom','hot','middle','cold','top']):
    for key in index.keys():
        ofs.write('group %8s id %d:%d\n'%(
            key, index[key][0]+1, index[key][1]+1))
    ofs.write("\n")
    for group in group_list:
        if group not in index.keys():
            print(" Error: %s region does not exist."%(group))
            exit()
    ofs.write("group const  union bottom top\n")
    ofs.write("group mobile union hot middle cold\n")
    ofs.write("group free   union middle\n")
    ofs.write("\n")

def _write_minimization(ofs, atoms, etol=1e-20, ftol=0.0001, 
        maxiter=50000, maxeval=50000, outfile='mini.dump'):
    """ Minimization is stoped by force tolerance with the default setting.
    """
    size = np.zeros(3)
    for j in range(3):
        size[j] = np.linalg.norm(atoms.cell[j])
    ofs.write("fix       1 all box/relax x %.2f y %.2f z %.2f\n"%(
        size[0]*0.1, size[1]*0.1, size[2]*0.1))
    if outfile is not None:
        ofs.write("dump         id1 all custom 100000 mini.dump id type x y z\n")
        ofs.write("dump_modify  id1 sort id\n") 
        ofs.write("dump_modify  id1 format float %13.8f \n") 
    ofs.write("thermo    10000\n")
    ofs.write("minimize  %.5e %.5e %d %d\n"%(etol, ftol, maxiter, maxeval))
    ofs.write("unfix     1\n")
    if outfile is not None:
        ofs.write("undump    id1\n")
    ofs.write("\n")
    
def _write_npt_simulation(ofs, outfile=None,
        temperature=300., pressure=1.0, pdamp=0.1):
    ofs.write("reset_timestep 0\n")
    ofs.write("\n")
    ofs.write("velocity  mobile create ${TEMP} %d rot yes dist gaussian\n"%(
        np.random.randint(1e7)))
    #ofs.write("fix       2 mobile nvt temp 300. 300. ${tdamp}\n")
    ofs.write("fix       2 all npt temp %f %f ${tdamp} iso %f %f %f\n"%(
        temperature, temperature, pressure, pressure, pdamp))
    ofs.write("\n")
    ofs.write("variable  temp atom c_ke*${ee}/1.5/${kb}\n")
    ofs.write("fix       tempave all ave/atom 1 ${DUMPSTEP} ${DUMPSTEP} v_temp\n")
    ofs.write("\n")
    ofs.write("thermo    ${THERSTEP}\n")
    if outfile is not None:
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
        names={'dump':'DUMPSTEP', 'nrun':'NRUN'}):
    
    ofs.write("fix frozen const setforce 0.0 0.0 0.0\n")
    if tinit is not None:
        ofs.write("fix t1 hot  nvt temp %.2f %.2f ${tdamp}\n"%(tinit, thot))
        ofs.write("fix t2 cold nvt temp %.2f %.2f ${tdamp}\n"%(tinit, tcold))
    else:
        ofs.write("fix t1 hot  nvt temp %.2f %.2f ${tdamp}\n"%(thot, thot))
        ofs.write("fix t2 cold nvt temp %.2f %.2f ${tdamp}\n"%(tcold, tcold))
    ofs.write("fix m1 free nve\n")
    ofs.write("\n")
    ofs.write("variable     temp atom c_ke*${ee}/1.5/${kb}\n")
    ofs.write("variable     diff equal f_t2-f_t1\n")
    ofs.write("fix          tempave all ave/atom "\
            "1 ${%s} ${%s} v_temp\n"%(names['dump'], names['dump']))
    ofs.write("\n")
    ofs.write("thermo_style custom step pe temp f_t1 f_t2 v_diff\n")
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
    
