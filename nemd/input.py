# -*- coding: utf-8 -*-
import numpy as np

def write_nemd_input(atoms, index_nemd, type='npt', 
        read_data='data.lammps', output='nemd.in', 
        dumpfile=None, restartfile='restart.nemd', 
        timestep=0.0005, damping_time=0.05, md_time=None, 
        ntotal_thermo=100, ntotal_dump=10,
        taverage=300., tcold=290., thot=310.,
        pressure=1.0, 
        lj_params={'cutoff':12.00, 'epsilon':0.00239, 'sigma':3.4}):
    """

    """
    ## get the list of tags which denote the layer index
    tags_list = list(set(atoms.get_tags()))
    tags_list.sort()
    
    ## output
    ofs = open(output, 'w')
    ofs.write("clear\n")
    ofs.write("\n")
    ofs.write("# ----- MD parameters -----\n")
    ofs.write("variable tstep   equal %f\n"%(timestep))
    ofs.write("variable tdamp   equal %f\n"%(damping_time))
    ofs.write("variable MDTIME  equal %f\n"%(md_time))
    ofs.write("\n")
    ofs.write("variable NRUN      equal ${MDTIME}/${tstep}\n")
    ofs.write("variable THERSTEP  equal ${NRUN}/%d\n"%(int(ntotal_thermo)))
    ofs.write("variable DUMPSTEP  equal ${NRUN}/%d\n"%(int(ntotal_dump)))
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
    ofs.write("read_data     %s\n"%(read_data))
    ofs.write("timestep      ${tstep}\n")
    ofs.write("\n")
    
    """ This part should be modified for more complex materials.
    """
    ofs.write("# ----- Define LJ potential -----\n")
    _write_lj_potential_graphite(ofs, lj_params, nlayers=len(tags_list))
    
    
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
                pressure=pressure, 
                pdamp=timestep*1000.,
                outfile=dumpfile
                )
        
        ofs.write("#--- prepare fore NEMD simulation\n")
        ofs.write("reset_timestep 0\n")
    
    elif type == 'nemd':
        ofs.write("# ----- NEMD simulation -----\n")
        _write_nemd_simulation(
                ofs, tcold=tcold, thot=thot, outfile=dumpfile)
    
    ##
    ofs.write("write_restart  %s\n"%(restartfile))
    ofs.write("\n")
    ofs.write("print \"All done!\"\n")

def _write_lj_potential_graphite(ofs, lj_params, nlayers=None):
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
        ofs.write("pair_coeff * * tersoff %d opt.tersoff"%(i+1))
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
    ofs.write("group mobile union hot middle cold full\n")
    ofs.write("group free   union middle full\n")
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
        ofs.write("dump      id1 all custom 100000 mini.dump id type x y z\n")
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
        ofs.write("dump      id2 all custom ${DUMPSTEP} %s "%(outfile))
        ofs.write("id type x y z vx vy vz f_tempave\n")
    ofs.write("run       ${NRUN}\n")
    ofs.write("unfix     2\n")
    ofs.write("unfix     tempave\n")
    if outfile is not None:
        ofs.write("undump    id2\n")
    ofs.write("\n")

def _write_nemd_simulation(ofs, thot=310, tcold=290, outfile=None):
    
    ofs.write("fix frozen const setforce 0.0 0.0 0.0\n")
    ofs.write("fix t1 hot  nvt temp %.2f %.2f ${tdamp}\n"%(thot, thot))
    ofs.write("fix t2 cold nvt temp %.2f %.2f ${tdamp}\n"%(tcold, tcold))
    ofs.write("fix m1 free nve\n")
    ofs.write("\n")
    ofs.write("variable temp atom c_ke*${ee}/1.5/${kb}\n")
    ofs.write("variable diff equal f_t2-f_t1\n")
    ofs.write("\n")
    ofs.write("variable   temp atom c_ke*${ee}/1.5/${kb}\n")
    ofs.write("fix        tempave all ave/atom 1 ${THERSTEP} ${THERSTEP} v_temp\n")
    
    ofs.write("#--- OUTPUT\n")
    ofs.write("thermo_style custom step pe temp f_t1 f_t2 v_diff\n")
    ofs.write("thermo       ${THERSTEP}\n")
    if outfile is not None:
        ofs.write("dump         nemd all custom ${DUMPSTEP} %s "%(outfile))
        ofs.write("id type x y z f_tempave\n")
    ofs.write("run          ${NRUN}\n")
    if outfile is not None:
        ofs.write("undump nemd\n")


