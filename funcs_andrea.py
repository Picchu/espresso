from model_def import *
from subtypes import *
from PPlibraries import *
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.orm.attributes import flag_modified
from sqlalchemy import create_engine
from wflib import *
from ase import io
from collections import OrderedDict
from atoms_info import *
import numpy as np
from saga_funcs import *
from espresso import *
import re, time, os
from dbinit import *

# warnings dictionary used to store information about possible errors from QE calc
WARNINGS = {'convergence NOT achieved after': 'Number of scf steps exceeded',
    'not orthogonal operation': 'Not orthogonal operation in checkallsym',
    'Not enough space allocated for radial FFT': 'Not enough space allocated for radial FFT: \
    increase cell_factor', 'Not a group': 'Error not a group in routine multable',
    'Aborted with signal SIGUSR2': 'Calculation interrupted',
    'lone vector': 'Error lone vector in sym_rho_init_shell',
    'too many bands are not converged': 'Too many not converged bands in c_bands',
    'The maximum number of steps has been reached': 'Number of BFGS steps exceeded',
    'invalid values for nk1, nk2, nk3': 'Invalid k-points grid',
    'MKL ERROR': 'MKL ERROR', 'problems computing cholesky': 'Error in Cholesky decomposition'}

def get_warnings(filename):
    """
    This function parses the output file of a QE calc, and stores the
    number of occurrences for a given error message.
    :param: the file to parse (full path)
    :return: dict where the keys are warning messages, and the
    values the number of occurrences in the file
    """
    fn = open(filename, 'r')
    output = fn.read()
    fn.close()

    return {value: len(re.findall(key, output)) for \
    key, value in WARNINGS.iteritems()}

def check_number(s):
    """
    Function used to interpret correctly values for different
    fields in a Quantum Espresso input file (i.e., if the user
    specifies the string '1.0d-4', the quotes must be removed
    when we write it to the actual input file)
    :param: a string
    :return: formatted string to be used in QE input file
    """
    try:
        float(s.replace('d','e'))
        return s
    except ValueError:
        if s == '.true.':
            return s
        else:
            return "'%s'" % s

def upload_PPs_to_db(*PPs):
    """
    Function that loads PseudoPotentials to the database. The input
    is a list of PPs in the format specified by file PPlibraries.py
    :param: a list of PPs in dict format
    :return: a dict of the PPs object added, {element: PP DB object}
    """
    session = Session()
    new_count = 0

    pseudo_list = {}

    for PP in PPs:
        name = PP['filename']
        label = 'PP'
        content = PP

        # upload to database only if we have a new PP object
        # avoid uploading duplicated to database
        # a unique PP is identified by its name and content
        pseudo, new_pseudo = get_or_create(session, PseudoPot, name=name, content=content)

        if new_pseudo == True:
            new_count += 1

        session.add(pseudo)

        pseudo_list[content['element']] = pseudo

    session.commit()

    session.close()

    print "Uploaded %d new PPs to the database." % new_count

    return pseudo_list

def import_structures_to_db(*cif_files):
    """
    This function upload Struc objects to the database,
    from a list of .cif files.
    :param: a list of .cif files
    :return: a dict of the Struc objects added, {struc_name: Struc DB object}
    """

    session = Session()
    new_count = 0

    struc_list = {}

    # use python ASE to extract information from cif file
    for cif_file in cif_files:
        # Atoms object in ASE
        at = io.read(cif_file, format='cif')
        name = cif_file.replace('.cif', '')
        chem_list = at.get_chemical_symbols()
        f = open(cif_file,'r')
        cif_content = f.read()
        f.close()
        content = {}
        composition = {}
        # the full content of the cif file used to generate the structure
        content['cif'] = cif_content
        # if periodic boundary conditions are applied (1=True, 0=False)
        content['pbc'] = [1 if x else 0 for x in at.pbc]
        # cell vectors in angstrom
        content['cell'] = at.cell.tolist()
        # absolute atomic positions in angstrom
        content['positions'] = at.get_positions().tolist()
        # relative atomic positions
        content['scaled_positions'] = at.get_scaled_positions().tolist()
        # species of each atom in the cell
        content['species'] = chem_list
        # nice compact formula representation
        content['formula'] = at.get_chemical_formula()
        # number of atoms in the cell
        content['nat'] = at.get_number_of_atoms()
        # volume of the cell
        content['volume'] = at.get_volume()
        # composition in dict format, i.e. for LiO2
        # {'Li': 1, 'O':2}
        for i in set(chem_list):
            composition[i] = chem_list.count(i)
        content['composition'] = composition

        # avoid duplicate structures in the database
        struc, new_struc = get_or_create(session, Struc, name=name, content=content)

        if new_struc == True:
            new_count += 1

        session.add(struc)

        struc_list[name] = struc

    session.commit()

    session.close()

    print "Uploaded %d new structures to the database." % new_count

    return struc_list


def import_qe_input_to_db(name, control = {}, system = {}, electrons = {},
    ions = {}, cell = {}, kspacing=0.07, kshift=[0,0,0]):
    """
    Function which generates a QEInput object starting from a set of dicts
    containing information about a QE calc (see test.py for examples).
    :param: dicts of parameters for a given field of a QE calc
    :return: dict for the QEInput added to db, {qeinput_name: QEInput DB object}
    """

    session = Session()
    new_count = 0

    content = {}
    # &CONTROL field in QE calc
    if control:
        content['control'] = control
    else:
        raise Exception("Control field not specified")
    # &SYSTEM field in QE calc
    if system:
        content['system'] = system
    else:
        raise Exception("System field not specified")
    # &ELECTRONS field in QE calc
    if electrons:
        content['electrons'] = electrons
    else:
        raise Exception("Electrons field not specified")
    # &IONS field in QE calc (may not be present)
    if ions:
        content['ions'] = ions
    # &CELL field in QE calc (may not be present)
    if cell:
        content['cell'] = cell
    # KPOINTS field in QE calc
    content['kspacing'] = kspacing
    content['kshift'] = kshift

    # avoid duplicates
    qeinp, new_qeinp = get_or_create(session, QEInput, name=name, content=content)

    if new_qeinp == True:
        new_count += 1

    session.add(qeinp)

    session.commit()

    session.close()

    print "Uploaded %d new QE input files to the database." % new_count

    return {name: qeinp}

def close_js(username, host, js_type):
    """
    Function to be called at the end of a workflow, in order
    to explicitly close a job service instance
    :param: username and hostname (i.e. 'gra2pal', 'pal-hpc-login')
    together with the job service type ('ssh', 'sge', 'local')
    """
    js = get_js(username, host, js_type)

    if (js is not None) and (js.valid == True):
        js.close()

    return


@threadit
@workit_factory(Combine)
def create_PP_set(self, code_name=None, **db_PPs):
    """
    This function takes as input a dict of PseudoPot objects, and creates
    a PPSubset object which contains the PseudoPots. Extra links are created
    so that the user can access to the PseudoPots directly from the set. (AG)
    :param: dict of PP DB objects, the code_name allows the user to specify
    an existing code in the database to be linked to this calc by the wrapper.
    """

    session = Session()

    uuid_list = []
    chemsys = []

    for el, PP in db_PPs.iteritems():
        # keep track of the UUIDs of the input PPs
        uuid_list.append(str(PP.uuid))
        # sanity check to ensure we have a PP consistent with the
        # elements we want
        if el != PP.content['element']:
            raise Exception("Discrepancy between element and PP")
        chemsys.append(el)

    # ensure we have at most one PP for a given element
    if len(chemsys) != len(set(chemsys)):
        raise Exception("Multiple PPs found for the same element")

    name = '_'.join(sorted(chemsys)) + '.pp_set'
    content = {'species': ','.join(sorted(chemsys)), 'PP_uuids': ','.join(sorted(uuid_list))}

    # create new PPSubset; since this is a work function, the existence
    # of a similar PPSubset will be taken care of by the wrapper itself
    # no need to use get or create here
    PP_set = PPSubset(name=name, content=content)

    # create extra links so that we can access the PPs directly from
    # the PP set without having to go through the CombineCalc.
    # This links are not set automatically by the wrapper, hence they
    # must be explicitly set by the user
    for el, PP in db_PPs.iteritems():
        PP_set.roots.append(PP)

    session.add(PP_set)

    session.commit()

    session.close()

    return {name: PP_set}


@threadit
@workit_factory(QEplugin)
def create_QEInput_file(self, code_name=None, **inp):
    """
    This function takes as input a dict of DB objects (one Struc, one QEInput, one PP set)
    and creates a QE Input file object which contains the full qe file.
    :param: dict with Struc obejct, PPSubset object, QEInput object; ,
    the code_name allows the user to specify an existing code in the
    database to be linked to this calc by the wrapper.
    :return: dict with the FileData instance created
    """
    uuid_list = []

    # currently we distinguish the inputs by type,
    # hence no need to use specific labels to retrieve
    # them. We will need to add checks to ensure these
    # are also the only inputs provided by the user
    for key, value in inp.iteritems():
        uuid_list.append(str(value.uuid))
        if isinstance(value, Struc):
            struc = value
        elif isinstance(value, QEInput):
            qeinp = value
        elif isinstance(value, PPSubset):
            ppsub = value

    # use ordered dict to ensure consistency of format between
    # the different input files created
    # info taken from QEInput
    control_od = OrderedDict(sorted(qeinp.content['control'].items()))
    system_od = OrderedDict(sorted(qeinp.content['system'].items()))
    electrons_od = OrderedDict(sorted(qeinp.content['electrons'].items()))
    ions_od = OrderedDict(sorted(qeinp.content['ions'].items()))
    cell_od = OrderedDict(sorted(qeinp.content['cell'].items()))
    kspacing = qeinp.content['kspacing']
    kshift = qeinp.content['kshift']
    title = control_od['title']
    work_dir = control_od['outdir']

    # info taken from Struc
    system_od['nat'] = str(struc.content['nat'])
    system_od['ntyp'] = str(len(struc.content['composition'].keys()))

    # used to have different labels for species of same type,
    # i.e. Li0, Li1, Li2.... useful when we want to distinguish
    # atoms of the same species, i.e. when doing U calculations
    # and we need to perturb only a specific site
    species_map = {}

    # spin pol calculation if at least one atom is in the list of magnetic ions
    if any(x in magnetic_ions for x in struc.content['composition'].keys()):
        system_od['nspin'] = '2'
    j = 0
    # loop over the species in the structure
    for el in sorted(struc.content['composition'].keys()):
        # mapping, i.e. Li -> Li0, Cl -> Cl1
        # need to change this when we differentiate atoms probably
        species_map[el] = '%s%d' % (el,j)
        # set starting magnetisation: currently we initialize spin polarised
        # ions to +1 and non-spin polarised ions to 0
        if el in magnetic_ions:
            pol = '1.0'
        else:
            pol = '0.0'
        system_od['starting_magnetization(' + str(j+1) + ')'] = pol
        j += 1

    system_od = OrderedDict(sorted(system_od.items()))
    species_map = OrderedDict(sorted(species_map.items()))

    ss = ''
    ss += "&CONTROL\n"
    for key, value in control_od.iteritems():
        ss += "%+25s = %s \n" % (key, check_number(value))
    ss += '/\n'
    ss += "&SYSTEM\n"
    for key, value in system_od.iteritems():
        ss += "%+25s = %s \n" % (key, check_number(value))
    ss += '/\n'
    ss += "&ELECTRONS\n"
    for key, value in electrons_od.iteritems():
        ss += "%+25s = %s \n" % (key, check_number(value))
    ss += '/\n'
    # only if dict is present
    if ions_od:
        ss += "&IONS\n"
        for key, value in ions_od.iteritems():
            ss += "%+25s = %s \n" % (key, check_number(value))
        ss += '/\n'
    # only if dict is present
    if cell_od:
        ss += "&CELL\n"
        for key, value in cell_od.iteritems():
            ss += "%+25s = %s \n" % (key, check_number(value))
        ss += '/\n'
    ss += "ATOMIC_SPECIES\n"
    for key, value in species_map.iteritems():
        # assume a set contains only one PP for a given element
        # this is ensured by the create_PP_set func
        # retrieve PP filename from the PP object
        pp = [x for x in ppsub.roots if x.content['element']==key][0]
        ss += "%+10s    %s    %s\n" % (value, atomic_mass[key], pp.content['filename'])
    # default: use relative positions and explicit cell vectors
    ss += "ATOMIC_POSITIONS crystal\n"
    for i in range(len(struc.content['species'])):
        # e.g.: Li0 x_coord y_coord z_coord
        ss += ' %s %.6f %.6f %.6f\n' % (species_map[struc.content['species'][i]],
            struc.content['scaled_positions'][i][0], struc.content['scaled_positions'][i][1],
            struc.content['scaled_positions'][i][2])

    # extract cell vectors from cell field in Struc object
    a_np = np.array(struc.content['cell'][0], dtype=float)
    b_np = np.array(struc.content['cell'][1], dtype=float)
    c_np = np.array(struc.content['cell'][2], dtype=float)

    # take magnitude of each cell vector
    a_mod = np.sqrt(a_np.dot(a_np))
    b_mod = np.sqrt(b_np.dot(b_np))
    c_mod = np.sqrt(c_np.dot(c_np))

    # determine k-point grid from spacing in reciprocal space
    KP_x = int(1./(a_mod*kspacing)+0.5)
    KP_y = int(1./(b_mod*kspacing)+0.5)
    KP_z = int(1./(c_mod*kspacing)+0.5)

    # ensures minimal dimensions for very large cells
    KP_x = KP_x if KP_x > 0 else 1
    KP_y = KP_y if KP_y > 0 else 1
    KP_z = KP_z if KP_z > 0 else 1

    ss += "K_POINTS automatic\n"
    ss += " %d %d %d    %d %d %d\n" % (KP_x, KP_y, KP_z, kshift[0], kshift[1], kshift[2])

    ss += "CELL_PARAMETERS angstrom\n"
    ss += "%.6f %.6f %.6f\n" % (a_np[0], a_np[1], a_np[2])
    ss += "%.6f %.6f %.6f\n" % (b_np[0], b_np[1], b_np[2])
    ss += "%.6f %.6f %.6f\n" % (c_np[0], c_np[1], c_np[2])

    content = {}
    content['filepath'] = work_dir + "/%s.in" % title
    content['filecontent'] = ss
    content['uuids'] =  ','.join(sorted(uuid_list))

    fileinp = FileData(name=title+'_input', content=content)

    return {title: fileinp}

@threadit
@workit_factory(QEplugin)
def run_QE_calc(self, code_name=None, **inp):
    """
    This function takes as input a QE input file DB object: it runs the calculation
    (locally or remotely), then returns a FileData object containing information
    about the output file generated by the calculation. The pasring of the file
    will be done by a separate routine in order to allow for different parsers to
    be used. (AG)
    :param: dict with FileData object; the code_name allows the user to specify
    an existing code in the database to be linked to this calc by the wrapper.
    :return: dict with the FileData instance created
    """
    uuid_list = []

    # TO-DO: need to check this is the only input
    for key, value in inp.iteritems():
        uuid_list.append(str(value.uuid))
        if isinstance(value, FileData):
            fileinp = value

    # full path of the raw input file (local or remote)
    filepath = fileinp.content['filepath']
    filecontent = fileinp.content['filecontent']

    # create input file for QE calculation
    # we will need to use saga for remote files
    fn = open(filepath, 'w')
    fn.write(filecontent)
    fn.close()

    # extract seed and working directory from filepath
    seed = filepath.replace('.in', '')
    cwd = '/'.join(seed.split('/')[:-1])

    # saga job states to check
    completed_states = [saga.job.DONE, saga.job.FAILED, saga.job.CANCELED]

    # probably better to move this information to the code instance or to
    # another data object which specifies executable path and other options
    # such potentially number of pools to use
    executable = 'mpirun -np 4 pw.x -i'

    # create saga job
    container, js = run_job_saga(hal_username,'', 'local', executable,
            inputs = [filepath], name = None, project = None, queue = None,
            wall_time_limit = None, total_cpu_count = None, spmd_variation = None,
            workdir = [cwd], outputs = [seed+'.out'], errors = [seed+'.err'])

    # run the saga job
    container.run()

    # get jobs id
    job_id_list = [job.id for job in container.jobs]

    # store jobs id into content of the calc object in database
    s = Session()
    self_copy = s.merge(self)
    # in case the content field is empty
    if not self_copy.content:
        self_copy.content = {}
    self_copy.content['jobs_id'] = job_id_list
    flag_modified(self_copy, "content")
    s.commit()

    # test: periodically check number of jobs running
    while any(x not in completed_states for x in container.get_states()) == True:
        running_jobs = [y for y in container.jobs if y.get_state() not in completed_states]
        print "Number of running jobs: %d" % len(running_jobs)

        time.sleep(11)

    res_list = []

    # when jobs are finished, get states of the jobs
    # this is important in order to determine if the
    # output will be in a Final or Partial state
    for state in container.get_states():
        if state == saga.job.DONE:
            res_list.append("Finished")
        elif state == saga.job.FAILED:
            res_list.append("Failed")
        elif state == saga.job.CANCELED:
            res_list.append("Canceled")

    print res_list

    content = {}
    content['filepath'] = seed+'.out'
    # currently assuming one job per saga task
    # convention: job finished -> data is final
    # otherwise data is partial
    # validity will be checked automatically by methods invoked by wrapper
    if res_list[0] == "Finished":
        state = "Final"
    else:
        state = "Partial"

    title = seed.split('/')[-1]
    # save FileData object to database; no need to check for duplicates since
    # the wrapper will take care of fast-forwarding if calculation with the same
    # inputs is found in the database
    fileout = FileData(name=title+'_output', state=state, content=content)

    return {title: fileout}

@threadit
@workit_factory(QEplugin)
def parsing_QE_output(self, code_name=None, **inp):
    """
    This function takes as input a FileData object with the location of the output
    file generated by QE. It then parses the file and store the relevant info into a
    Data object. The use of the code_name allows to specify the parser to be used.
    Currently only parsing using espresso python module is allowed. (AG)
    :param: FileData object created by QE calc
    :return: dict with QEOutput object containing information about parsed QE output.
    TO-DO: implement different parsing options
    """

    uuid_list = []

    # TO-DO: check these are the only inputs
    for key, value in inp.iteritems():
        uuid_list.append(str(value.uuid))
        if isinstance(value, Observed):
            out = value
        if isinstance(value, Parser):
            parser = value

    # currently this is the only parsing option supported
    if parser.name != 'Espresso wrapper by Zhongnan Xu':
        raise Exception("{} not currently supported".format(parser))
    else:
        print "Analyzing QE output with {}".format(parser)

    title = out.content['filepath'].split('/')[-1]

    seed = out.content['filepath'].replace('.out','')

    content = {}

    try:
        # using espresso module with some changes done by me
        # in order to initialise Espresso instance from filename
        e = Espresso(filename=seed)
        # True if calculation is finished
        content['calc_finished'] = e.calc_finished
        # True if calculation convergenced
        content['converged'] = e.converged
        # check for electronic convergence
        content['electronic_converged'] = e.electronic_converged
        content['total_energy'] = e.energy_free
        content['energy_atom'] = e.energy_free/e.int_params['nat']
        if e.int_params['nspin'] == 2:
            spin_pol = True
        else:
            spin_pol = False
        # save some of the input parameters for easier queries
        content['spin_pol'] = spin_pol
        content['calc_type'] = e.string_params['calculation']
        content['title'] = e.string_params['title']
        content['ecutwfc'] = e.real_params['ecutwfc']
        content['ecutrho'] = e.real_params['ecutrho']
        content['formula'] = e.atoms.get_chemical_formula()
        content['nat'] = e.int_params['nat']
        content['kpts'] = e.input_params['kpts'].tolist()
        content['offset'] = e.input_params['offset']
        # initial cell vectors
        content['cell_in'] = e.all_cells[0].tolist()
        # initial atomic positions
        content['pos_in'] = e.all_pos[0].tolist()
        # final cell vectors
        if len(e.all_cells) > 1:
            content['cell_fin'] = np.vstack((e.all_cells[-1][0],e.all_cells[-1][1],e.all_cells[-1][2])).tolist()
        else:
            content['cell_fin'] = e.all_cells[0].tolist()
        # final atomic positions
        content['pos_fin'] = e.all_pos[-1].tolist()
        chem_list = e.atoms.get_chemical_symbols()
        composition = {}
        for i in set(chem_list):
            composition[i] = chem_list.count(i)
        content['composition'] = composition
        # keep track of UUIDs of inputs
        content['uuids'] =  ','.join(sorted(uuid_list))
        # output file full path
        content['filepath'] = seed + '.out'
        # warnings from output file
        content['warnings'] = get_warnings(seed + '.out')
        # U values in DFT+U calculation
        content['Hubbard_U'] = e.list_params['Hubbard_U']
        # J values in DFT+U+J calculation
        content['Hubbard_J'] = e.list_params['Hubbard_J']
        # alpha values in linear response DFT+U calculation
        content['Hubbard_alpha'] = e.list_params['Hubbard_alpha']
        # applied electric field
        content['efield_cart'] = e.list_params['efield_cart']
        # starting magnetization for each atom
        content['starting_magnetization'] = e.list_params['starting_magnetization']
        # fixed magnetization (if present) for each atom
        content['fixed_magnetization'] = e.list_params['fixed_magnetization']

    except:
        # default if file not found and/or other parsing problems
        content['calc_finished'] = False

    res = QEOutput(name=title, content=content)

    return {title: res}

@threadit
@workit_factory(QEplugin)
def observeQE(self, code_name=None, **inp):
    """
    This function is a prototype for an observer calc for QE. It periodically checks
    the calculation output and stops it if something goes wrong according to
    the logic set by the user.
    :param: dict with FileData object; the code_name allows the user to specify
    an existing code in the database to be linked to this calc by the wrapper.
    :return: dict with the Observed instance created
    """
    uuid_list = []
    content = {}

    for key, value in inp.iteritems():
        uuid_list.append(str(value.uuid))
        if isinstance(value, FileData):
            fileinp = value

    filepath = fileinp.content['filepath']
    filecontent = fileinp.content['filecontent']

    out_file = filepath.replace('.in', '.out')

    # the observer runs the QE calculation, and then periodically looks at it
    res = run_QE_calc(self, 'Quantum Espresso 5.1.2 Mac OS X 10.9', **inp)

    # in threading, res is now a Thread object because we haven't blocked
    # with w4 yet, so we can access the output before the function actually
    # finishes the execution
    # ensures we don't check anymore if the thread finishes without problems
    while res.is_alive():
        # necessary because we need first to wait for the output
        # file to be created, before examining it
        if os.path.isfile(out_file):
            fn = open(out_file, 'r')
            hh = fn.readlines()
            fn.close()
            # dummy prototype for checking: print last 10 lines
            # in production, call the function that does the required
            # checks on the file
            print "Last 10 output lines:\n%s" % "".join(hh[-10:])
            # assume one slave per observer
            slave = self.slaves[0]
            slave.kill_itself('gra2pal', '', 'local')
        else:
            print "Waiting for output file to be created....."
        time.sleep(30)

    # get the result only when the thread has finished/has been killed
    final_res = w4(res)

    content['uuids'] =  ','.join(sorted(uuid_list))
    content['filepath'] = out_file

    # TO-DO: maybe change state to Partial according to what happened to
    # the thread. We may have a flag to change the state accordingly
    observer_res = Observed(name='observed', state='Final', content=content)

    return {observer_res.name: observer_res}


@threadit
@workit_factory(QEplugin)
def QE_workflow(self, code_name=None, **inp):
    """
    Example of workflow to run an energy calculation in Quantum Espresso. (AG)
    :param: Struc, QEInput, PPSubset, Parser.
    :return: dict with QEOutput, i.e. parsed QE output file
    """

    for key, value in inp.iteritems():
        if isinstance(value, Parser):
            parser = value

    res = w4(create_QEInput_file(self, **inp))

    # not observed version
    #out = w4(run_QE_calc(self, **res))

    # version with observer
    qe_out = w4(observeQE(self, **res))

    # parser option specified by user as input
    qe_out['parser'] = parser

    qe_out_parsed = w4(parsing_QE_output(self, 'Espresso wrapper by Zhongnan Xu', **qe_out))

    return qe_out_parsed

# for testing the workflow
if __name__ == "__main__":

    s = Session()
    # retrieve entries that are already in the DB
    # we can have functions that generate this entries automatically
    struc = s.query(Struc).filter(Struc.content['formula'].cast(String) == 'Li3').first()
    qeinp = s.query(QEInput).filter(QEInput.content['control','title'].\
        cast(String) == 'Li_C19_coarse_50Ry').first()
    ppset = s.query(PPSubset).filter(PPSubset.name == 'Li.pp_set').first()
    parser = s.query(Parser).filter(Parser.name == 'Espresso wrapper by Zhongnan Xu').first()

    c = Calc(name='QE_workflow')

    # call the workflow
    res = w4(QE_workflow(c, struc=struc, qeinp=qeinp, ppset=ppset, parser=parser))

    s.close()

    # explicitly close the job service (should be closed
    # automatically by saga anyway)
    close_js(hal_username, '', 'local')
