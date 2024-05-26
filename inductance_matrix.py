from parameters import *



# folders containing getdp results and the one using for sweeping
getdp_res_dir = 'getdp_results'
sweep_script = 'parametric_sweep'

# lists to declare variable parameter

# operating frequency
#freq_list = [1e7, 2e8, 3e9]
#freq_list = np.linspace(0.01e6, 100e6, 60).tolist()

#frequency sweep

freq_list = []
f_initial = 1
f_last = 35
f_step = 0.5


for i in range(f_last):
    freq_list += [f_initial]
    f_initial *= 1+f_step


#freq_list = [1e5]

radius_list = [1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3]

#distance_list = [0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.00205]

distance_list = [0.008, 0.006, 0.004, 0.0038, 0.0036, 0.0034, 0.0032, 0.003, 0.0028, 0.0026, 0.0024, 0.0022, 0.0021, 0.00205]
#distance_list = [0.004, 0.00205]
#distance_list = [0.004, 0.00205]


radius_divisions = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]

# ----------------------------------------------------------------------------------
# sweep control parameters
# -----------------------------------------------------------------------------------

#

#sweep_variable = 'frequency'
rs = 0.001
# -----------------------------------
# DON'T CHANGE ANYTHING BELOW
# -----------------------------------

if sweep_variable == 'distance':

    # model dimension for meshing (2 or 3)
    #str_model_dimension = '2'
    remesh = True # if sweep doesn't need remeshing False otherwise True
                # e.g. Frequency sweep

    # static value during sweep (if not setup then it .pro values will be used)
    input_static_name = 'Frequency'
    input_static_units = '[Hz]'
    input_static_servername = servername_freq
    input_static_value = int(1)  #freq_list[0]

    # input variables to be swept
    str_varname = 'Distance'
    unit_scale = 1e-3
    str_var_units = '[mm]'
    input_servername = servername_disty
    sweep_list = distance_list
else:

    # model dimension for meshing (2 or 3)
    #str_model_dimension = '2'
    # if sweep doesn't need remeshing False otherwise True
    # e.g. Frequency sweep
    if adaptive_mesh is True:
        remesh = True
    else:
        remesh = False
    sweep_list = freq_list
    str_varname = 'Frequency'
    str_var_units = '[Hz]'
    input_servername = servername_freq

    input_static_servername = 'Input/0Geometry Parameters/3Sleeve rs Radius/[m]'
    input_static_value = rs  # freq_list[0]


    '''
    # static value during sweep (if not setup then it .pro values will be used)
    input_static_name = 'Sleeve radius'
    input_static_units = '[m]'
    

    # input variables to be swept
    str_varname = 'Frequency'
    unit_scale = 1
    str_var_units = '[Hz]'
    '''

str_ovarname = 'Impedance'



# --------------------------------------------------------------------------------------------




# -----------------functions------------------------------------------------------------------

def get_onelab(gmsh_path, mygmsh, getdp_path, mygetdp):
    # create a new onelab client
    c = onelab.client(__file__)

    # get Gmsh and GetDP locations from Gmsh options
    gmsh  = os.path.join(gmsh_path, mygmsh)
    getdp = os.path.join(getdp_path, mygetdp)
    c.sendInfo('Using packages: \n gmsh={0} \n getdp={1}'.format(gmsh, getdp))
    return c, gmsh, getdp



def get_model(c, gmsh, model_path, model_name):
    # create a onelab variable for the model name
    model = c.defineString(model_name, value=model_name)

    # we're done if we don't do the actual calculation
    if c.action == 'check':
        exit(0)

    # get model file names with correct path
    model_geo = c.getPath(model_path + model_name + '.geo')
    model_msh = c.getPath(model_path + model_name + '.msh')
    model_pro = c.getPath(model_path + model_name + '.pro')


    # check if model exists
    print("--------------------Active Model-------------------------")
    print('Geometry file {0} found: {1}'.format(model_name + '.geo', os.path.isfile(model_geo)))
    print('Solve file {0} found: {1}'.format(model_name + '.pro', os.path.isfile(model_pro)))
    print('Mesh file {0} found: {1}'.format(model_name + '.msh', os.path.isfile(model_msh)))
    print("---------------------------------------------------------")

    # if .pro or .geo are not found exit simulation
    if os.path.isfile(model_geo) is False or os.path.isfile(model_pro) is False:
        print('No input model. Please make sure to have a .pro and a .geo file on path {0}'.format(model_path))
        sys.exit(0)

    # mesh if .msh not available
    if os.path.isfile(model_msh) is False:
        print("---------------------------------------------------------")
        print("Meshing {0}D geometry {1}.geo with Gmsh ...".format(str(model_dimension),model_name))
        print("---------------------------------------------------------")
        gmsh_options = get_gmsh_options()
        run_gmsh(c, gmsh, model_geo, gmsh_options)

    return model_geo, model_msh, model_pro

# Gmsh related functions


def get_gmsh_options():
    gmsh_opt_mesh = '-' + str(model_dimension)

    if str(gmsh_verbosity) not in '':
        gmsh_opt_verbosity = '-v ' + str(gmsh_verbosity)
    else:
        gmsh_opt_verbosity = str(gmsh_verbosity)

    gmsh_options = gmsh_opt_mesh + ' ' + gmsh_opt_verbosity
    return gmsh_options


def run_gmsh(c, gmsh, model_geo, gmsh_options):
    c.runSubClient('myGmsh', gmsh + ' ' + model_geo + ' ' + gmsh_options)


def remesh_gmsh(c, gmsh, model_geo, model_msh):

    # remove current mesh if exists
    if os.path.isfile(model_msh) is True:
        os.system('rm -r ' + model_msh)
    # mesh
    gmsh_options = get_gmsh_options()
    run_gmsh(c, gmsh, model_geo, gmsh_options)

# GetDP related functions


def run_getdp(c, getdp, model_pro, model_msh):
    getdp_options = get_getdp_options()
    c.runSubClient('myGetDP', getdp + ' ' + model_pro + ' -msh ' + model_msh + ' ' + getdp_options)


def get_getdp_options_resolution():
    getdp_opt_cmd = '-' + str(getdp_command)
    getdp_opt_resolution = ' ' + str(getdp_resolution)

    if str(getdp_verbosity) not in '':
        getdp_opt_verbosity = ' -v ' + str(getdp_verbosity)
    else:
        getdp_opt_verbosity = ' ' + str(getdp_verbosity)

    getdp_options = getdp_opt_cmd + getdp_opt_resolution + getdp_opt_verbosity
    return getdp_options


def get_getdp_options_post():
    getdp_opt_post_cmd = '-' + str(getdp_post_cmd)
    getdp_opt_post = ' ' + str(getdp_postoperation)

    getdp_options = getdp_opt_post_cmd + getdp_opt_post
    return getdp_options


def get_getdp_options():
    getdp_options = get_getdp_options_resolution()
    getdp_options += get_getdp_options_post()
    return getdp_options

def set_variable(c, var_servername, var_value):
    c.setNumber(var_servername, value=var_value)


def set_method(c, solve_method):

    if solve_method is not 'thin':
        print('Solving:', str_varname, 'sweep with Full Model')
        isthin_val = 0
    else:
        print('Solving:', str_varname, 'sweep with Thin Wire Model')
        isthin_val = 2

    isthin = c.defineNumber(servername_model, value=isthin_val, choices=[0, 1, 2])
    #set_variable(c, servername_model, isthin_val)
    #print(c.getNumber(servername_model))


def set_numwires(c, nrows, ncolumns):
    set_variable(c, servername_columns, ncolumns)
    set_variable(c, servername_rows, nrows)
    numwires = nrows*ncolumns
    print('This model has', numwires, 'wires')
    return numwires

def set_connectiontype(c, connection_type):

    if connection_type == 'parallel':
        con_type = 0
    else:
        con_type = 1
    set_variable(c, servername_connectiontype, con_type)
    print('Connection type: ', connection_type)


def set_model_dimension(c, model_dimension):

    if (model_dimension == '2'):
        print('Setting up 2D model')
        model_menu_val = 0
    else:
        print('Setting up 3D model')
        model_menu_val = 1

    set_variable(c, servername_dimension, model_menu_val)



def calculate_impedance_matrix(c, getdp, model_pro, model_msh):

    # create dictionaries to store impedance labels and vales
    resdict = {}
    resdict_prox = {}
    inddict = {}
    ################
    resistance_labels_str = []
    resistance_prox_labels_str = []
    inductance_labels_str = []
    ################

    # set number of columns in model (if connected in series it's turns)
    numwires = set_numwires(c, nrows, ncolumns)

    # create input server names
    for j in range(1, numwires + 1):

        # activate and deactivate terminals
        set_variable(c, servername_excitation + str(j), 1)
        for k in range(1, numwires + 1):
            if k != j:
                print('Setting feed:', k, 'to zero')
                set_variable(c, servername_excitation + str(k), 0)

        run_getdp(c, getdp, model_pro, model_msh)

        for i in range(1, numwires + 1):

            ext_ij = str(i) + str(j)

            # Get self part of the matrix (Z11, Z22... i=j)
            if i == j:
                resistance = c.getNumber(servername_resistance + ext_ij)
            elif (i != j) and (solve_method == 'thin'):
                resistance_prox = c.getNumber(servername_resistance + '-Prox'+ ext_ij)
                # print(servername_resistance_p + smethod_str + '-Prox'+ ext_ij)

            inductance = c.getNumber(servername_inductance + ext_ij)


            for k in range(1, numwires + 1):
                print('Feed', k, ':',
                        c.getNumber(servername_excitation + str(k)))

            ###################
            print('Calculating', str_ovarname.lower())  # 'Z = R + jX '
            if i == j:
                print('R' + ext_ij + ':', resistance)
            elif (i != j) and (solve_method == 'thin'):
                print('R' + ext_ij + '(prox):', resistance_prox)

            print('L' + ext_ij + ':', inductance)
            print('---------------------------------------------------------')
            #################

            # populate dictionaries with their labels
            if i == j:
                resistance_key = 'R' + str(j)
                if resistance_key not in resdict.keys():
                    resdict[resistance_key] = resistance
                    ################
                    resistance_labels_str.append(resistance_key)
                    #################
            elif (i != j) and (solve_method == 'thin'):
                resistance_prox_key = 'R' + ext_ij + '(prox)'
                if resistance_prox_key not in resdict_prox.keys():
                    resdict_prox[resistance_prox_key] = resistance_prox
                #################
                resistance_prox_labels_str.append(resistance_prox_key)
                #################

            inductance_key = 'L' + str(i) + str(j)
            if inductance_key not in inddict.keys():
                inddict[inductance_key] = inductance
            inductance_labels_str.append(inductance_key)
            ################


    # we are returning a a string to be pasted onto the output file not the actual value.
    # the actual values if needed are in the dictionary

    resistance_value_str = ''
    resistance_label_str = ''
    for val in resistance_labels_str:
        resistance_label_str += '\t' + str(val)
        resistance_value_str += '\t' + str(resdict[val])

    inductance_value_str = ''
    inductance_label_str = ''
    for val in inductance_labels_str:
        inductance_label_str += '\t' + str(val)
        inductance_value_str += '\t' + str(inddict[val])

    resistance_prox_value_str = ''
    resistance_prox_label_str = ''
    if resistance_prox_labels_str:
        for val in resistance_prox_labels_str:
            resistance_prox_label_str += '\t' + str(val)
            resistance_prox_value_str += '\t' + str(resdict_prox[val])


    return resistance_label_str, resistance_value_str, inductance_label_str, inductance_value_str, \
           resistance_prox_label_str, resistance_prox_value_str


def variable_sweep_getdp(c, gmsh, getdp, model_geo, model_pro, model_msh, variable_list):

    print('Solving platform: GetDP')
    print(str_varname, 'sweep for range ', [(str(var) + str_var_units) for var in variable_list])

    for n, var in enumerate(variable_list):

        print('---------------------------------------------------------')
        print('Simulation: {0}/{1}: {2} {3}{4}'.format(n+1,len(variable_list), str_varname,var, str_var_units))


        #set variable user parameter
        set_variable(c, input_servername, var)

        #equal spacing in y (top) and x (bottom)
        if sweep_variable == 'distance':
            set_variable(c, servername_distx, var)


        #set static parameter
        set_variable(c, input_static_servername, input_static_value)

        if remesh is True:
            # remesh according to solve methods
            remesh_gmsh(c, gmsh, model_geo, model_msh)

        # run_getdp(c, getdp, model_pro, model_msh)

        # retrieve inductance matrix (note that these are strings values are inside function)

        # Impedance matrix has been commented out!!!!

        if connection_type == 'parallel':
            resistance_label_str, resistance_value_str, inductance_label_str, inductance_value_str, \
                resistance_prox_label_str, resistance_prox_value_str = \
            calculate_impedance_matrix(c, getdp, model_pro, model_msh)
        else:
            resistance_label_str = '\t' + 'Resitance'
            inductance_label_str = '\t' + 'Inductance'
            run_getdp(c, getdp, model_pro, model_msh)

            resistance_prox_label_str = ''
            resistance_prox_value_str = ''
            resistance_value = c.getNumber(servername_resistance)
            inductance_value = c.getNumber(servername_inductance)
            resistance_value_str = '\t' + str(resistance_value)
            inductance_value_str = '\t' + str(inductance_value)

            print('Resistance p.u.', resistance_value)
            print('Inductance p.u.', inductance_value)


        # get frequency, R_c, r_s and distance to use on file name
        active_freq = c.getNumber(servername_freq)
        active_rs = c.getNumber(servername_rs)
        active_disty = c.getNumber(servername_disty)

        # create strings to label output file
        if sweep_variable is not 'frequency':
            freq_str = '_freq_' + str(int(active_freq))
        else:
            freq_str = ''

        if sweep_variable is not 'distance':
            dist_str = '_dist_' + str(active_disty)
        else:
            dist_str = ''

        if sweep_variable is not 'sleeve':
            sleeve_str = '_sleeve_' + str(active_rs)
        else:
            sleeve_str = ''


        # create folder called python_results if not there
        python_results_dir = 'python_results'
        if os.path.isdir(python_results_dir) is False:
            os.mkdir(python_results_dir)
        # create file to store getdp post processing
        # if it's thin the file will append the type of sleeve (structured/unstructured)
        if solve_method is not 'thin':
            o_filename = os.path.join(python_results_dir, str(solve_method) + '_' +
                                      (sweep_variable + '-sweep').lower() +
                                      freq_str +
                                      dist_str + '_' + connection_type + '_' + str(nrows*ncolumns) + 'turns' +
                                      '.dat')
        else:
            o_filename = os.path.join(python_results_dir, str(solve_method) + '_' +
                                      (sweep_variable + '-sweep').lower() +
                                      freq_str +
                                      dist_str +
                                      '_' + sleeve_type + '_' + connection_type + '_' + str(nrows*ncolumns) + 'turns' +
                                      '.dat')

        if n == 0 and os.path.isfile(o_filename) is True:
            os.system('rm -r ' + o_filename)


        if n == 0:
            file2write = open(o_filename, 'w')

            if resistance_prox_label_str:
                print('#' + 'Frequency' + '\t' + 'Distance' + '\t' + 'Sleeve Radius'
                      + resistance_label_str
                      + inductance_label_str
                      + resistance_prox_label_str, file=file2write)
            else:
                print('#' + 'Frequency' + '\t' + 'Distance' + '\t' + 'Sleeve Radius'
                      + resistance_label_str
                      + inductance_label_str, file=file2write)


            file2write.close()


        file2write = open(o_filename,'a')

        if resistance_prox_label_str:
            print(str(active_freq) + '\t' + str(active_disty) + '\t' + str(active_rs) +
                  resistance_value_str + inductance_value_str + resistance_prox_value_str, file=file2write)
        else:
            print(str(active_freq) + '\t' + str(active_disty) + '\t' + str(active_rs) +
                  resistance_value_str + inductance_value_str, file=file2write)



    file2write.close()


# --------------------------------------------------------------------------------------------
def main(argv=None):

    # get onelab instance and paths for gmsh and getdp
    c, gmsh, getdp = get_onelab(gmsh_path, mygmsh, getdp_path, mygetdp)

    # setup paths and files for model to solve
    model_geo, model_msh, model_pro = get_model(c, gmsh, model_path, model_name)

    # set model dimensions (2D or 3D)
    set_model_dimension(c, model_dimension)

    # set connection type: parallel or series (multi-turn)
    set_connectiontype(c, connection_type)

    # setup solve method: full or thin
    set_method(c, solve_method)

    # set number of wires in model
    numwires = set_numwires(c, nrows, ncolumns)

    # remesh according to solve methods
    remesh_gmsh(c, gmsh, model_geo, model_msh)

    # perform parametric sweep over user defined values values
    #calculate_impedance_matrix(c, getdp, model_pro, model_msh)
    variable_sweep_getdp(c, gmsh, getdp, model_geo, model_pro, model_msh, sweep_list)


if __name__ == "__main__":
    sys.exit(main() or 0)

# --------------------------------------------------------------------------------------------
