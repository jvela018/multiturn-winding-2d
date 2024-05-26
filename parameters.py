# ---------------------------------------paths to software packages-----------------------------------------------------

# getdp and gmsh paths
getdp_path = '/Users/jonathanvelasco/Documents/ONELAB/getdp-git-MacOSX/bin/'
gmsh_path = '/Users/jonathanvelasco/Documents/ONELAB/gmsh/build/'
onelab_path = '/Users/jonathanvelasco/Documents/ONELAB'
mygmsh  = 'gmsh'
mygetdp = 'getdp'

# ---------------------------------------------import ONELAB------------------------------------------------------------

import os
import sys
# import onelab module
sys.path.append(onelab_path)
import onelab

# -----------------------------------------paths to model to solve------------------------------------------------------

# model name and path
model_name = 'multiturn_2D'
model_path = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/'
matrix_path = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/'
sweep_script = 'inductance_matrix'

# ------------------------------------------------parameters------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# Gmsh and GetDp settings
# --------------------------------------------------------------------------------------------------------------
# gmsh meshing options
gmsh_verbosity = ''

# getdp options: resolution
getdp_command = 'solve'
getdp_resolution = 'Resolutions'
getdp_verbosity = '0'

# post-processing commands
getdp_post_cmd = ''
getdp_postoperation = 'PostOperations'

# --------------------------------------------------------------------------------------------------------------
# Model settings
# --------------------------------------------------------------------------------------------------------------
model_dimension = '2'         # 2 or 3
sweep_variable = 'frequency'   # frequency or distance
source_type = 'voltage'       # current or voltage
solve_method = 'full'         # full or thin
sleeve_type = 'structured'  # structured or unstructured
connection_type = 'series'    # series or parallel (series = conductors together as turns, parallel = individual wires
model_shape = 'cylinder'      # cylinder or sphere
mesh_size = 'coarse'          # coarse, medium, fine
adaptive_mesh = True          # this is not tied parameters.pro manually change freqRefinement.

# number of rows and columns populated by wires
nrows    = 5  # 10 full  5
ncolumns = 1   # 2 full  1,2,3

# --------------------------------------------------------------------------------------------------------------
# ONELAB sever names
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------static severnames-------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# Geometry parameters
# --------------------------------------------------------------------------------------------------------------
servername_Rc = "Input/0Geometry Parameters/2Conductor Radius (R_c)/(m)"
servername_rs = "Input/0Geometry Parameters/3Sleeve Radius (rs)/(m)"
servername_distx = "Input/0Geometry Parameters/4Distance Between Conductors in X/(m)"
servername_disty = "Input/0Geometry Parameters/5Distance Between Conductors in Y/(m)"
servername_rows = "Input/0Geometry Parameters/7Coil Bundle Array/0Rows"
servername_columns = "Input/0Geometry Parameters/7Coil Bundle Array/1Columns"
servername_nturns = "Input/0Geometry Parameters/7Coil Bundle Array/2Number of Turns"

# --------------------------------------------------------------------------------------------------------------
# Material Properties and Frequency parameters
# --------------------------------------------------------------------------------------------------------------

servername_sig = "Input/1Conductor Material/ Electrical Conductivity (S\m)/"
servername_mur = "Input/1Conductor Material/ Relative Permeability/"
servername_epr = "Input/1Conductor Material/ Relative Permittivity/"
servername_freq = "Input/2Model Excitation/3Frequency (Hz)/ "
servername_Vp = "Input/2Model Excitation/1Voltage (V)/Wire %g"
servername_Ip = "Input/2Model Excitation/2Current (A)/Wire %g"
servername_Vg = "Input/2Model Excitation/1Voltage Ground (V)/Wire %g"
servername_excitation = "Input/2Model Excitation/0Excitation Type"
servername_connectiontype = "Input/2Model Excitation/4Feature Checkboxes/Connect in Series"

# --------------------------------------------------------------------------------------------------------------
# Model Options
# --------------------------------------------------------------------------------------------------------------
servername_model = "Input/0Model Configuration/0Type"

# --------------------------------------------------------------------------------------------------------------
# Particular test cases:
# --------------------------------------------------------------------------------------------------------------
servername_condOutline = "Input/1Geometry Parameters/2Feature Checkboxes/Show Conductor Outline"
servername_dimension = "Input/1Geometry Parameters/0Geometrical Dimension"
servername_curve = "Input/1Geometry Parameters/1Model Shape"
servername_ied = "Input/1Geometry Parameters/2Feature Checkboxes/Use Shell Transformation"
servername_bdl = "Input/1Geometry Parameters/2Feature Checkboxes/Use Boundary Layer in Conductor"
servername_meshtype = "Input/2Meshing Parameters/0Mesh Type"
servername_meshsize = "Input/2Meshing Parameters/2Mesh Size"

# --------------------------------------------------------------------------------------------------------------
# Output variables
# --------------------------------------------------------------------------------------------------------------

if solve_method == 'thin':
    smethod_str = '1Corr'

    if connection_type == 'series':
        servername_resistance = 'Output/0Impedance/0Resistance at Source ' \
                                '(' + model_dimension + 'D)/' + smethod_str + '-' + sleeve_type.lower()
        servername_inductance = 'Output/0Impedance/1Inductance at Source ' \
                                '(' + model_dimension + 'D)/' + smethod_str + '-' + sleeve_type.lower()
    else:
        servername_resistance = 'Output/0Impedance/0Resistance ' \
                                '(' + model_dimension + 'D)/' + smethod_str + '-' + sleeve_type.lower()
        servername_inductance = 'Output/0Impedance/1Inductance ' \
                                '(' + model_dimension + 'D)/' + smethod_str + '-' + sleeve_type.lower()


else:
    smethod_str = '0Full'

    if connection_type == 'series':
        servername_resistance = 'Output/0Impedance/0Resistance at Source (' + model_dimension + 'D)/' + smethod_str
        servername_inductance = 'Output/0Impedance/1Inductance at Source (' + model_dimension + 'D)/' + smethod_str
    else:
        servername_resistance = 'Output/0Impedance/0Resistance (' + model_dimension + 'D)/' + smethod_str
        servername_inductance = 'Output/0Impedance/1Inductance (' + model_dimension + 'D)/' + smethod_str

# --------------------------------------------dynamic severnames-------------------------------------------------

if source_type in 'current':
    servername_excitation = 'Input/2Model Excitation/2Current (A)/Wire '
else:
    servername_excitation = 'Input/2Model Excitation/1Voltage (V)/Wire '



''' THESE ARE FOR IMPEDANCE MATRIX
servername_resistance = '}Resistance/'  # note space after Wire
servername_inductance = '}Inductance/'
'''''


servername_resistance_p = '}Resistance/'



###