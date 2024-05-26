import os
import sys

#series
#path1 = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/python_results/multi-turn-voltage/'
#path2 = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/python_results/multi-turn-voltage/'
#parallel
path1 = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/python_results/three_thin_wire/'
path2 = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/python_results/three_thin_wire/'


# For frequency sweep files

str_freq = '1'
str_distance = ' '
str_sleeve = '0.001'
str_radius = '0.001'
local = False

#file_to_read = 'impedance_frequency'
#file_reference = path1 + 'full_' + file_to_read + '_Rc_' + str_radius + '_rs_' + str_sleeve + '.dat'
#file_comparison = path2 + 'thin_' + file_to_read + '_Rc_' + str_radius + '_rs_' + str_sleeve + '.dat'
#unchanged_column = 5

#file_to_write = 'error_impedance_frequency'
#ofile    = path1 + file_to_write + '_Rc_' + str_radius + '_rs_' + str_sleeve + '.dat'

#file_to_read = 'full_distance'

#Distance
'''
freq_param = '1000000'
sleeve_type = 'unstructured'
connection_type = 'series'
nturns = '15'

sleeve_param = '_0.5mm'

if connection_type == 'parallel':
    file_reference = path1 + 'full_distance-sweep_freq_' + freq_param + '_' + connection_type + '.dat'
    file_comparison = path2 + 'thin_distance-sweep_freq_' + freq_param + '_' + sleeve_type + '_' + connection_type + sleeve_param +'.dat'
    file_to_write = 'error_distance_freq_' + freq_param + '_' + sleeve_type + '_' + connection_type + sleeve_param
else:
    file_reference = path1 + 'full_distance-sweep_freq_' + freq_param + '_' + connection_type + '_' + nturns + 'turns.dat'
    file_comparison = path2 + 'thin_distance-sweep_freq_' + freq_param + '_' + sleeve_type + '_' + connection_type + '_' + nturns + 'turns'+ sleeve_param +'.dat'
    file_to_write = 'error_distance_freq_' + freq_param + '_' + sleeve_type + '_' + connection_type + '_' + nturns + 'turns'+ sleeve_param

unchanged_column = 1   # this is for series connection change column manually
ofile    = path1 + file_to_write + '.dat'
'''


#Frequency
'''
dist_param = '0.008'
sleeve_type = 'unstructured'
connection_type = 'series'
nturns = '5'

sleeve_param = '_0.5mm'

if connection_type == 'parallel':
    file_reference = path1 + 'full_frequency-sweep_dist_' + dist_param + '_' + connection_type + '_' + nturns + 'turns.dat'
    file_comparison = path2 + 'thin_frequency-sweep_dist_' + dist_param + '_' + sleeve_type + '_' + connection_type + '_' + nturns + 'turns'+  sleeve_param +'.dat'
    file_to_write = 'error_distance_dist_' + dist_param + '_' + sleeve_type + '_' + connection_type + sleeve_param
else:
    file_reference = path1 + 'full_frequency-sweep_dist_' + dist_param + '_' + connection_type + '_' + nturns + 'turns.dat'
    file_comparison = path2 + 'thin_frequency-sweep_dist_' + dist_param + '_' + sleeve_type + '_' + connection_type + '_' + nturns + 'turns'+  sleeve_param +'.dat'
    file_to_write = 'error_distance_dist_' + dist_param + '_' + sleeve_type + '_' + connection_type + '_turns' + nturns + sleeve_param

unchanged_column = 0
ofile    = path1 + file_to_write + '.dat'
'''

'''
# comparison for impedance vs distance
str_freq = '1000000'
str_distance = ' '
str_sleeve = '0.001'
str_radius = '0.001'
local = False

file_to_read = 'impedance_distance'
file_reference = path1 + 'full_' + file_to_read + '_freq_' + str_freq + '_Rc_' + str_radius + '_rs_' + str_sleeve + '.dat'
file_comparison = path2 +'thin_' + file_to_read + '_freq_' + str_freq + '_Rc_' + str_radius + '_rs_' + str_sleeve + '.dat'
unchanged_column = 5

file_to_write = 'error_impedance_distance'
ofile    = path1 + file_to_write + '_freq_' + str_freq + '_Rc_' + str_radius + '_rs_' + str_sleeve + '.dat'
'''



# ---------------------------------
#For local comparison
# ---------------------------------

#str_freq = '1'
#str_freq = '1e+06'
#str_distance = '0.04'
#str_sleeve = '0.003' # line
#str_sleeve = '0.001' # noline
#str_radius = '0.001'
#str_outline = '_line'
local = True
'''
#for finite element model
'''
#path1_local = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/getdp_results_full_2D/'
#path2_local = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/getdp_results_thin_2D/'
#path3_local ='/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/Bundles/Bundle_2D_clean/python_results/three_thin_wire/'

path1_local = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/skin_prox/getdp_results_full_2D/'
path2_local = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/skin_prox/getdp_results_thin_2D/'
path3_local ='/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/skin_prox/python_results/three_thin_wire/'

dist_parame = '0.00205' #'0.008' or '0.00205'
freq_parame = '1'   #'1' or '1e6'
sleeve_type = 'structured'
sleeve_radpar = '_0.5mm'  # _0.5mm or _3mm and for 1mm just ''

file_reference = path1_local  + 'Cut_az_dist_' + dist_parame + '_' + freq_parame + 'Hz.dat'
file_comparison = path2_local + 'Cut_az_thin-' + sleeve_type + '_dist_' + dist_parame + '_' + freq_parame + 'Hz' + sleeve_radpar + '.dat'
unchanged_column = 0

file_to_write = 'error_az_' + sleeve_type + '_dist_' + dist_parame + '_' + freq_parame + 'Hz' + sleeve_radpar
ofile    = path3_local + file_to_write + '.dat'

# remove file if already exists
if os.path.isfile(ofile) is True:
    os.system('rm -r ' + ofile)

'''
#analytical numerical
path1_local = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/skin_proximity/coeff/'
path2_local = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/skin_proximity/coeff/'
path3_local = '/Users/jonathanvelasco/Documents/Validation_Models/my_phd_models/skin_proximity/coeff/'

file_reference = path1_local  + 'jz_line_numerical.dat'
file_comparison = path2_local + 'jz_line_analytical.dat'
unchanged_column = 1

file_to_write = 'error_analytical_numerical'
ofile    = path3_local + file_to_write + '.dat'
'''


def get_error(fref,fcomp, maxval_list):

    print('Comparing files...')
    print('File1:', file_reference)
    print('File2:', file_comparison)

    file2read_1 = open(fref, 'r')
    file2read_2 = open(fcomp, 'r')

    for n, (line1, line2) in enumerate(zip(file2read_1, file2read_2)):

        line1_len = len(line1.strip().split())

        # write headers onto file
        if n == 0 and local is False:
            file2write = open(ofile, 'w')
            file2write.write(line1)
            file2write.close()

        # read values and calculate error
        if n>0 or local is True:
            error_list = []
            for i in range(line1_len):
                val_ref = float(line1.strip().split()[i])
                val_comp = float(line2.strip().split()[i])
                val_refmax = maxval_list[i]
                #print(i, val_refmax > val_ref)
                try:
                    error_val = str((abs(val_comp - val_ref)/abs(val_refmax))*100)
                except:
                    error_val = 'nan'


                print('These are the values (refval,compval,maxref,error):', val_ref, val_comp, val_refmax, error_val)

                # store errors per line into list and append to file
                if(i == unchanged_column ):
                    error_list.append(str(val_ref))  # keep last value as it's the sweeping value
                else:
                    error_list.append(error_val)
            file2write = open(ofile, 'a')
            print('\t'.join(error_list), file=file2write)


    file2write.close()
    file2read_2.close()
    file2read_1.close()
    print("Data has been written in file", ofile)

def get_maxref(fref):
    print('Maximum reference value in...')
    print('File:', file_reference)

    # get length of line and initialize zeros list
    file2read_1 = open(fref, 'r')
    file2read_1.readline().strip().split()
    line1_len = len(file2read_1.readline().strip().split())
    val_maxlist = [0.0]*line1_len
    file2read_1.close()

    # read line by line and get the max value of each column
    file2read_1 = open(fref, 'r')

    for n,line1 in enumerate(file2read_1):

        line1_len = len(line1.strip().split())
        # read values and calculate error
        if n > 0 or local is True:
            for i in range(line1_len):
                val_ref = float(line1.strip().split()[i])
                max_valref = val_maxlist[i]
                if (val_ref > max_valref):
                    val_maxlist[i] = val_ref
    file2read_1.close()

    # Verify that these are actually max values
    '''
    file2read_1 = open(fref, 'r')
    print(val_maxlist)
    for n,line1 in enumerate(file2read_1):
        line1_len = len(line1.strip().split())
        # read values and calculate error
        if n > 0 or local is True:
            for i in range(line1_len):
                val_ref = float(line1.strip().split()[i])
                print(i,val_maxlist[i] >= val_ref)
    file2read_1.close()
    '''

    return val_maxlist


def main(argv=None):

    maxval_list = get_maxref(file_reference)
    get_error(file_reference, file_comparison, maxval_list)


if __name__ == "__main__":
    sys.exit(main() or 0)