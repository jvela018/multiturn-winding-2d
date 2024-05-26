import os

def run_onelab_parametric_sweep(gmsh, sweep_script):
    os.system(gmsh + ' ' + sweep_script+'.py' + ' - -v 0 ')

    del_script = sweep_script + '.db'
    if os.path.isfile(del_script) is True:
        os.system('rm -r '+ del_script)

