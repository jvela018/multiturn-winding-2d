from run_onelab import run_onelab_parametric_sweep
from parameters import *



def get_gmsh(path, which_gmsh):
    gmsh = os.path.join(path, which_gmsh)
    return gmsh


def parametric_sweep_getdp():
    # frequency sweep with onelab
    gmsh = get_gmsh(gmsh_path, mygmsh)
    run_onelab_parametric_sweep(gmsh, sweep_script)




#--------------------------------------------------------------------------------------
def main(argv=None):

    #return solve_rom()
    return parametric_sweep_getdp()

if __name__ == "__main__":
    sys.exit(main() or 0)