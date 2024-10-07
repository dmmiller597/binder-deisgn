import os, sys

my_path = "/gpfs/home/sangiole/Github-repos/Programmable_Allostery/" 
sys.path.append( my_path )

# Check the number of command-line arguments
if len(sys.argv) != 2:
    print( f"No path specified, using: {esm_dir}" )

try:
  # Get the input value from the command line
  esm_path = sys.argv[1]
  print( f"PATH to Allostery: {esm_path}" )
  sys.path.append( esm_path )
except:
  pass

print( f"ALL path {sys.path}" )

from copy import deepcopy
from biotite.database.rcsb import fetch
from biotite.structure import AtomArray

from language.folding_callbacks import EsmFoldv1
from language.optimize import run_simulated_tempering
from BinderLocalPLDDT import create_binder

#Decide which model you want to use and load
folding_callback = EsmFoldv1()

#In principle you can choose 8,35,150,650 and 3B (default)
#folding_callback.load(device="cpu", my_model = "8M") 
#folding_callback.load(device="cuda:0", my_model = "8M") 
#folding_callback.load(device="cuda:0", my_model = "650M") 
folding_callback.load(device="cuda:0")

#Sequence of the CD28 protein (Tcell marker)
sequence_A = "FMRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAGIVE"

dataB = 20 #This or you can provide string with sequence to start with Or integer for the string length
groupA = range( 113, 132 ) #This is just to pick a specific part, requires a list of AA
# Define the program, that is, what the final protein will optimise for:
program_Binder =create_binder( sequence_A, groupA, dataB )

#Ok, now here is the real optimisation
optimized_Binder = run_simulated_tempering(
    program = program_Binder,
    initial_temperature=1.0,
    folding_callback=folding_callback,
    display_progress=True,
    sampling_interval = 200,
    maxChanges = 2,
    pdb_file_name = "CARTCD20_",
    low_temperature = 0.1,
    high_temperature = 1.0,
    tot_sweeps = 2000,
    low_t_steps = 200, 
    high_t_steps = 200, 
    dirSave = "",
    restart_from_best = True,
    log_all_energies = True,
)

