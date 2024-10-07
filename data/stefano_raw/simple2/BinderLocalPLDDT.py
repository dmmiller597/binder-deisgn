## Implementing the code to find binders for Albumin

from typing import List, Union
from copy import deepcopy
from biotite.database.rcsb import fetch
from biotite.structure import AtomArray

from language.constants import MULTIMER_RESIDUE_INDEX_SKIP_LENGTH

from language.program import ProgramNode
from language.sequence import ConstantSequenceSegment, FixedLengthSequenceSegment, VariableLengthSequenceSegment
from language.energy import EnergyPLDDT, EnergyPTM, EnergyBinderPAE, AverageDistance, MinimizeSurfaceHydrophobics, EnergyLocalPLDDT, EnergyHydrophobics

def create_binder( seqA: str, groupA: list, dataB: Union[str, int] ) -> ProgramNode:
    '''INPUT:
    groupA: is a list of AminoAcid on the target protein that are used to
    identify the binding spot on the protein. 
    dataB: either a starting sequence
    for the peptide (e.g., to continue optimisation) or just the length of the 
    sequence (in which case a random sequence with that length is used to start
    OUTPUT:
    A ProgramNode object specifying how to build the loss function for optimising
    the system
    '''

    #Constant part is target protein, the other is the peptide
    sequence_A = ConstantSequenceSegment( seqA )
    sequence_B = FixedLengthSequenceSegment( dataB )

    #Now basically define binding site AND the binder
    NA = len( seqA )
    NB = len( sequence_B.get() )
    print( f"Length of protein: {NA}" )
    print( f"Length of binder: {NB}" )
    first_res_B_ID = NA + MULTIMER_RESIDUE_INDEX_SKIP_LENGTH

    #This is because ESM adds "MULTIMER_RESIDUE_INDEX_SKIP_LENGTH"
    #to generate the index used for positional encoding of proteins
    #that are not attached but instead form a multimer
    groupB = range( first_res_B_ID, first_res_B_ID + NB )

    groupC = list( set( range( NA ) ) - set( groupA ) ) #Basically all AA in the protein that are not part of the binding spot

    w_ptm = 0.2 
    w_protein_plddt = 0.2 
    w_spot_plddt = 1.0
    w_pae = 2.0
    w_distance = 1.0
    w_hydrophobic = 1.0
 
    # Chains are combined into a multimer.
    return ProgramNode(
            energy_function_terms=[
            EnergyPTM(), #Maximises PTM of the whole construct
            EnergyLocalPLDDT(  group = groupC, chain_index = 0 ),#Minimise error on groupC local environment 
            EnergyLocalPLDDT(  group = groupA, chain_index = 0, verbose = False ), #Minimise error on groupA
            EnergyLocalPLDDT(  group = groupB, chain_index = 1 ), #Minimise error on groupB
            EnergyBinderPAE( groupA, groupB, mode = "limited" ), #Increases confidence in PAE between binder and binding spot
            AverageDistance( groupA, groupB ), #Penalise distance between binding spot vs binder
            EnergyHydrophobics( groupB )
            ],
            energy_function_weights=[ w_ptm, 
                                      w_protein_plddt, 
                                      w_spot_plddt, 
                                      w_spot_plddt, 
                                      w_pae, 
                                      w_distance,
                                      w_hydrophobic, 
                                      ],
            root = True,
            children=[ ProgramNode( sequence_segment = sequence_A ),
                       ProgramNode( sequence_segment = sequence_B ),
                                    ], 
            children_are_different_chains=True,
            )
