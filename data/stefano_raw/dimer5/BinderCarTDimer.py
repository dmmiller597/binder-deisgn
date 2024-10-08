## Implementing the code to find binders for Albumin

from typing import Union

from language.constants import MULTIMER_RESIDUE_INDEX_SKIP_LENGTH

from language.program import ProgramNode
from language.sequence import ConstantSequenceSegment, FixedLengthSequenceSegment
from language.energy import EnergyPTM, GeneralEnergyBinderPAE, AverageDistance, EnergyLocalPLDDT, EnergyHydrophobics

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

    #Constant part is target protein and part of CART. Binder is the programmable part of CART
    sequence_A1 = ConstantSequenceSegment( seqA )
    sequence_A2 = ConstantSequenceSegment( seqA )
    sequence_B = FixedLengthSequenceSegment( dataB )

    #Now basically define binding site AND the binder
    NA1 = len( sequence_A1.get() )
    NA2 = len( sequence_A2.get() )
    NB = len( sequence_B.get() )
    print( f"Length of protein monomer: {NA1}" )
    print( f"Length of binder: {NB}" )
    first_res_dimer2_ID= NA1 + MULTIMER_RESIDUE_INDEX_SKIP_LENGTH
    first_res_B_ID = first_res_dimer2_ID + MULTIMER_RESIDUE_INDEX_SKIP_LENGTH + NA2

    #This is because ESM adds "MULTIMER_RESIDUE_INDEX_SKIP_LENGTH"
    #to generate the index used for positional encoding of proteins
    #that are not attached but instead form a multimer
    group_dimer_1 = list( range( NA1 ) ) 
    group_dimer_2 = list( range( first_res_dimer2_ID, first_res_dimer2_ID + NA2 ) ) 
    group_dimer1_epitope = list( groupA )
    group_dimer2_epitope = list( range( first_res_dimer2_ID + groupA[ 0 ], first_res_dimer2_ID + groupA[ -1 ] ) )
    group_dimer_epitope = group_dimer1_epitope + group_dimer2_epitope
    groupB = range( first_res_B_ID, first_res_B_ID + NB )

    w_ptm = 0.2 
    w_spot_plddt = 1.0
    w_pae = 2.0
    w_distance = 1.0
    w_hydrophobic = 1.0
 
    # Chains are combined into a multimer.
    return ProgramNode(
            energy_function_terms=[
            EnergyPTM(), #Maximises PTM of the whole construct
            EnergyLocalPLDDT(  group = group_dimer_1, chain_index = 0, verbose = False ),#Minimise error on group_dimer1_epitope local environment 
            EnergyLocalPLDDT(  group = group_dimer_2, chain_index = 1, verbose = False ),#Minimise error on group_dimer2_epitope local environment 
            EnergyLocalPLDDT(  group = groupB, chain_index = 2, verbose = False ), #Minimise error on groupB / binder region
            GeneralEnergyBinderPAE( group1 = group_dimer_1, group2 = group_dimer_2, chain_1_num = 0, chain_2_num = 1, mode = "limited"),
            GeneralEnergyBinderPAE( group1 = group_dimer1_epitope, group2 = groupB, chain_1_num = 0, chain_2_num = 2, mode = "limited"),
            GeneralEnergyBinderPAE( group1 = group_dimer2_epitope, group2 = groupB, chain_1_num = 1, chain_2_num = 2, mode = "limited"),
            AverageDistance( group_dimer_epitope, groupB, debug = True ), #Penalise distance between binding spot vs binder
            EnergyHydrophobics( groupB )
            ],
            energy_function_weights=[ w_ptm, 
                                      w_spot_plddt, 
                                      w_spot_plddt, 
                                      w_spot_plddt, 
                                      w_pae, 
                                      w_pae, 
                                      w_pae, 
                                      w_distance,
                                      w_hydrophobic, 
                                      ],
            root = True,
            children=[  ProgramNode( sequence_segment = sequence_A1 ),
                        ProgramNode( sequence_segment = sequence_A2 ),
                        ProgramNode( sequence_segment = sequence_B ),
                                                  ],
            children_are_different_chains=True,
            )
