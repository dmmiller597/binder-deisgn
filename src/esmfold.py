import torch
from transformers import EsmForProteinFolding, EsmTokenizer, AutoTokenizer

import time
from transformers import AutoTokenizer
from typing import List
from dataclasses import dataclass

from abc import ABC, abstractmethod

from biotite.structure import AtomArray

from io import StringIO

import numpy as np

from .protein import output_to_pdb, pdb_file_to_atomarray

from .residue_constants import atom_order

@dataclass
class FoldingResult:
    atoms: AtomArray
    ptm: float
    plddt: float
    pae: float #Added by Stefano
    residue_index: list #Added by Stefano
    pdb: str #Added by Stefano
    local_plddt: float


class FoldingCallback(ABC):
    "Interface for running ESMFold and other folding methods."

    def __init__(self) -> None:
        pass

    @abstractmethod
    def load(self, device: str) -> None:
        pass

    @abstractmethod
    def fold(self, sequence: str, residue_indices: List[int]) -> FoldingResult:
        pass



#Modified by Stefano Angioletti-Uberti to include additional outputs/features
class EsmFoldv1(FoldingCallback):
    "Runs ESMFold v1.0."

    def __init__(self) -> None:
        super().__init__()

        self.model = None
        self.i = 0

    def load( self, device: str, memory_light = True, backend: str = None) -> None:
        """Load the correct model from pretrained one. if memory_light is True, then
        also uses an optimisation that makes the model potentially slower but
        uses less memory"""

        self.tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")

        self.model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")

        if backend is not None:
            import time
            t0 = time.time()
            self.model = torch.compile(self.model, backend=backend) #This is to optimise to use pytorch2.0, if 1.*, comment this line
            t1 = time.time()
            print( f"Time to compile model: {t1-t0} for backend {backend}")

        if memory_light:
          #Some optimisations to reduce memory usage
          self.model.esm = self.model.esm.half()

        self.model.to(device)
        self.model.eval() #This sets the model in evaluation mode
      
        # Check the total GPU memory allocated in bytes
        allocated_memory_gib = torch.cuda.memory_allocated()/ (1024 ** 3)
        print(f"Total GPU memory allocated: {allocated_memory_gib:.2f} GiB")

    def fold(self, 
             sequence: str, 
             residue_indices: List[int], 
             print_pdb : bool = True,
             pdb_file : str = "Prediction.pdb", 
             my_device = "cuda:0",
             ) -> FoldingResult:
        assert self.model is not None, "Must call load() before fold()."

        with torch.no_grad():
          device = torch.device( my_device )

          # Stefano Added comment:
          # the method to set residue_indices in Program makes them start from one, instead of zero.
          # The reason why I subtract one here is that I want to make sure I have exactly the positional encodings
          # used by Meta's team and don't want to change it even one bit because this is how the ESMFold model was trained,
          # with positional encodings starting from zero and not 1 (at least this is the case based on Meta's own code)

          # The below comment from Meta's team is somewhat misplaced. Yes, if residx start from 1 output_to_pdb starts from 2. But the 
          # fact remains that even correcting this (which is done in my version of output_to_pdb found in external.py)
          # you still want positional encodings starting from 0 as in the original trained network. Not sure how big of 
          # a difference it would make, but still, we leave it as in the original case.

          # ORIGINAL COMMENT BY META's TEAM
          # TODO: Current `esm.esmfold.v1.misc.output_to_pdb()` adds 1 to the `residx`
          # mistakenly, just subtract 1 for now but fix in a later version.
          residue_indices = np.array(residue_indices) - 1

          #position_ids are the ids for positional encodings. 
          position_ids = torch.Tensor(residue_indices).long().reshape(1, -1)
          position_ids = position_ids.to(device) # the .to method only works with model, for a torch tensor it must be re-assigned

          # if parallel:

          #create tokenized version of the sequence for input in the model
          # inputs = self.tokenizer([sequence, sequence, sequence], return_tensors="pt", add_special_tokens=False)['input_ids']
          inputs = self.tokenizer([sequence], return_tensors="pt", add_special_tokens=False)['input_ids']
          inputs = inputs.to(device) # the .to method only works with model, for a torch tensor it must be re-assigned

          # print(f"{inputs.shape=}")
          # print(f"{inputs=}")

          # extend position_ids to the same length as the inputs
          # position_ids = position_ids.expand(inputs.shape[0], -1)

          import time
          t0 = time.time()
          #calculate output of the model given tokenized sequence and position_id for embedding
          raw_output = self.model(
            input_ids = inputs, 
            position_ids = position_ids
            ) 
          t1 = time.time()
          print( f"Iteration {self.i}: Time to run model: {t1-t0}")
          # Check the total GPU memory allocated in bytes
          allocated_memory_gib = torch.cuda.memory_allocated()/ (1024 ** 3)
          print(f"Iteration {self.i}: Total GPU memory allocated: {allocated_memory_gib:.2f} GiB")
          self.i += 1

          # HACK: the below code didn't work for me for some reason
          # for each attribute, if it's a pytorch tensor, move it to cpu
          for attr in dir(raw_output):
            if isinstance(getattr(raw_output, attr), torch.Tensor):
              setattr(raw_output, attr, getattr(raw_output, attr).to("cpu"))

          #Move the output to cpu. Not sure why it is needed to be honest
          # raw_output = tree_map(lambda x: x.to("cpu"), raw_output)

          pdb_string = output_to_pdb(raw_output)[0]

          if print_pdb:
            with open( pdb_file, "w+" ) as f:
                f.write( pdb_string )

          atoms: AtomArray = pdb_file_to_atomarray(StringIO(pdb_string))
          
          #print( f"AtomArray later used for calculations of energies" )
          #print( f"{atoms}" )

          plddt = raw_output["plddt"]
          plddt = plddt[0, ...].detach().numpy()
          plddt = plddt.transpose()
          
          #I also want to save the local plddt instead of averaging so I can use it later
          #for calculations if needed.
          local_plddt = plddt[atom_order["CA"], :]

          plddt = float(local_plddt.mean()) 

          ptm = float(raw_output["ptm"])
          pae = raw_output["predicted_aligned_error"]
          pae = pae.detach().numpy() #Corrected so it can be used by numpy later to calculate energy/calculate_PAE. 
          residue_index = raw_output["residue_index"]

          #print( f"Residue index from raw_output of esmfold: {residue_index}" )

          return FoldingResult( atoms=atoms, ptm=ptm, plddt=plddt, pae=pae, residue_index=residue_index, 
                                pdb=pdb_string, local_plddt = local_plddt )
