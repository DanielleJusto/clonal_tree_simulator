import random
import copy
from chain import Chain
import uuid


class Cell:

    def __init__(self, barcode: str, heavy_chain: Chain, light_chain: Chain, generation: int, parent = None, sampling_weight=1):
       self.barcode = barcode
       self.clone_id =  uuid.uuid1() 
       self.heavy_chain = heavy_chain
       self.light_chain = light_chain
       self.heavy_mutations = 0
       self.light_mutations = 0
       self.generation = generation
       self.parent = parent
       self.children = []
       self.sampling_weight = sampling_weight

    """
    Apply a random point mutation to either the light or heavy chain
    TODO: Make point mutation only applied to V section for SHM
    """
    def mutate(self):
        target_chain = random.choice(['heavy_chain', 'light_chain'])
        
        if target_chain == 'heavy_chain':
            self.heavy_chain = self.heavy_chain.mutate()
            self.heavy_mutations += 1
        else:
            self.light_chain = self.light_chain.mutate()
            self.light_mutations += 1
        
    """
    Create up to four clones of a given cell.
    """
    def divide(self, max_clones):
        num_clones = random.randint(2, max_clones)
        clones = []

        for i in range(num_clones):
            clone = Cell(
                barcode=self.barcode,
                heavy_chain=copy.deepcopy(self.heavy_chain),
                light_chain=copy.deepcopy(self.light_chain),
                generation=self.generation + 1,
                parent=self,
                sampling_weight=self.sampling_weight + 1
            )
            clone.clone_id = f"{self.clone_id}_{i}"
            self.children.append(clone)
            clones.append(clone)

        return clones

    