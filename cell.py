import random
import copy
from chain import Chain

"""
TODO: Add parent and children for tree construction. 
"""

class Cell:

    def __init__(self, barcode: str, heavy_chain: Chain, light_chain: Chain, generation: int, root : str, parent, sampling_weight=1):
       self.barcode = barcode
       self.heavy_chain = heavy_chain
       self.light_chain = light_chain
       self.generation = generation
       self.root = root 
       self.parent = parent
       self.children = []
       self.sampling_weight = sampling_weight

    def mutate(self, d_prob):
        # 1. Select chain
        target_chain = random.choice(['heavy_chain', 'light_chain'])

        # 2. Mutate chain
        mutated_chain = getattr(self, target_chain).mutate()

        #3. Create new mutated cell
        updated_weight = self.sampling_weight + d_prob
        mutated_cell = Cell(self.barcode, 
                            mutated_chain if target_chain == 'heavy_chain' else copy.deepcopy(self.heavy_chain),
                            copy.deepcopy(self.light_chain) if target_chain == 'heavy_chain' else mutated_chain,
                            self.generation,
                            self.root,
                            self.parent,
                            updated_weight)
        return mutated_cell
    