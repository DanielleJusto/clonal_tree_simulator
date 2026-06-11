from chain import *

class Cell:

    def __init__(self, heavy_chain: Chain, light_chain: Chain, generation: int, parent = None, mutation_prob=1):
       self.heavy_chain = heavy_chain
       self.light_chain = light_chain
       self.generation = generation
       self.parent = parent
       self.mutation_probability = mutation_prob

    def mutate(self, d_prob):
        # 1. Select chain
        target_chain = random.choice(['heavy_chain', 'light_chain'])

        # 2. Mutate chain
        mutated_chain = getattr(self, target_chain).mutate()

        #3. Create new mutated cell
        updated_probability = self.mutation_probability + d_prob
        updated_generation = self.generation + 1
        if target_chain == 'heavy_chain':
            mutated_cell = Cell(mutated_chain, self.light_chain, updated_generation, self, self.mutation_probability)
        else:
            mutated_cell = Cell(self.heavy_chain, mutated_chain, updated_generation, self, self.mutation_probability)
        self.mutation_probability = updated_probability
        return mutated_cell
    