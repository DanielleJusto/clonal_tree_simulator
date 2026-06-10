from chain import *

class Cell:

    def __init__(self, heavy_chain: Chain, light_chain: Chain, generation: int, mutation_prob=0.25):
       self.heavy_chain = heavy_chain
       self.light_chain = light_chain
       self.generation = generation,
       self.mutation_probability = mutation_prob

    def mutate(self, time_step, d_prob):
        # 1. Select chain
        # target_chain = random.choice(['heavy_chain', 'light_chain'])
        target_chain = 'heavy_chain'

        # 2. Mutate chain
        mutated_chain = getattr(self, target_chain).mutate()

        #3. Create new mutated cell
        updated_probability = self.mutation_probability + d_prob
        if target_chain == 'heavy_chain':
            mutated_cell = Cell(mutated_chain, self.light_chain, time_step, updated_probability)
        else:
            mutated_cell = Cell(self.heavy_chain, mutated_chain, time_step, updated_probability)

        return mutated_cell
    