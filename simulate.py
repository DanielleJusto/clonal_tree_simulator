from cell import *
from chain import *
import random
import math
import pandas as pd
import matplotlib.pyplot as plt

def generate_light_chains():
    # 1. Define germline sequences 
    vL_sequences = [
    'NNNNNNNNNNTTACTCAGCCAAGCTCTGTGTCTACGTCTCTAGGAAGCACAGTCAAACTGCCTTGCAAGCGCAGCACTGGTAACATTGGAAACGACTATGTGCACTGGTACCAGCAACACATGGGAAGATCTCCCACCAATATGATCTATAGAGATGATCAGCGACCATCTGGAGTTTCTGATAGGTTCTCTGGCTCCATTGACAGCTCTTCCAACTCAGCCTTCCTGACAATCAATAATGTGCAGGCTGAGGATNNNNNNNNNNNNNNNNNN',
    'NNNNNNNNNNTTACTCAGCCAAGCTCTGTGTCTACGTCTCTAGGAAGCACAGTCAAACTGTCTTGCAAGCGCAGCACTGGTAACATTGGAAACAACTATGTGCACTGGTACCAGCAGTACATGGGAAGATCTCCCACCAATATGATCTATGATGATAATAAGCGACCATCTGGAGTTTCTGATAGGTTCTCTGGCTCCATTGACAGCTCTTCCAACTCAGCCTTCCTGACAATCAATAATGTGCAGGCTGAGGATNNNNNNNNNNNNNNNNNN',
    'NNNNNNNNNNTTACTCAGCCAAGCTCTGTGTCTACATCTCTAGGAAGCACAGTCAAACTGTCTTGCAAGCGCAGCACTGGTAACATTGGAAACAACTATGTGAACTGGTACCAGCAATACATGGGAAGATCTCCCACTAATATGATCTATGGAGATGATCAGCGACCAACTGGAGTTTCTGATAGGTTCTCTGGCTCCATTGACAGCTCTTCCAACTCAGCCTTCCTGACAATCAATAATGTGCAGGCTGAGGATNNNNNNNNNNNNNNNNNN',
    'NNNNNNNNNNTTACTCAGCCAAGCTCTGTGTCTACGTCTCTAGGAAGCACAGTCAAACTGTCTTGCAAGCCCAGCACTGGTAAAATTGGAAATTACTTTATGAGCTGGTACCAGCAACACATGGGAAGATCTCCCACAAATATGATCTATAGAGATGATCTCCGACCATCTGGAGTTTCTGATAGGTTCTCTGGCTCCATTGACAGCTCTTCCAACTCAGCCTTCCTGACAATCAATAATGTGCAGGCTGAGGATNNNNNNNNNNNNNNNNNN',
    'NNNNNNNNNNTTACTCAGCCAAGCTCTGTGTCTACATCTCTAGGAAGCACAGTCAAACTGCCTTGCAAGTGCAGCACTGGTAATATTGGAAGCTACTATGTGCATTGGTACCAGCAACATATGGGAAGATCTCCTACCAATATGATCCATAGTGATGATCAGCGACCATCTGGAGTTTCTGATAGGTTCTCTGGCTCCATTGACAGCTCTTCCAACTCAGCCTTCCTGACAATCAATAATGTGCAGGCTGAGGATNNNNNNNNNNNNNNNNNN',
    'CAGGCTGTTGTGACTCAGGAATCTGCACTCACCACATCACCTGGTGAAACAGTCACACTCACTTGTCGCTCAAGTACTGGGGCTGTTACAACTAGTAACTATGCCAACTGGGTCCAAGAAAAACCAGATCATTTATTCACTGGTCTAATAGGTGGTACCAACAACCGAGCTCCAGGTGTTCCTGCCAGATTCTCAGGCTCCCTGATTGGAGACAAGGCTGCCCTCACCATCACAGGGGCACAGACTGAGGATGAGGCAATATATTTCTGTGCTCTATGGTACAGCAACCATT',
    'CAGGCTGTTGTGACTCAGGAATCTGCACTCACCACATCACCTGGTGGAACAGTCATACTCACTTGTCGCTCAAGTACTGGGGCTGTTACAACTAGTAACTATGCCAACTGGGTTCAAGAAAAACCAGATCATTTATTCACTGGTCTAATAGGTGGTACCAGCAACCGAGCTCCAGGTGTTCCTGTCAGATTCTCAGGCTCCCTGATTGGAGACAAGGCTGCCCTCACCATCACAGGGGCACAGACTGAGGATGATGCAATGTATTTCTGTGCTCTATGGTACAGCACCCATT',
    'CAACTTGTGCTCACTCAGTCATCTTCAGCCTCTTTCTCCCTGGGAGCCTCAGCAAAACTCACGTGCACCTTGAGTAGTCAGCACAGTACGTACACCATTGAATGGTATCAGCAACAGCCACTCAAGCCTCCTAAGTATGTGATGGAGCTTAAGAAAGATGGAAGCCACAGCACAGGTGATGGGATTCCTGATCGCTTCTCTGGATCCAGCTCTGGTGCTGATCGCTACCTTAGCATTTCCAACATCCAGCCTGAAGATGAAGCAATATACATCTGTGGTGTGGGTGATACAATTAAGGAACAATTTGTG'
    ]

    jL_sequences = [
    'CTGGGTGTTCGGTGGAGGAACCAAACTGACTGTCCTAG',
    'TTATGTTTTCGGCGGTGGAACCAAGGTCACTGTCCTAG',
    'GTTTATTTTCGGCAGTGGAACCAAGGTCACTGTCCTAG'
    ]

    # 2. Generate combinations
    light_chains = []
    for v in vL_sequences:
        for j in jL_sequences:
            lc_sequence = ''.join(v+j)
            light_chain = Chain(lc_sequence, False)
            light_chains.append(light_chain)
    return light_chains

def generate_initial_population(light_chains):
    population = []
    
    # 1. Create heavy chain
    vH = 'CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA'
    dH = 'CAACTGGGAC' # DQ52-BALB/c : CAACTGGGA  DQ52-C57BL/6 : CTAACTGGGAC
    jH = 'ACTACTTTGACTACTGGGGCCAAGGCACCACTCTCACAGTCTCCTCA'
    heavy_chain_seq = ''.join(vH+dH+jH)
    heavy_chain = Chain(heavy_chain_seq, True)

    # 2. Combine all light chains with heavy chain
    for light_chain in light_chains:
        cell = Cell(heavy_chain, light_chain, 0)
        population.append(cell)

    return population

def simulate(population, time_period, dP):
        t = 0

        while t < time_period:

            # 1. Sample
            probabilities = [i.mutation_probability for i in population]
            sample = random.choices(population, probabilities, k=5) # number of cells sampled per iteration fixed

            # 2. Mutate
            for target_cell in sample:
                mutated_cell = target_cell.mutate(dP)
                # 3. Divide
                for i in range(0,2):
                    # 4. Update population
                    population.append(mutated_cell)
            t += 1

        return population


def make_sparse(population, portion_removed):
    # 1. Randomly sample some portion of the cells
    sample = random.choices(population, k = math.floor(len(population)*portion_removed))
    # 2. Remove them from the population
    for i in sample:
        if i in population:
            population.remove(i)
    return population

def main():

    # 1. Initialize population
    light_chain_options = generate_light_chains()
    population = generate_initial_population(light_chain_options)

    # 2. Simulate mutation and division
    simulation_length = 100
    d_prob = 1 # change in probability each mutation

    new_population = simulate(population, 100, d_prob)

    # # UNCOMMENT TO MAKE SPARSE DATA ------------------------------
    # # X3. Make data sparse -- randomly remove cells from population
    # new_population = make_sparse(new_population, 0.25)
    # # ------------------------------------------------------------

    # 4. Create a dataframe of cells
    heavy_chain = []
    light_chain = []
    generation = []
    mutation_prob = []
    for clone in new_population:
        heavy_chain.append(clone.heavy_chain.sequence)
        light_chain.append(clone.light_chain.sequence)
        generation.append(clone.generation)
        mutation_prob.append(clone.mutation_probability)

    data = {'heavy' : heavy_chain,
            'light' : light_chain,
            'gen' : generation, # generation in lineage
            'mutation_prob' : mutation_prob}
    
    df = pd.DataFrame.from_dict(data)

    # 5. Plot of clone generations to verify that the population is made mostly of clones
    #    Checking for power law distribution
    plt.hist(df['gen'])

    plt.xlabel("Generation")
    plt.ylabel("Frequency")

    plt.show()
       
if __name__ == "__main__":
    main()