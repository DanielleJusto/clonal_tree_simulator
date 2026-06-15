from cell import *
from chain import *
import random
import math
import pandas as pd
import pprint
import matplotlib.pyplot as plt
import copy

"""
TODO: 
- Cells divide at different rates based on cell type. Implement cell type. Division rate differs by cell type.
"""

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

def generate_10x_r1():
    cell_barcode = "".join(random.choices(["A", "C", "G", "T"], k=16))
    umi = "".join(random.choices(["A", "C", "G", "T"], k=10))
    return cell_barcode + umi

def generate_initial_population(light_chains):
    population = []
    barcodes = []
    
    # 1. Create heavy chain
    vH = 'CAGGTCCAACTGCAGCAGCCTGGGGCTGAGCTTGTGAAGCCTGGGGCTTCAGTGAAGCTGTCCTGCAAGGCTTCTGGCTACACCTTCACCAGCTACTGGATGCACTGGGTGAAGCAGAGCCTGGACGAGGCCTTGAGTGGATTGGAAGGATTGATCCTAATAGTGGTGGTACTAAGTACAATGAGAAGTTCAAGAGCAAGGCCACACTGACTGTAGACAAACCCTCCAGCACAGCCTACATGCAGCTCAGCAGCCTGACATCTGAGGACTCTGCGGTCTATTATTGTGCAAGA'
    dH = 'CAACTGGGAC' # DQ52-BALB/c : CAACTGGGA  DQ52-C57BL/6 : CTAACTGGGAC
    jH = 'ACTACTTTGACTACTGGGGCCAAGGCACCACTCTCACAGTCTCCTCA'
    heavy_chain_seq = ''.join(vH+dH+jH)
    heavy_chain = Chain(heavy_chain_seq, True)

    # 2. Generate unique barcodes for every cell
    while len(barcodes) < len(light_chains):
        barcode = generate_10x_r1()
        if not barcode in barcodes:
            barcodes.append(barcode)

    # 3. Combine all light chains with heavy chain
    for i in range(len(light_chains)):
        light_chain = light_chains[i]
        barcode = barcodes[i]
        cell = Cell(barcode=barcode, heavy_chain=heavy_chain, light_chain=light_chain, generation=0, root=barcode, parent=None, sampling_weight=1)
        population.append(cell)

    return population

def simulate(population, time_period, dP, division_rate, mutation_rate, death_rate):
        t = 0

        new_population = population.copy()

        while t < time_period and new_population:

            # 1. Mutate
            sampling_weights = [i.sampling_weight for i in new_population]
            sample_size = math.floor(mutation_rate * len(new_population))

            if sample_size > 0:
                sample = random.choices(new_population, weights=sampling_weights, k=sample_size)
                for target_cell in sample:
                    if target_cell in new_population:
                        mutated_cell = target_cell.mutate(dP)
                        # Replace 
                        new_population.remove(target_cell)
                        new_population.append(mutated_cell)
                
            # 2. Divide
            sample_size = math.floor(division_rate * len(new_population))
            if sample_size > 0:
                sample = random.sample(new_population, k=sample_size)
                for target_cell in sample:
                    if target_cell in new_population:
                        for i in range(2):
                            daughter_cell = Cell(
                                barcode=f"{target_cell.barcode}_{i}",
                                heavy_chain=Chain(target_cell.heavy_chain.sequence, target_cell.heavy_chain.is_heavy),
                                light_chain=Chain(target_cell.light_chain.sequence, target_cell.light_chain.is_heavy),
                                generation=target_cell.generation+1,
                                root=target_cell.root,
                                parent=target_cell.barcode,
                                sampling_weight=target_cell.sampling_weight
                            )
                            target_cell.children.append(daughter_cell.barcode)
                            new_population.append(daughter_cell)


            # 3. Cell death
            sample_size = math.floor(death_rate * len(new_population))
            if sample_size > 0:
                sample = random.sample(new_population, k=sample_size)
                for target_cell in sample:
                    if target_cell in new_population:
                        new_population.remove(target_cell)

            t += 1

        return new_population


def make_sparse(population, portion_removed):
    # 1. Randomly sample some portion of the cells
    sample = random.choices(population, k = math.floor(len(population)*portion_removed))
    # 2. Remove them from the population
    for i in sample:
        if i in population:
            population.remove(i)
    return population


def format_airr(population):
    pass

def format_csv(population):
    # 4. Create a dataframe of cells
    barcode = []
    heavy_chain = []
    light_chain = []
    generation = []
    root = []
    parent = []
    children = []
    sampling_weight = []

    for clone in population:
        barcode.append(clone.barcode)
        heavy_chain.append(clone.heavy_chain.sequence)
        light_chain.append(clone.light_chain.sequence)
        generation.append(clone.generation)
        root.append(clone.root)
        parent.append(clone.parent)
        children.append(clone.children)
        sampling_weight.append(clone.sampling_weight)

    data = {'barcode' : barcode,
            'heavy' : heavy_chain, # heavy chain sequence
            'light' : light_chain, # light chain sequence
            'gen' : generation,    # generation in lineage
            'root' : root,         # root node barcode
            'parent' : parent,
            'children' : children,
            'sampling_weight' : sampling_weight} # sampling weight
    
    df = pd.DataFrame.from_dict(data)

    return df


def main():

    # 1. Initialize population
    light_chain_options = generate_light_chains()
    population = generate_initial_population(light_chain_options)

    # 2. Simulate mutation and division
    simulation_length = 20              # total number of time steps in the simulation
    d_weight = 1                        # change in weight each mutation
    division_rate = 0.05                # portion of population that divides each iteration
    mutation_rate = 0.10                # portion of population that mutates each iteration
    death_rate = 0.02                   # poriton of population that dies each iteration

    new_population = simulate(population, simulation_length, d_weight, division_rate, mutation_rate, death_rate)

    # Save data to csv
    df = format_csv(new_population)
    df.to_csv('population.csv', index=False) # Save to csv

    # 5. Plot parent counts to verify that most of the clones are from the original population
    counts = df['root'].value_counts()

    counts.plot(kind='bar', color='skyblue')
    plt.xlabel('root')
    plt.ylabel('Count')
    plt.title('Size of each parent cell')
    plt.show()

       
if __name__ == "__main__":
    main()