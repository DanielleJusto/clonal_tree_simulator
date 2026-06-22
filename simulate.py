from cell import *
from chain import *
import random
import math
import pandas as pd
from Bio import Phylo
from io import StringIO
import pprint
import matplotlib.pyplot as plt
from tree import Node

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
            light_chain = Chain(v, None, j, False)
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
    heavy_chain = Chain(vH, dH, jH, True)

    # 2. Generate unique barcodes for every cell
    while len(barcodes) < len(light_chains):
        barcode = generate_10x_r1()
        if not barcode in barcodes:
            barcodes.append(barcode)

    # 3. Combine all light chains with heavy chain
    for i in range(len(light_chains)):
        light_chain = light_chains[i]
        barcode = barcodes[i]
        cell = Cell(barcode=barcode, 
                    heavy_chain=heavy_chain, 
                    light_chain=light_chain, 
                    generation=0, 
                    parent=None, 
                    sampling_weight=1)
        cell.root = cell
        population.append(cell)

    return population

def simulate(population, time_period, division_rate, death_rate, max_clones):
        t = 0

        new_population = population.copy()

        while t < time_period and new_population:
            
            # 1. Proliferate
            sampling_weights = [cell.sampling_weight for cell in new_population]
            sample_size = math.floor(division_rate * len(new_population))

            if sample_size > 0:
                sample = random.choices(new_population, weights=sampling_weights, k=sample_size)
                for target_cell in sample:
                    if target_cell in new_population:
                        clones = target_cell.divide(max_clones)
            
            # 2. Mutate
                for clone in clones:
                     clone.mutate()
            
            # 3. Add to population
            new_population.extend(clones)

            t += 1

        return new_population

def format_population_data(population):
    barcode = []
    clone_id = []
    heavy_chain = []
    light_chain = []
    heavy_mutations = []
    light_mutations = []
    generation = []

    for cell in population:
        barcode.append(cell.barcode)
        clone_id.append(cell.clone_id)
        heavy_chain.append(cell.heavy_chain.sequence)
        light_chain.append(cell.light_chain.sequence)
        heavy_mutations.append(cell.heavy_mutations)
        light_mutations.append(cell.light_mutations)
        generation.append(cell.generation)

    table = {
        'barcode' : barcode,
        'clone_id' : clone_id,
        'heavy_chain' : heavy_chain,
        'light_chain' : light_chain,
        'heavy_mutations' : heavy_mutations,
        'light_mutations' : light_mutations,
        'generation' : generation
    }

    df = pd.DataFrame(table)
    df.to_csv('data/population.csv', index=False)

    return table

def format_training_data(old_population, new_population):
    barcodes = [cell.barcode for cell in old_population]
    lineage_table = {key: [] for key in barcodes}

    for cell in new_population:
        lineage_table[cell.barcode].append(cell.clone_id)

    lineages = list(lineage_table.values())

    root_nodes = [cell for cell in old_population]
    newick_trees = []
    for root in root_nodes:
        root_as_node = Node(root)
        newick_trees.append(root_as_node.write_newick())

    table = {
        'lineage' : lineages,
        'tree' : newick_trees
    }

    df = pd.DataFrame(table)
    df.to_csv('data/trees.csv', index=False)

    return table
    

def main():

    # 1. Initialize population
    light_chain_options = generate_light_chains()
    population = generate_initial_population(light_chain_options)

    # 2. Simulate mutation and division
    simulation_length = 10             # total number of time steps in the simulation
    division_rate = 0.25               # portion of population that divides each iteration
    death_rate = 0.02                  # poriton of population that dies each iteration
    max_clones = 4                     # maximum number of clones that result from proliferation step

    new_population = simulate(population, simulation_length, division_rate, death_rate, max_clones)

    # 3. Format data into csvs
    population_dict = format_population_data(new_population)
    tree_dict = format_training_data(population, new_population)

    trees = tree_dict['tree']
 
    for newick in trees:
        tree = Phylo.read(StringIO(newick), "newick")
        Phylo.draw(tree)

       
if __name__ == "__main__":
    main()