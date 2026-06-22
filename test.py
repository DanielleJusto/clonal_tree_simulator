from simulate import *
from cell import Cell
from chain import Chain
from tree import Node
import pprint

"""
Testing integration of simulated data with Simble tree script
"""
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

def main():

    cell_1 = Cell('1111', Chain('TA', 'CG', 'AG', True), Chain('AC', None, 'GT', False), 0, None)
    cell_1.clone_id = "A"


    roots = [cell_1]
    population = copy.deepcopy(roots)

    # Step 1: Proliferate and mutate cell 1
    clones = cell_1.divide(max_clones=2)
    for clone in clones:
         clone.mutate()
    
    population.extend(clones)

    # Step 2: Proliferate and mutate cell one of cell 1's children
    clones = cell_1.children[0].divide(max_clones=2)
    for clone in clones:
         clone.mutate()

    population.extend(clones)

    # Step 3: Do that again to the new clones

    clones = cell_1.children[0].children[0].divide(max_clones=2)
    for clone in clones:
         clone.mutate()

    population.extend(clones)
    print(" ")
    pp = pprint.PrettyPrinter(indent=0, width=150) # Table pretty printer

    for i in roots:

        root = Node(i)
        newick = root.write_newick()
        
        pp.pprint(newick)
        print(" ")
        print(" ")
        
        tree = Phylo.read(StringIO(newick), "newick")
        Phylo.draw(tree)

main()