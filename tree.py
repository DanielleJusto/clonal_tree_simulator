from chain import Chain
"""
TODO:
- For each parent cell:
    - Collect all cells with the same root
    - Find root
    - Find next immediate children, add to tree.
    - Find children of children, add to tree.
    - Continue recursively until reach tip node. 
"""

class Node:

    def __init__(self, cell):
       self.heavy_chain = cell.heavy_chain
       self.light_chain = cell.light_chain
       self.generation = cell.generation
       self.root = cell.root 
       self.sampling_weight = cell.sampling_weight


