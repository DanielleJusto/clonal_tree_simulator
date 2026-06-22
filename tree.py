from cell import Cell
import pandas as pd 


"""
Script to convert simulated data to Newick Tree format. Code adapted from Simble. 
"""
class Node:

    def __init__(self, cell, parent = None):
        self.cell = cell
        self.parent = parent
        self.children = cell.children
        self.heavy_mutations = cell.heavy_mutations
        self.light_mutations = cell.light_mutations
        self.generation = cell.generation
        

    def write_newick(self, time_tree=False):
        """Writes the node and its children in Newick format.

        Args:
            time_tree (bool): Whether to write the tree with time information.
        Returns:
            str: The Newick representation of the node and its children.
        """
        newick = self._write_newick_iteratively(time_tree=time_tree)
        
        # If the root node has a parent, include in the newick tree
        if self.parent is not None:
            parent_node = Node(self.cell.parent, parent=None)
            parent_newick = parent_node.write_newick_node(time_tree=time_tree, subtrees=[newick])
            return parent_newick
        
        return newick
    
    def write_newick_node(self, time_tree=False, subtrees=None):
        """Writes the node in Newick format.

        Args:
            time_tree (bool): Whether to write the tree with time information.
            subtrees (list): A list of Newick strings for the children.
        Returns:
            str: The Newick representation of the node and, 
                if subtrees' Newick strings are provided, its children.
        """
        name = f"{self.cell.clone_id}"
        labels = (
            f"cell_root={self.cell.barcode}",
            f"generation={self.generation}",
            f"heavy_chain={self.cell.heavy_chain.sequence}",
            f"light_chain={self.cell.light_chain.sequence}"
        )
        if time_tree:
            branch_length = str(self.time_since_last_split)
        else:
            branch_length = str(self.heavy_mutations+self.light_mutations)
        branch = f':{branch_length}'
        labels = f"[&{labels}]"
        if len(self.children)==0:
            children = ""
        elif subtrees is None or len(subtrees) == 0:
            children = ""
        else:
            children = "(" + ",".join(subtrees) + ")"
        return children + name + labels + branch


    def _write_newick_iteratively(self, time_tree=False):
        """Writes the tree in Newick format iteratively.

        Args:
            tree (Node): The root node of the tree/subtree.
            time_tree (bool): Whether to write the tree with time information.
        Returns:
            str: The Newick representation of the tree.
        """
        stack = [self]
        children_newick = {}
        newick = ""

        def add_to_newick_dict(node, newick):
            # for memory efficiency, once we've added this node's newick to its parent
            # we can remove it from the dict
            if node.parent is None:
                return newick
            if node.parent not in children_newick:
                children_newick[node.parent] = []
            children_newick[node.parent].append(newick)
            if node in children_newick:
                children_newick.pop(node)
            return ""

        while len(stack) > 0:
            current = stack.pop()
            number_of_children = len(current.children)
            child_newicks = children_newick.get(current, [])
            number_of_child_newicks = len(child_newicks)
            if len(current.children) == 0:
                # leaf node
                # this should be handled by number_of_children == number_of_child_newicks
                curr_newick = current.write_newick_node(time_tree=time_tree)
                newick += add_to_newick_dict(current, curr_newick)
            elif number_of_children == number_of_child_newicks:
                # all children have been processed
                curr_newick = current.write_newick_node(time_tree=time_tree, subtrees=child_newicks)
                # add the newick string to the parent and if there is no parent
                # i.e. we have the root, then we can just write the newick string
                newick += add_to_newick_dict(current, curr_newick)
            else:
                # not all children have been processed
                # this node can't be processed yet, so push it back onto the stack
                stack.append(current)
                # then push all children onto the stack so the children are above current node
                for child in current.children:
                    stack.append(Node(child, parent=current))
        return newick
 
