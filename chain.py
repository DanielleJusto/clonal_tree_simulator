import random

class Chain:

        def __init__(self, v_sequence : str , d_sequence : str | None , j_sequence : str, is_heavy: bool):
            self.v_sequence = v_sequence
            self.d_sequence = d_sequence
            self.j_sequence = j_sequence
            self.is_heavy = is_heavy
            self.sequence = "".join(self.v_sequence + self.d_sequence + self.j_sequence) if self.is_heavy else "".join(v_sequence + self.j_sequence)


        """
        Apply a point mutation to the variable region of a chain's sequence
        """
        def mutate(self):
            seq = self.v_sequence

            # 1. Choose a random position to mutate
            random_ind = random.randint(0, len(seq) - 1)

            # 2. Choose mutation
            mutation = random.choices(list({"A", "C", "G", "T"} - {seq[random_ind]}))[0]

            # 3. Apply mutation
            mutated_seq = seq[:random_ind] + str(mutation) + seq[random_ind+1:]
            mutated_chain = Chain(mutated_seq, self.d_sequence, self.j_sequence, self.is_heavy)

            return mutated_chain