import random

class Chain:

        def __init__(self, sequence: str, is_heavy: bool):
            self.sequence = sequence
            self.is_heavy = is_heavy

        def mutate(self):
            seq = self.sequence

            # 1. Choose a random position to mutate
            # random_ind = random.randint(0,len(seq)-1)
            random_ind = 0

            # 2. Choose mutation
            mutation = random.choices(list({"A", "C", "G", "T"} - {seq[random_ind]}))[0]

            # 3. Apply mutation
            mutated_seq = seq[:random_ind] + str(mutation) + seq[random_ind+1:]
            mutated_chain = Chain(mutated_seq, self.is_heavy)

            return mutated_chain