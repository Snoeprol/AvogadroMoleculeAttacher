from itertools import chain, combinations
def powerset(iterable):
    s = list(iterable)
    print(s)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

list_atoms = []
for i in range(3):
    list_atoms.append(i)


print(list(powerset(list_atoms)))