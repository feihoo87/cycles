import operator
from functools import reduce
from itertools import islice, product

from ..cycles import Cycles, find_permutation
from ..group import PermutationGroup
from .chp import CZ, H, P
from .paulis import encode_paulis


def cliffordOrder(n: int) -> int:
    """
    Order of complex Clifford group of degree 2^n arising in quantum coding theory.
    
    Sloane, N. J. A. (ed.). "Sequence A003956 (Order of Clifford group)".
    The On-Line Encyclopedia of Integer Sequences. OEIS Foundation.
    https://oeis.org/A003956
    """
    return reduce(operator.mul, (((1 << (2 * j)) - 1) << 2 * j + 1
                                 for j in range(1, n + 1)), 1)


def make_clifford_generators(N: int,
                             graph: tuple[int, int] | None = None
                             ) -> dict[tuple, Cycles]:
    if graph is None:
        graph = []
        for i in range(N - 1):
            graph.append((i, i + 1))

    stablizers = []
    for s in islice(product('IXYZ', repeat=N), 1, None):
        n = encode_paulis(''.join(s))
        stablizers.append(n)
        stablizers.append(n | 2)

    generators = {}

    for i in range(N):
        generators[('H', i)] = find_permutation(stablizers,
                                                [H(s, i) for s in stablizers])
        generators[('S', i)] = find_permutation(stablizers,
                                                [P(s, i) for s in stablizers])

    for i, j in graph:
        generators[('CZ', i,
                    j)] = find_permutation(stablizers,
                                           [CZ(s, i, j) for s in stablizers])

    return generators


class CliffordGroup(PermutationGroup):

    def __init__(self, N: int, graph: tuple[int, int] | None = None):
        self.N = N
        generators = make_clifford_generators(N, graph)
        super().__init__(list(generators.values()))
        self.reversed_map = {v: k for k, v in generators.items()}
        self.gate_map = generators
        self.gate_map_inv = {v: k for k, v in generators.items()}

    def __len__(self):
        return cliffordOrder(self.N)
