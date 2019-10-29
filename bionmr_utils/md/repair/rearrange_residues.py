from typing import List, Iterator
from bionmr_utils.md import (Frame, ChainName)


def rearrange_residues_in_chains(frame: Frame,
                                 residues_per_chain: List[int],
                                 chain_names: Iterator[str] = None,
                                 preserve_number_of_residues: bool = True
                                 ):
    """

    Returns new frame with residues rearranged into different chains

    :param frame: input frame
    :param residues_per_chain: number of residues in each chain
    :param chain_names: new chain names
    :param preserve_number_of_residues: enables check that initial and resulting frame have same number of residues
    :return: new frame with residues rearranged into different chains
    """
    from itertools import cycle

    if chain_names is None:
        chain_names = cycle("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

    assert isinstance(chain_names, str) and len(chain_names) == len(residues_per_chain)

    n = 0

    residues = frame.asResidues

    assert sum(residues_per_chain) <= len(residues)

    if preserve_number_of_residues:
        assert sum(residues_per_chain) == len(residues)

    new_frame = Frame(frame.index)
    for n_res, chain_name in zip(residues_per_chain, chain_names):
        chain = new_frame.emplace(ChainName(chain_name))
        for r in residues[n:n + n_res]:
            chain.emplace(r)
        n += n_res

    return new_frame
