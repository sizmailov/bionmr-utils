from typing import List, Iterator
from bionmr_utils.md import Frame


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

    residues = frame.residues

    assert sum(residues_per_chain) <= len(residues)

    if preserve_number_of_residues:
        assert sum(residues_per_chain) == len(residues)

    new_frame = Frame()
    for n_res, chain_name in zip(residues_per_chain, chain_names):
        new_molecule = new_frame.add_molecule()
        new_molecule.name = chain_name
        for r in residues[n:n + n_res]:
            new_residue = new_molecule.add_residue()
            new_residue.name = r.name
            new_residue.id = r.id
            for atom in r.atoms:
                new_atom = new_residue.add_atom()
                new_atom.id = atom.id
                new_atom.name = atom.name
                new_atom.mass = atom.mass
                new_atom.r = atom.r
                new_atom.vdw_radius = atom.vdw_radius

        n += n_res

    return new_frame
