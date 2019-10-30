import pyxmolpp2
from bionmr_utils.md import Frame, Atom, AtomName, rName
from typing import List, Tuple


def get_methyl_vectors(frame: Frame) -> List[Tuple[Atom, Atom]]:
    """

    :param frame: Frame
    :return: list of all methyl C-H atom pairs of given frame.
    """
    CH3_dict = {('ALA', 'CB'): ['CB', 'HB1'],
                ('VAL', 'CG1'): ['CG1', 'HG11'],
                ('VAL', 'CG2'): ['CG2', 'HG21'],
                ('THR', 'CG2'): ['CG2', 'HG21'],
                ('LEU', 'CD1'): ['CD1', 'HD11'],
                ('LEU', 'CD2'): ['CD2', 'HD21'],
                ('ILE', 'CD1'): ['CD1', 'HD11'],
                ('ILE', 'CG2'): ['CG2', 'HG21'],
                ('MET', 'CE1'): ['CE1', 'HE1']}

    atom_pairs = []

    for r in frame.asResidues:
        for atom in r.asAtoms:
            if (r.name.str, atom.name.str) in CH3_dict.keys():
                C = r[AtomName(CH3_dict[(r.name.str, atom.name.str)][0])]
                H = r[AtomName(CH3_dict[(r.name.str, atom.name.str)][1])]
                atom_pairs.append((C, H))

    return atom_pairs


def get_NH_vectors(frame: Frame) -> List[Tuple[Atom, Atom]]:
    """

    :param frame: Frame
    :return: list of all backbone atom pairs of given frame.
    """
    atom_pairs = []

    for r in frame.asResidues:
        try:
            atom_pairs.append((r[AtomName("N")], r[AtomName("H")]))
        except pyxmolpp2.polymer.OutOfRangeResidue:
            pass

    return atom_pairs


class UnpairedElectron:
    def __init__(self, nitrogen, oxygen):
        self.nitrogen = nitrogen
        self.oxygen = oxygen

    @property
    def rId(self):
        return self.nitrogen.rId

    @property
    def aName(self):
        return AtomName("e")

    @property
    def r(self):
        return (self.nitrogen.r + self.oxygen.r) / 2


def get_mtsl_vectors(frame: Frame) -> List[Tuple[UnpairedElectron, Atom]]:
    """

    :param frame: Frame
    :return: list of all backbone atom pairs of given frame.
    """
    r = frame.asResidues.filter(rName == "CML")[0]
    electron = UnpairedElectron(r[AtomName("N1")], r[AtomName("O1")])
    atom_pairs = []

    for r in frame.asResidues:
        try:
            atom_pairs.append((r[AtomName("H")], electron))
        except pyxmolpp2.polymer.OutOfRangeResidue:
            pass

    return atom_pairs
