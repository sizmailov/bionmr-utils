from pyxmolpp2 import AtomSelection, Frame, rId, aName, rName
from typing import Tuple


def get_methyl_selection(frame: Frame) -> Tuple[AtomSelection, AtomSelection]:
    """

    :param frame: Frame
    :return: tuple of all methyl C-H atom pairs of given frame.
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

    C_atoms = []
    H_atoms = []

    for r in frame.residues:
        for atom in r.atoms:
            if (r.name, atom.name) in CH3_dict.keys():
                c_name, h_name = CH3_dict[r.name, atom.name]
                C_atoms.append(r[c_name])
                H_atoms.append(r[h_name])

    return AtomSelection(C_atoms), AtomSelection(H_atoms)


class UnpairedElectronSelection:
    def __init__(self, nitrogen_selection, oxygen_selection):
        self.nitrogen_selection = nitrogen_selection
        self.oxygen_selection = oxygen_selection

    class ElectronCoords:
        def __init__(self, nitrogen_coords, oxygen_coords):
            self.nitrogen_coords = nitrogen_coords
            self.oxygen_coords = oxygen_coords

        @property
        def values(self):
            return (self.nitrogen_coords.values + self.oxygen_coords.values) / 2

    @property
    def coords(self):
        return self.ElectronCoords(self.nitrogen_selection.coords, self.oxygen_selection.coords)

    @property
    def size(self):
        return self.nitrogen_selection.size

    @property
    def index(self):
        return self.nitrogen_selection.index


def get_NH_selection(frame: Frame) -> Tuple[AtomSelection, AtomSelection]:
    """

    :param frame: Frame
    :return: tuple of all backbone atom pairs of given frame.
    """
    return (frame.atoms.filter((rName != "PRO") & (aName == "N") & (rId > 1)),
            frame.atoms.filter((rName != "PRO") & (aName == "H") & (rId > 1)))


def get_mtsl_selection(frame: Frame,
                       label_name="CML",
                       ):
    """

    :param frame: Frame
    :return: tuple of all backbone atom pairs of given frame.
    """
    label_atoms = frame.atoms.filter(rName == label_name)
    electron = UnpairedElectronSelection(label_atoms.filter(aName == "N1"), label_atoms.filter(aName == "O1"))
    assert electron.size == 1

    atom_pairs_selection = (electron, frame.atoms.filter(aName == "H"))

    return atom_pairs_selection
