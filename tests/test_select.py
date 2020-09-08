import numpy as np
from pyxmolpp2 import PdbFile, rName, aName
from bionmr_utils.md.process.select import get_NH_selection, get_methyl_selection
from bionmr_utils.md.process.select import UnpairedElectronSelection, get_mtsl_selection


def test_get_NH_selection():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    N_1 = np.array([-9.920, -8.301, 5.541])
    N_56 = np.array([9.382, 4.941, 1.349])

    H_1 = np.array([-9.909, -9.029, 6.241])
    H_56 = np.array([8.943, 5.351, 2.161])

    frame = PdbFile(path).frames()[0]
    NH_pairs_selection = get_NH_selection(frame)

    np.testing.assert_allclose(NH_pairs_selection[0].coords.values[0], N_1)
    np.testing.assert_allclose(NH_pairs_selection[0].coords.values[-1], N_56)

    np.testing.assert_allclose(NH_pairs_selection[1].coords.values[0], H_1)
    np.testing.assert_allclose(NH_pairs_selection[1].coords.values[-1], H_56)

    assert NH_pairs_selection[0].size == 55
    assert NH_pairs_selection[1].size == 55


def test_get_methyl_selection():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    CG2_THR2 = np.array([-8.902, -5.237, 7.570])
    CG2_THR55 = np.array([10.884, 2.547, 3.218])

    HG21_THR2 = np.array([-9.476, -4.531, 6.970])
    HG21_THR55 = np.array([10.525, 3.199, 4.014])

    frame = PdbFile(path).frames()[0]
    СН3_pairs_selection = get_methyl_selection(frame)

    np.testing.assert_allclose(СН3_pairs_selection[0].coords.values[0], CG2_THR2)
    np.testing.assert_allclose(СН3_pairs_selection[0].coords.values[-1], CG2_THR55)

    np.testing.assert_allclose(СН3_pairs_selection[1].coords.values[0], HG21_THR2)
    np.testing.assert_allclose(СН3_pairs_selection[1].coords.values[-1], HG21_THR55)

    assert СН3_pairs_selection[0].size == 33
    assert СН3_pairs_selection[1].size == 33


def test_UnpairedElectronSelection():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    nitrogen = np.array([-7.476, 2.988, 1.402])
    oxygen = np.array([-8.401, 3.498, 0.733])

    electron = (nitrogen + oxygen) / 2

    frame = PdbFile(path).frames()[0]
    nitrogen_selection = frame.atoms.filter((rName == "CML") & (aName == "N1"))
    oxygen_selection = frame.atoms.filter((rName == "CML") & (aName == "O1"))
    unpaired_electron = UnpairedElectronSelection(nitrogen_selection, oxygen_selection)

    np.testing.assert_allclose(unpaired_electron.coords.values[0], electron)


def test_get_mtsl_selection():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"

    H_1 = np.array([-9.909, -9.029, 6.241])
    H_56 = np.array([8.943, 5.351, 2.161])

    frame = PdbFile(path).frames()[0]
    mtsl_pairs_selection = get_mtsl_selection(frame)

    np.testing.assert_allclose(mtsl_pairs_selection[-1].coords.values[0], H_1)
    np.testing.assert_allclose(mtsl_pairs_selection[-1].coords.values[-1], H_56)

    assert mtsl_pairs_selection[0].size == 1
    assert mtsl_pairs_selection[1].size == 55


def test_proline_get_NH_selection():

    path = "tests_dataset/gromacs/xtc/1am7_protein.pdb"

    frame = PdbFile(path).frames()[0]
    N_selections, H_selections = get_NH_selection(frame)

    assert N_selections.size == H_selections.size
