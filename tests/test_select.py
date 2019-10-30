from bionmr_utils.md import PdbFile, rName, AtomName
from bionmr_utils.md.process.select import get_NH_vectors, get_methyl_vectors, UnpairedElectron, get_mtsl_vectors
import numpy as np


def test_get_NH_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    N_1 = [-9.920, -8.301, 5.541]
    N_56 = [9.382, 4.941, 1.349]

    H_1 = [-9.909, -9.029, 6.241]
    H_56 = [8.943, 5.351, 2.161]

    frame = PdbFile(path).get_frame()
    NH_pairs = get_NH_vectors(frame)

    np.testing.assert_allclose(NH_pairs[0][0].r.to_np, N_1)
    np.testing.assert_allclose(NH_pairs[-1][0].r.to_np, N_56)

    np.testing.assert_allclose(NH_pairs[0][1].r.to_np, H_1)
    np.testing.assert_allclose(NH_pairs[-1][1].r.to_np, H_56)

    assert len(NH_pairs) == 55


def test_get_methyl_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    CD1_LEU5 = np.array([0.504, 0.021, 1.995])
    CG2_THR55 = [10.884, 2.547, 3.218]

    HD11_LEU5 = [-0.010, -0.523, 1.202]
    HG21_THR55 = [10.525, 3.199, 4.014]

    frame = PdbFile(path).get_frame()
    СН3_pairs = get_methyl_vectors(frame)

    np.testing.assert_allclose(СН3_pairs[1][0].r.to_np, CD1_LEU5)
    np.testing.assert_allclose(СН3_pairs[-1][0].r.to_np, CG2_THR55)

    np.testing.assert_allclose(СН3_pairs[1][1].r.to_np, HD11_LEU5)
    np.testing.assert_allclose(СН3_pairs[-1][1].r.to_np, HG21_THR55)

    assert len(СН3_pairs) == 33


def test_UnpairedElectron():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    nitrogen = np.array([-7.476, 2.988, 1.402])
    oxygen = np.array([-8.401, 3.498, 0.733])

    electron = (nitrogen + oxygen) / 2

    frame = PdbFile(path).get_frame()
    r = frame.asResidues.filter(rName == "CML")[0]
    unpaired_electron = UnpairedElectron(r[AtomName("N1")], r[AtomName("O1")])

    np.testing.assert_allclose(unpaired_electron.r.to_np, electron)


def test_get_mtsl_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"

    H_1 = [-9.909, -9.029, 6.241]
    H_56 = [8.943, 5.351, 2.161]

    frame = PdbFile(path).get_frame()
    mtsl_pairs = get_mtsl_vectors(frame)

    np.testing.assert_allclose(mtsl_pairs[0][1].r.to_np, H_1)
    np.testing.assert_allclose(mtsl_pairs[-1][1].r.to_np, H_56)

    assert len(mtsl_pairs) == 55
