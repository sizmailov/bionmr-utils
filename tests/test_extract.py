from bionmr_utils.md import PdbFile
from bionmr_utils.md.process.select import get_mtsl_vectors, get_NH_vectors, get_methyl_vectors
from bionmr_utils.md.process.extract import extract_vectors
import numpy as np


def test_extract_NH_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    N_1 = np.array([-9.920, -8.301, 5.541])
    N_56 = np.array([9.382, 4.941, 1.349])

    H_1 = np.array([-9.909, -9.029, 6.241])
    H_56 = np.array([8.943, 5.351, 2.161])

    NH_1 = [N_1 - H_1]
    NH_56 = [N_56 - H_56]

    frame = PdbFile(path).get_frames()
    vectors = extract_vectors(frame, get_vectors=get_NH_vectors)

    assert (len(vectors)) == 55

    keys = sorted(list(vectors.keys()))
    np.testing.assert_allclose(vectors[keys[0]].to_numpy(), NH_1)
    np.testing.assert_allclose(vectors[keys[-1]].to_numpy(), NH_56)

    assert str(keys[0][0]) == '2'
    assert str(keys[-1][0]) == '56'


def test_extract_methyl_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    CD1_LEU2 = np.array([0.504, 0.021, 1.995])
    CG2_THR55 = np.array([10.884, 2.547, 3.218])

    HD11_LEU2 = np.array([-0.010, -0.523, 1.202])
    HG21_THR55 = np.array([10.525, 3.199, 4.014])

    CH3_LEU2 = [CD1_LEU2 - HD11_LEU2]
    CH3_THR55 = [CG2_THR55 - HG21_THR55]

    frame = PdbFile(path).get_frames()
    vectors = extract_vectors(frame, get_vectors=get_methyl_vectors)

    assert (len(vectors)) == 33

    keys = sorted(list(vectors.keys()), key=lambda x: x[0])

    np.testing.assert_allclose(vectors[keys[1]].to_numpy(), CH3_LEU2)
    np.testing.assert_allclose(vectors[keys[-1]].to_numpy(), CH3_THR55)

    assert str(keys[0][0]) == '2'
    assert str(keys[-1][0]) == '55'


def test_extract_mtsl_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    nitrogen = np.array([-7.476, 2.988, 1.402])
    oxygen = np.array([-8.401, 3.498, 0.733])

    electron = (nitrogen + oxygen) / 2

    H_1 = np.array([-9.909, -9.029, 6.241])
    H_56 = np.array([8.943, 5.351, 2.161])

    mtsl_H_1 = [H_1 - electron]
    mtsl_H_56 = [H_56 - electron]

    frame = PdbFile(path).get_frames()
    vectors = extract_vectors(frame, get_vectors=get_mtsl_vectors)

    assert (len(vectors)) == 55

    keys = sorted(list(vectors.keys()))
    np.testing.assert_allclose(vectors[keys[0]].to_numpy(), mtsl_H_1)
    np.testing.assert_allclose(vectors[keys[-1]].to_numpy(), mtsl_H_56)

    assert str(keys[0][0]) == '2'
    assert str(keys[-1][0]) == '56'
