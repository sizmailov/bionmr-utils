from pyxmolpp2 import PdbFile
from bionmr_utils.md.process.select import get_mtsl_selection, get_NH_selection, get_methyl_selection
from bionmr_utils.md.process.extract import extract_vectors
import numpy as np


def test_extract_NH_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    N_1 = np.array([-9.920, -8.301, 5.541])
    N_56 = np.array([9.382, 4.941, 1.349])

    H_1 = np.array([-9.909, -9.029, 6.241])
    H_56 = np.array([8.943, 5.351, 2.161])

    NH_1 = N_1 - H_1
    NH_56 = N_56 - H_56

    frames = PdbFile(path).frames()
    rids_anames_pairs, vectors_generator = extract_vectors(frames, get_selection=get_NH_selection)
    vectors = list(vectors_generator)[0]

    assert (len(vectors)) == 55

    np.testing.assert_allclose(vectors[0], NH_1)
    np.testing.assert_allclose(vectors[-1], NH_56)

    assert rids_anames_pairs[0][0] == 2
    assert rids_anames_pairs[-1][0] == 56


def test_extract_methyl_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    CG2_THR2 = np.array([-8.902, -5.237, 7.570])
    CG2_THR55 = np.array([10.884, 2.547, 3.218])

    HG21_THR2 = np.array([-9.476, -4.531, 6.970])
    HG21_THR55 = np.array([10.525, 3.199, 4.014])

    CH3_THR2 = CG2_THR2 - HG21_THR2
    CH3_THR55 = CG2_THR55 - HG21_THR55

    frames = PdbFile(path).frames()
    rids_anames_pairs, vectors_generator = list(extract_vectors(frames, get_selection=get_methyl_selection))
    vectors = list(vectors_generator)[0]

    assert (len(vectors)) == 33

    np.testing.assert_allclose(vectors[0], CH3_THR2)
    np.testing.assert_allclose(vectors[-1], CH3_THR55)

    assert rids_anames_pairs[0][0] == 2
    assert rids_anames_pairs[-1][0] == 55


def test_extract_mtsl_vectors():
    path = "tests_dataset/amber/GB1_F30C_MTSL/box.pdb"
    nitrogen = np.array([-7.476, 2.988, 1.402])
    oxygen = np.array([-8.401, 3.498, 0.733])

    electron = (nitrogen + oxygen) / 2

    H_1 = np.array([-9.909, -9.029, 6.241])
    H_56 = np.array([8.943, 5.351, 2.161])

    mtsl_H_1 = electron - H_1
    mtsl_H_56 = electron - H_56

    frame = PdbFile(path).frames()
    rids_anames_pairs, vectors_generator = list(extract_vectors(frame, get_selection=get_mtsl_selection))
    vectors = list(vectors_generator)[0]

    assert (len(vectors)) == 55

    np.testing.assert_allclose(vectors[0], mtsl_H_1)
    np.testing.assert_allclose(vectors[-1], mtsl_H_56)

    assert rids_anames_pairs[0][0] == 2
    assert rids_anames_pairs[-1][0] == 56
