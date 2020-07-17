from bionmr_utils.md.process.calc import calc_mean_square_displacement, calc_mean_square_displacement_by_axes
from bionmr_utils.md.process.calc import calc_autocorr, calc_inertia_tensor_vectors_autocorr
import numpy as np


def test_calc_mean_square_displacement():
    cm = np.array([[1, 1, 1],
                   [2, 2, 2],
                   [4, 6, 5],
                   [7, 11, 9],
                   ])
    time = np.linspace(0, 4, 4, endpoint=False)
    lag_index = [1, 2]
    time_lag, msds = calc_mean_square_displacement(time, cm, lag_index)

    msd_first = (sum((cm[1] - cm[0]) ** 2) + sum((cm[2] - cm[1]) ** 2) + sum(
        (cm[3] - cm[2]) ** 2)) / 3
    msd_second = (sum((cm[2] - cm[0]) ** 2) + sum((cm[3] - cm[1]) ** 2)) / 2
    np.testing.assert_allclose(msd_first, msds[0])
    np.testing.assert_allclose(msd_second, msds[1])


def test_calc_mean_square_displacement_by_axes():
    cm = np.array([[1, 1, 1],
                   [2, 2, 2],
                   [4, 6, 5],
                   [7, 11, 9],
                   ])
    time = np.linspace(0, 4, 4, endpoint=False)
    lag_index = [1, 2]
    time_lag, msds_x, msds_y, msds_z = calc_mean_square_displacement_by_axes(time, cm, lag_index)
    msd_first_x = ((cm[0][0] - cm[1][0]) ** 2 + (cm[2][0] - cm[1][0]) ** 2 + (cm[3][0] - cm[2][0]) ** 2) / 3
    msd_second_x = ((cm[2][0] - cm[0][0]) ** 2 + (cm[3][0] - cm[1][0]) ** 2) / 2
    np.testing.assert_allclose(msd_first_x, msds_x[0])
    np.testing.assert_allclose(msd_second_x, msds_x[1])


def test_calc_autocorr():
    vectors = {(1, "tests"): np.array([[1, 2, 3], [1.1, 2.1, 3.3], [0.9, 2.4, 3.5]])}
    acorr_order_2 = np.array([[1.0, 0.99566433, 0.99407925]])
    acorr_bionmr_utils = calc_autocorr(vectors)

    assert (1, "tests") == acorr_bionmr_utils.keys()
    np.testing.assert_allclose(acorr_order_2, list(acorr_bionmr_utils.values()))


def test_calc_inertia_tensor_vectors_autocorr():
    rotation_matrix = np.array([0.589, -0.263, 0.477, -0.452, -0.588, 0.688, -0.246, 0.915, -0.403, -0.659,
                                -0.193, 0.153, 0.099, -0.373, -0.756, -0.147, 0.513, -0.153, -0.011, 0.269,
                                -0.692, 0.649, 0.496, 0.061, 0.554, 0.638, -0.19])

    rotation_matrix = rotation_matrix.reshape(3, 3, 3)
    rotation_axes = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    rotation_axes_weights = np.array([0.8, 0.1, 0.1])

    acorr_order_2 = np.array([0.07957747, -0.01486589, -0.0245047])
    acorr_bionmr_utils = calc_inertia_tensor_vectors_autocorr(rotation_matrices=rotation_matrix,
                                                              rotation_axes=rotation_axes,
                                                              rotation_axes_weights=rotation_axes_weights)

    np.testing.assert_allclose(acorr_order_2, acorr_bionmr_utils, rtol=1e-7, atol=1e-5)
