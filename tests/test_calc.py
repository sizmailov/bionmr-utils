from bionmr_utils.md.process.calc import calc_mean_square_displacement, calc_mean_square_displacement_by_axes
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
