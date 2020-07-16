import numpy as np
from typing import Tuple, Dict, List, Callable
from pyxmolpp2 import calc_autocorr_order_2
from tqdm import tqdm


def calc_autocorr(vectors: Dict[Tuple[str, str], np.ndarray],
                  calc_autocorr_func: Callable[[np.ndarray], List[float]] = calc_autocorr_order_2,
                  limit=-1
                  ) -> Dict[Tuple[str, str], List[float]]:
    """
    Get auto-correlation from trajectory

    :param vectors: dict of (rid, aname): VectorXYZ
    :calc_autocorr_func: function calculates auto-correlation
    :limit: length of auto-correlation function
    :return dict of (rid, aname): auto-correlation
    """
    autocorr = {(rid, aname): calc_autocorr_func(vector, limit)
                for (rid, aname), vector in vectors.items()
                }
    return autocorr


def calc_inertia_tensor_vectors_autocorr(rotation_matrices: np.array,
                                         rotation_axes: np.array,
                                         rotation_axes_weights: np.array,
                                         limit: int = -1,
                                         ) -> np.array:
    """
    Get auto-correlation for inertia vectors

    :rotation matrices: rotation matrices as numpy array of shape (N, 3, 3), where N is length of trajectory
    :rotation_axes: (x, y, z) coordinates of axes as numpy array of shape (N, 3), where N is number of axes
                    corresponding to the directions in space along which rotation is considered
    :rotation_axes_weights: significance of the rotation axes as numpy array of shape (N, 1)
    :limit: length of auto-correlation function
    :return: auto-correlation function
    """

    number_of_vectors = rotation_matrices.shape[0]
    sum_acorr = np.zeros(number_of_vectors)

    for r, w in tqdm(zip(rotation_axes, rotation_axes_weights), desc="calc autocorr"):
        vectors = r.reshape((3, 1)).dot(rotation_matrices)
        autocorr = np.array(calc_autocorr_order_2(vectors, limit=limit))
        sum_acorr += autocorr * w

    return np.array(sum_acorr / 4 / np.pi)


def calc_mean_square_displacement(time: np.array,
                                  mass_centers: np.array,
                                  lag_index: np.array
                                  ) -> Tuple[np.array, np.array]:
    """
    :param time: time-points between mass_centers
    :param mass_centers: N*3 array [x,y,z] columns
    :param lag_index: int array of time lags between mass_centers
    :return: tuple of two arrays (time_lag, msd)
    """

    assert time.shape[0] == mass_centers.shape[0]
    assert mass_centers.shape[1] == 3

    time_lags = []
    msds = []
    for int_lag in lag_index:
        lag_time = time[int_lag]
        msd = ((mass_centers[int_lag:] - mass_centers[:-int_lag]) ** 2).sum(axis=1).mean()
        time_lags.append(lag_time)
        msds.append(msd)

    return np.array(time_lags), np.array(msds)


def calc_mean_square_displacement_by_axes(time: np.array,
                                          mass_centers: np.array,
                                          lag_index: np.array
                                          ) -> Tuple[np.array, np.array, np.array, np.array]:
    """
    :param time: time-points between mass_centers
    :param mass_centers: N*3 array [x,y,z] columns
    :param lag_index: int array of time lags between mass_centers
    :return: tuple of four arrays (time_lag, msd_x, msd_y, msd_z)
    """

    assert mass_centers.shape[1] == 3

    mass_centers_x = mass_centers[:, 0]
    mass_centers_y = mass_centers[:, 1]
    mass_centers_z = mass_centers[:, 2]

    time_lags = []
    msds_x = []
    msds_y = []
    msds_z = []

    for int_lag in lag_index:
        lag_time = time[int_lag]
        msd_x = ((mass_centers_x[int_lag:] - mass_centers_x[:-int_lag]) ** 2).mean()
        msd_y = ((mass_centers_y[int_lag:] - mass_centers_y[:-int_lag]) ** 2).mean()
        msd_z = ((mass_centers_z[int_lag:] - mass_centers_z[:-int_lag]) ** 2).mean()
        time_lags.append(lag_time)
        msds_x.append(msd_x)
        msds_y.append(msd_y)
        msds_z.append(msd_z)

    return np.array(time_lags), np.array(msds_x), np.array(msds_y), np.array(msds_z)
