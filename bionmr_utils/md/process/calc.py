import numpy as np
from typing import Tuple, Dict, List, Callable
from pyxmolpp2.geometry import calc_autocorr_order_2, VectorXYZ
from pyxmolpp2.polymer import ResidueId, ResidueName


def calc_autocorr(vectors: Dict[Tuple[ResidueId, ResidueName], VectorXYZ],
                  calc_autocorr_func: Callable[[VectorXYZ], List[float]] = calc_autocorr_order_2,
                  ) -> Dict[Tuple[ResidueId, ResidueName], List[float]]:
    """
    Get auto-correlation from trajectory

    :param vectors: dict of (rid, aname): VectorXYZ
    :calc_autocorr_func: function calculates auto-correlation
    :return dict of (rid, aname): auto-correlation
    """
    autocorr = {(rid, aname): calc_autocorr_func(vector)
                for (rid, aname), vector in vectors.items()
                }

    return autocorr


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
