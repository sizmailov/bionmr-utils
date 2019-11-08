from bionmr_utils.md.process.fit import fit_auto_correlation, decorated_fit_auto_correlation, fit_msd, fit_limit
import numpy as np


def test_fit_limit():
    data1 = np.array([8, 7, 6, 5, 4, 5, 6, 7, 8, 3, 2, -1])
    data2 = np.array([8, -7, 6, 5, 4, 5, 6, 7, 8, 3, 2, -1])

    assert fit_limit(data1, window_size=5) == 4
    assert fit_limit(data2, window_size=5) == 1


def test_fit_autocorrelation():
    x = np.linspace(1, 1000, 1000)
    y = 0.5 * np.exp(-x / 0.5) + 0.5 * np.exp(-x / 5)
    result = fit_auto_correlation(x, y, bounds=[[0, 0, 0, 1], [1, 1, 1, 10]])

    np.testing.assert_allclose(result, np.array([0.5, 0.5, 0.5, 5]), rtol=1e-03)


def test_decorated_fit_auto_correlation():
    x = np.linspace(1, 1000, 1000)
    y = 0.5 * np.exp(-x / 0.5) + 0.5 * np.exp(-x / 5)
    result = decorated_fit_auto_correlation(x, y, bounds=[[0, 0, 0, 1], [1, 0.5, 1, 6]])[1]

    np.testing.assert_allclose(result, np.array([0.5, 0.5, 0.5, 5]), rtol=1e-03)


def test_fit_msd():
    x = np.linspace(0, 1500, 500)
    y = 5 * x

    delta_x = np.random.random(len(x))
    x_delta = x
    x_delta[::10] += delta_x[::10]
    result = fit_msd(y, x_delta)

    np.testing.assert_allclose(result, 5, rtol=1e-02)
