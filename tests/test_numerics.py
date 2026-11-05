import numpy as np
import pytest
from rds2023.statistics import Moments, directional_flux

@pytest.mark.parametrize("ddof", [0, 1])
def test_streaming_moments(ddof):
    x = np.random.default_rng(4).normal(size=(107, 4, 5, 3))+10
    m = Moments()
    for block in np.array_split(x, 7):
        m.update(block)
    np.testing.assert_allclose(m.mean, x.mean(0), atol=1e-13)
    centered = x-x.mean(0)
    expected = np.einsum("tzxi,tzxj->zxij", centered, centered)/(len(x)-ddof)
    np.testing.assert_allclose(m.covariance(ddof), expected, atol=1e-13)


def test_flux_sign_split_is_instantaneous():
    q = directional_flux(np.array([[2., -2.], [-2., 2.]]), np.ones(2), reference_speed=2, reference_measure=2)
    np.testing.assert_allclose(q, [0, .5, -.5])
    assert q[0] == pytest.approx(q[1]+q[2])
