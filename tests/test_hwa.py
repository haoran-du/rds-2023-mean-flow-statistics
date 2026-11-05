import numpy as np
import pytest
from scipy.signal import welch
from scipy.signal.windows import hamming
from rds2023.hwa import welch_stream, corrected

def test_streaming_welch_matches_explicit_scipy():
    x = np.random.default_rng(12).normal(size=(1027, 15))
    channels = [6, 2, 0, 4, 8]
    mean = x[:, channels].mean(0)
    freq, p = welch_stream(x, "1b", channels, mean, nperseg=128, fs=15000)
    ref_freq, ref_p = welch(x[:, channels]-mean, fs=15000, window=hamming(128, sym=True),
                          nperseg=128, noverlap=64, nfft=128, detrend=False, axis=0)
    np.testing.assert_allclose(freq, ref_freq)
    np.testing.assert_allclose(p, ref_p, atol=1e-16)


def test_correction_never_changes_source():
    x = np.arange(300, dtype=float).reshape(20, 15)
    before = x.copy()
    result = corrected(x, "1b", 0, 4, [0, 2])
    result[:] = 0
    np.testing.assert_array_equal(x, before)


@pytest.mark.parametrize("policy", ["saved-zero", "recorded-replacement"])
@pytest.mark.parametrize("offsets", [(-2, 2), (8999, 9003), (-2, 9003)])
def test_b6_version_and_inclusive_interval(policy, offsets):
    class VirtualRecording:
        # Deterministic read-only data: no allocation of a 27-million-row file.
        def __getitem__(self, key):
            rows, columns = key
            return np.arange(rows.start, rows.stop, dtype=float)[:, None] + np.arange(15)[columns][None, :]/20

    data = VirtualRecording()
    lo, hi, source = 19375499, 19384500, 19385999
    start, stop = lo+offsets[0], lo+offsets[1]
    channels = [0, 6, 14]  # Camera trigger must receive the SAME intervention.
    actual = corrected(data, "6b", start, stop, channels, b6_policy=policy)
    expected = data[start:stop, :][:, channels]
    left, right = max(lo, start), min(hi, stop)
    expected[left-start:right-start] = (0 if policy == "saved-zero" else
        data[source+left-lo:source+right-lo, :][:, channels])
    np.testing.assert_array_equal(actual, expected)


def test_unknown_hwa_policy_rejected():
    with pytest.raises(ValueError, match="policy"):
        corrected(np.zeros((10, 15)), "6b", 0, 10, [0], b6_policy="guess")
