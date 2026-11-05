from pathlib import Path
import h5py
import numpy as np
import pytest
from scipy.io import savemat
from rds2023.io import dantec, one_mat, mat_blocks, premultiplied_spectrum


@pytest.mark.parametrize("frequency_column", [0, 1])
def test_spectrum_column_order_is_explicit(frequency_column):
    pair = np.array([[0., 4.], [1., 2.], [2., 1.]])
    if frequency_column == 1:
        pair = pair[:, ::-1]
    f, spectrum = premultiplied_spectrum(pair, 2., frequency_column=frequency_column)
    np.testing.assert_array_equal(f, [1., 2.])
    np.testing.assert_array_equal(spectrum, [1., 1.])
    with pytest.raises(ValueError, match="ordering"):
        premultiplied_spectrum(pair, 2., frequency_column=1-frequency_column)


def test_mat_column_major_roundtrip(tmp_path):
    a = np.arange(24).reshape(2, 3, 4)
    p = tmp_path/"v73.mat"
    with h5py.File(p, "w") as f:
        f.create_dataset("u_contour", data=a.transpose(2, 1, 0))
    np.testing.assert_array_equal(one_mat(p), a)
    np.testing.assert_array_equal(np.concatenate(list(mat_blocks(p, "u_contour", "XZ", 1))), a)


def test_one_mat_rejects_ambiguous_variables(tmp_path):
    p = tmp_path/"v5.mat"
    savemat(p, {"a": np.ones(2), "b": np.zeros(3)})
    with pytest.raises(ValueError):
        one_mat(p)


def test_dantec_header_and_grid(tmp_path):
    p = tmp_path/"sample.txt"
    p.write_text(">>*HEADER*<<\nGridSize:{Width=2, Height=1}\n>>*DATA*<<\nx\ty\n1\t2\n3\t4\n")
    meta, a = dantec(p)
    assert meta["nx"] == 2 and meta["nz"] == 1
    np.testing.assert_array_equal(a, [[1, 2], [3, 4]])
