"""Reusable statistics. Historical discretizations stay in paper adapters."""
import numpy as np


class Moments:
    """Mergeable stable mean/variance/covariance, with time on first axis.

    Missing/invalid samples must be handled by the caller using an explicit
    policy; this accumulator rejects NaNs instead of silently imputing them.
    """
    def __init__(self):
        self.n = 0
        self.mean = None
        self.m2 = None

    def update(self, batch):
        a = np.asarray(batch, dtype=float)
        if a.ndim < 2 or not len(a) or not np.isfinite(a).all():
            raise ValueError("Expected nonempty finite (time,...,component) batch")
        n = len(a)
        mean = a.mean(axis=0)
        centered = a - mean
        m2 = np.einsum("t...i,t...j->...ij", centered, centered)
        if self.n == 0:
            self.n, self.mean, self.m2 = n, mean, m2
            return
        if mean.shape != self.mean.shape:
            raise ValueError("Batch shape changed")
        delta = mean - self.mean
        total = self.n + n
        self.m2 += m2 + np.einsum("...i,...j->...ij", delta, delta) * self.n * n / total
        self.mean += delta * n / total
        self.n = total

    def covariance(self, ddof=0):
        if ddof not in (0, 1) or self.n <= ddof:
            raise ValueError("ddof must be 0 or 1 and less than sample count")
        return self.m2 / (self.n - ddof)


def directional_flux(velocity, weights, *, reference_speed, reference_measure, density=1.0):
    """Instantaneous sign split before temporal averaging; negative is SIGNED.

    velocity: (time,...), weights: matching spatial line/area weights in SI.
    Reference measure is nominal width or area, not automatically observed area.
    This distinction must be configured for the future idealized experiment.
    """
    a = np.asarray(velocity, float)
    w = np.asarray(weights, float)
    if a.ndim < 2 or a.shape[1:] != w.shape or len(a) == 0:
        raise ValueError("Velocity/weight shape mismatch")
    if not np.isfinite(a).all() or not np.isfinite(w).all() or (w < 0).any():
        raise ValueError("Finite velocities and nonnegative weights required")
    if reference_speed <= 0 or reference_measure <= 0 or density <= 0:
        raise ValueError("Positive reference scales required")
    scale = density / (reference_speed * reference_measure)
    return np.array([np.mean(np.sum(b * w, axis=tuple(range(1, a.ndim)))) * scale
                     for b in (a, np.maximum(a, 0), np.minimum(a, 0))])
