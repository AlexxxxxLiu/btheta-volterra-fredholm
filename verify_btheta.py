#!/usr/bin/env python3
"""Numerical regression checks for the centered Volterra manuscript.

The calculations check normalization and finite-dimensional convergence. They
do not test the Riemann hypothesis and are not used in the analytic proof.
"""

from __future__ import annotations

import json

import mpmath as mp
import numpy as np


mp.mp.dps = 60


def xi(z: mp.mpc) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * z
    return (
        mp.mpf("0.5")
        * s
        * (s - 1)
        * mp.power(mp.pi, -s / 2)
        * mp.gamma(s / 2)
        * mp.zeta(s)
    )


def theta_density(x: mp.mpf, terms: int = 12) -> mp.mpf:
    x = abs(x)
    return 2 * sum(
        (
            2 * mp.pi**2 * n**4 * mp.exp(mp.mpf("4.5") * x)
            - 3 * mp.pi * n**2 * mp.exp(mp.mpf("2.5") * x)
        )
        * mp.exp(-mp.pi * n**2 * mp.exp(2 * x))
        for n in range(1, terms + 1)
    )


def fourier_check(z: mp.mpf) -> dict[str, str]:
    integral = 2 * mp.quad(
        lambda x: theta_density(x) * mp.cos(z * x),
        [0, mp.mpf("0.5"), 1, 2, 4],
    )
    target = xi(z)
    return {
        "z": mp.nstr(z, 8),
        "thetaFourier": mp.nstr(integral, 35),
        "xi": mp.nstr(target, 35),
        "absoluteError": mp.nstr(abs(integral - target), 8),
    }


def finite_determinant_check(size: int, z: float = 1.0) -> dict[str, float]:
    """Midpoint Nystrom test using p(x) proportional to exp(-x^4)."""

    radius = 3.2
    dx = 2 * radius / size
    x = -radius + (np.arange(size) + 0.5) * dx
    density = np.exp(-(x**4))
    density /= density.sum() * dx

    cdf = (np.cumsum(density) - 0.5 * density) * dx
    step = np.tril(np.ones((size, size)), -1) + 0.5 * np.eye(size)
    kernel = (
        1j
        * np.sqrt(density[None, :] / density[:, None])
        * (step - cdf[:, None])
    )
    matrix = dx * kernel

    sign, log_abs_det = np.linalg.slogdet(np.eye(size) + z * matrix)
    det2 = sign * np.exp(log_abs_det - z * np.trace(matrix))
    characteristic = np.sum(density * np.exp(1j * z * x)) * dx
    return {
        "gridSize": size,
        "absoluteError": float(abs(det2 - characteristic)),
        "det2Real": float(det2.real),
        "characteristicReal": float(characteristic.real),
    }


def main() -> None:
    r = mp.e ** (-3 * mp.pi)
    delta = 2 * ((1 + 11 * r + 11 * r**2 + r**3) / (1 - r) ** 5 - 1)
    d_star = (2 - 3 / mp.pi) * (1 - 5 / (4 * mp.pi))
    payload = {
        "precisionDigits": mp.mp.dps,
        "fourierChecks": [fourier_check(mp.mpf(v)) for v in (0, 1, 2)],
        "finiteDeterminantChecks": [
            finite_determinant_check(size) for size in (100, 200, 400, 800)
        ],
        "explicitConstants": {
            "Delta": mp.nstr(delta, 20),
            "dStar": mp.nstr(d_star, 20),
            "HilbertSchmidtSquaredUpperBound": mp.nstr(
                (1 + delta) / (mp.pi * d_star), 20
            ),
        },
        "interpretation": (
            "Regression checks only; the manuscript's proof is analytic."
        ),
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
