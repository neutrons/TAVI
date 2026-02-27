import numpy as np
from scipy.constants import hbar, m_n


def SE2K(e: float):
    """Convert energy E to momentum transfer k"""
    E2K = np.sqrt(2e-3 * e * m_n) * 1e-10 / hbar
    return np.sqrt(e) * E2K
