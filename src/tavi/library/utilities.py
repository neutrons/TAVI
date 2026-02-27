from scipy.constants import e, hbar, m_n
import numpy as np

def SE2K(e:float):
    """Convert energy E to momentum transfer k"""
    E2K=np.sqrt(2e-3 * e * m_n) * 1e-10 / hbar
    return np.sqrt(e)*E2K