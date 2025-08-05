import time
from functools import lru_cache

import numpy as np


# --- Cached Helper Function ---
# This function works on individual (scalar) Q-points, which are hashable.
# Its results will be cached. We set maxsize=None for an unbounded cache.
# If memory is a concern, you can set it to a specific number (e.g., maxsize=4096).
@lru_cache(maxsize=None)
def _calculate_single_point(q1, q2, q3):
    """
    Helper function to calculate dispersion for a single Q-point.
    The results of this function are cached.
    DO NOT CALL DIRECTLY. Use model_disp_cached instead.
    """
    # To see the cache in action, you can uncomment the following line:
    # print(f"Cache miss! Calculating for point: ({q1:.2f}, {q2:.2f}, {q3:.2f})")

    sj = 5
    gamma_q = (np.cos(2 * np.pi * q1) + np.cos(2 * np.pi * q2) + np.cos(2 * np.pi * q3)) / 3
    disp = 2 * sj * (1 - gamma_q)
    time.sleep(0.2)
    # Return a tuple of the two energy bands for this point
    return (disp - 2, disp + 2)


# --- Main Function with Caching ---
def model_disp_cached(vq1, vq2, vq3):
    """
    Return energy for given Q points, with caching for individual points.
    3d FM J=-1 meV S=1, formula: E = 6*S*J*(1-cos(Q))

    Args:
        vq1 (array-like): Vector of Qx coordinates.
        vq2 (array-like): Vector of Qy coordinates.
        vq3 (array-like): Vector of Qz coordinates.

    Returns:
        np.ndarray: An array of shape (2, N) where N is the number of Q-points.
                    The first row is the lower energy band, the second is the upper.
    """
    # Ensure inputs are NumPy arrays for consistent processing
    vq1 = np.asarray(vq1)

    # Handle the case where a single point (scalar) is passed as input
    if vq1.ndim == 0:
        result = _calculate_single_point(vq1.item(), np.asarray(vq2).item(), np.asarray(vq3).item())
        # Reshape to the standard (2 bands, 1 point) output format
        return np.array(result).reshape(2, 1)

    # If we have vectors, iterate over each point, call the cached helper,
    # and collect the results. The list comprehension is a concise way to do this.
    results = [_calculate_single_point(q1, q2, q3) for q1, q2, q3 in zip(vq1, vq2, vq3)]

    # The `results` is a list of tuples: [(e1_p1, e2_p1), (e1_p2, e2_p2), ...]
    # Convert this list to a NumPy array and transpose it to get the
    # desired (2, N) shape, where each column is a point.
    disp = np.array(results).T

    return disp


# --- Original Function (for comparison) ---
def model_disp_original(vq1, vq2, vq3):
    """return energy for given Q points
    3d FM J=-1 meV S=1, en=6*S*J*(1-cos(Q))
    """
    sj = 5
    gamma_q = (np.cos(2 * np.pi * vq1) + np.cos(2 * np.pi * vq2) + np.cos(2 * np.pi * vq3)) / 3
    disp = 2 * sj * (1 - gamma_q)
    disp = np.array((disp - 2, disp + 2))
    num_disp = len(disp.shape)
    if num_disp == 1:
        disp = np.reshape(disp, (1, np.size(disp)))
    time.sleep(0.2)
    return disp


# --- Verification Example ---
if __name__ == "__main__":
    # Define some Q points, notice the duplicates: (0,0,0) and (0.5,0,0)
    q1_vec = np.array([0.0, 0.5, 0.0, 0.25, 0.5] * 10)
    q2_vec = np.array([0.0, 0.0, 0.0, 0.25, 0.0] * 10)
    q3_vec = np.array([0.0, 0.0, 0.0, 0.25, 0.0] * 10)
    q1_prime_vec = np.array([0.1, 0.5, 0.0, 0.25, 0.5] * 10)

    print("--- Original function ---")
    # On the first run, all unique points will be calculated and their results stored.
    # If you uncomment the `print` in _calculate_single_point, you will see "Cache miss!" here.
    start_time = time.perf_counter()
    result1 = model_disp_original(q1_vec, q2_vec, q3_vec)
    end_time = time.perf_counter()
    print(f"Original function took: {end_time - start_time:.6f} seconds")

    print("--- First call with cached function ---")
    # On the first run, all unique points will be calculated and their results stored.
    # If you uncomment the `print` in _calculate_single_point, you will see "Cache miss!" here.
    start_time = time.perf_counter()
    result1 = model_disp_cached(q1_vec, q2_vec, q3_vec)
    end_time = time.perf_counter()
    print(f"First call took: {end_time - start_time:.6f} seconds")

    print("\n--- Second call with the same data ---")
    # On the second run, all results are retrieved from the cache instantly.
    # No "Cache miss!" messages will be printed.
    start_time = time.perf_counter()
    result2 = model_disp_cached(q1_prime_vec, q2_vec, q3_vec)
    end_time = time.perf_counter()
    print(f"Second (cached) call took: {end_time - start_time:.6f} seconds")

    print("\n--- Cache Information ---")
    # You can inspect the cache's performance
    # hits: times a result was retrieved from cache
    # misses: times a result had to be computed
    # currsize: number of items currently in the cache
    print(_calculate_single_point.cache_info())
    # Expected output: CacheInfo(hits=2, misses=3, maxsize=None, currsize=3)
    # Hits=2 because (0,0,0) and (0.5,0,0) were repeated.
    # Misses=3 for the unique points (0,0,0), (0.5,0,0), and (0.25,0.25,0.25).

    # You can also clear the cache if needed
    _calculate_single_point.cache_clear()
    print("\nCache cleared!")
    print(_calculate_single_point.cache_info())
    # Expected output: CacheInfo(hits=0, misses=0, maxsize=None, currsize=0)
