import numpy as np
from tavi.sample.sample import Sample


def test_b_matrix(a, b, c, alpha, beta, gamma, b_matrix):

    sample = Sample(lattice_params=(a, b, c, alpha, beta, gamma))
    print(sample.b_mat())
    print(np.allclose(sample.b_mat(), b_matrix, atol=1e-4))


def test_ub_matrix_to_uv(a, b, c, alpha, beta, gamma, ub_matrix):

    sample = Sample(lattice_params=(a, b, c, alpha, beta, gamma))
    (u, v) = sample.ub_matrix_to_uv(ub_matrix)
    print(f"u={u}")
    print(f"v={v}")


def test_ub_matrix_to_lattice_params(ub_matrix):
    sample = Sample()
    (a, b, c, alpha, beta, gamma) = sample.ub_matrix_to_lattice_params(ub_matrix)
    print(f"a={np.round(a,5)}, b={np.round(b,5)},c={np.round(c,5)}")
    print(f"alpha={int(alpha)}, beta={int(beta)}, gamma={int(gamma)}")


def test_uv_to_ub_matrix(u, v, lattice_params):
    sample = Sample(lattice_params=lattice_params)
    u_matrix, b_matrix, ub_matrix = sample.uv_to_ub_matrix(u, v)
    print(f"u matrix = {u_matrix}")
    print(f"b matrix = {b_matrix}")
    print(f"ub matrix = {ub_matrix}")


if __name__ == "__main__":

    a = 3.574924
    b = 3.574924
    c = 5.663212
    alpha = 90
    beta = 90
    gamma = 120
    ub_matrix = np.array(
        [[0.053821, 0.107638, 0.166485], [0.272815, -0.013290, 0.002566], [0.164330, 0.304247, -0.058788]]
    )
    b_matrix = np.array([[0.3230, 0.1615, 0.0000], [0.0000, 0.2797, -0.0000], [0.0000, 0.0000, 0.1766]])

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]

    # test_b_matrix(a, b, c, alpha, beta, gamma, b_matrix)
    # test_ub_to_uv(a, b, c, alpha, beta, gamma, ub_matrix)
    # test_ub_matrix_to_lattice_params(ub_matrix)
    test_uv_to_ub_matrix(u, v, (a, b, c, alpha, beta, gamma))
