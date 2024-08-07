import numpy as np
from tavi.sample.sample import Sample
from tavi.sample.xtal import Xtal


def test_b_matrix(lattice_params, b_matrix):

    sample = Sample(lattice_params=lattice_params)
    print(sample.b_mat())
    print(np.allclose(sample.b_mat(), b_matrix, atol=1e-4))


def test_ub_matrix_to_uv(lattice_params, ub_matrix):

    xtal = Xtal(lattice_params=lattice_params)
    (u, v) = xtal.ub_matrix_to_uv(ub_matrix)
    print(f"u={u}")
    print(f"v={v}")


def test_ub_matrix_to_lattice_params(ub_matrix):
    sample = Xtal()
    (a, b, c, alpha, beta, gamma) = sample.ub_matrix_to_lattice_params(ub_matrix)
    print(f"a={np.round(a,5)}, b={np.round(b,5)},c={np.round(c,5)}")
    print(f"alpha={np.round(alpha,3)}, beta={np.round(beta,3)}, gamma={np.round(gamma,3)}")


def test_uv_to_ub_matrix(u, v, lattice_params):
    sample = Xtal(lattice_params=lattice_params)
    ub_matrix = sample.uv_to_ub_matrix(u, v)

    print(f"ub matrix = {ub_matrix}")


if __name__ == "__main__":

    a = 3.574924
    b = 3.574924
    c = 5.663212
    alpha = 90
    beta = 90
    gamma = 120

    # ub_matrix = np.array(
    #     [
    #         [0.053821, 0.107638, 0.166485],
    #         [0.272815, -0.013290, 0.002566],
    #         [0.164330, 0.304247, -0.058788],
    #     ]
    # )
    ub_matrix = np.array(
        [
            [0.0538, 0.1076, 0.1665],
            [0.2728, -0.0133, 0.0026],
            [0.1643, 0.3042, -0.0588],
        ]
    )
    b_matrix = np.array(
        [
            [0.3230, 0.1615, 0.0000],
            [0.0000, 0.2797, -0.0000],
            [0.0000, 0.0000, 0.1766],
        ]
    )

    u = [0.15623, 2.83819, -1.88465]
    v = [-0.00060, 1.03219, 5.33915]

    lattice_params = (a, b, c, alpha, beta, gamma)

    # test_b_matrix(lattice_params, b_matrix)
    test_ub_matrix_to_uv(lattice_params, ub_matrix)
    # test_ub_matrix_to_lattice_params(ub_matrix)
    test_uv_to_ub_matrix(u, v, lattice_params)
