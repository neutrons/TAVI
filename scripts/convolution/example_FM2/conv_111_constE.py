from concurrent.futures import ProcessPoolExecutor
from functools import partial
from time import time

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes

from tavi.instrument.resolution.convolution import convolution
from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample
from tavi.ub_algorithm import (
    UBConf,
    b_mat_from_ub_matrix,
    mantid_to_spice,
    plane_normal_from_two_peaks,
    u_mat_from_ub_matrix,
    uv_to_ub_matrix,
)


def model_disp(qh, qk, ql):
    """return energy for given Q points
    3d FM J=-1 meV S=1, en=6*S*J*(1-cos(Q))
    """

    sj = 5
    # gamma_q = np.cos(2 * np.pi * vq1)
    gamma_q = (np.cos(2 * np.pi * qh) + np.cos(2 * np.pi * qk) + np.cos(2 * np.pi * ql)) / 3

    disp = 2 * sj * (1 - gamma_q)
    disp = np.array((disp - 2, disp + 2))

    # reshape if only one band
    num_disp = len(disp.shape)
    if num_disp == 1:
        disp = np.reshape(disp, (1, np.size(disp)))
    return disp


def model_inten(qh, qk, ql):
    """return intensity for given Q points
    3d FM J=-1 meV S=1, inten = S/2 for all Qs
    """
    inten = np.ones_like(qh, dtype=float) / 2
    inten = np.array((inten, inten))

    # reshape if only one band
    num_inten = len(inten.shape)
    if num_inten == 1:
        inten = np.reshape(inten, (1, np.size(inten)))

    return inten


if __name__ == "__main__":
    # setup instrument
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb3.json"
    hb3 = TAS(fixed_ef=14.7)
    hb3.load_instrument_params_from_json(instrument_config_json_path)
    # set up a cubic sample
    lattice_params = (10, 10, 10, 90, 90, 90)
    sample = Sample(lattice_params)
    sample.set_mosaic(30, 30)
    # set up sample orientation, u along ki, v in plane, u cross v is up
    # when all goniometer angles are zeros
    u = (1, 1, 0)
    v = (0, 0, 1)
    # set up the scattering plane, (100) and (010) in two peaks in plane
    peaks = ((1, 1, 0), (0, 0, 1))
    ub_matrix_mantid = uv_to_ub_matrix(u, v, lattice_params)
    plane_normal_mantid, in_plane_ref_mantid = plane_normal_from_two_peaks(
        u_mat_from_ub_matrix(ub_matrix_mantid),
        b_mat_from_ub_matrix(ub_matrix_mantid),
        *peaks,
    )
    sample.ub_conf = UBConf(
        ub_mat=mantid_to_spice(ub_matrix_mantid),
        plane_normal=mantid_to_spice(plane_normal_mantid),
        in_plane_ref=mantid_to_spice(in_plane_ref_mantid),
    )
    hb3.mount_sample(sample)

    # ----------------------------------------------------
    # points being measured
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = -3, 3, 0.05
    q2_min, q2_max, q2_step = -3, 3, 0.05
    en = 10

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    q2 = np.linspace(q2_min, q2_max, int((q2_max - q2_min) / q2_step) + 1)
    qe_list = np.array([(h, h, l, en) for l in q2 for h in q1])

    reso_params = [
        (reso.hkl, reso.en, reso.r0, reso.mat) if reso is not None else None
        for reso in hb3.cooper_nathans(hkle=qe_list)
    ]
    conv_model = partial(convolution, model_disp=model_disp, model_inten=model_inten)
    t0 = time()
    num_worker = 8
    with ProcessPoolExecutor(max_workers=num_worker) as executor:
        results = executor.map(conv_model, reso_params)
    measurement_inten = np.asarray(list(results))

    print(f"Convolution completed in {(t1 := time()) - t0:.4f} s")

    # ----------------------------------------------------
    # plot 2D contour
    # ----------------------------------------------------
    # calculate and plot resolution
    q1_rez = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / (q1_step * 5)) + 1)
    q2_rez = np.linspace(q2_min, q2_max, int((q2_max - q2_min) / (q2_step * 5)) + 1)
    qe_rez = np.array([(h, h, l, en) for l in q2_rez for h in q1_rez])
    rez_list = hb3.cooper_nathans(hkle=qe_rez, axes=((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en"))

    p = Plot2D()
    for rez in rez_list:
        if rez is None:
            continue
        e_co = rez.get_ellipse(axes=(0, 2), PROJECTION=False)
        e_inco = rez.get_ellipse(axes=(0, 2), PROJECTION=True)
        p.add_reso(e_co, c="w", linestyle="solid")
        p.add_reso(e_inco, c="w", linestyle="dashed")
    # create plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, axes_class=Axes)
    p.plot(ax)

    # overplot contour
    vq1, vq2 = np.meshgrid(q1, q2, indexing="ij")
    im = ax.pcolormesh(
        vq1,
        vq2,
        measurement_inten.reshape(np.shape(vq1)),
        cmap="turbo",
        vmin=0,
        vmax=0.005,
    )
    fig.colorbar(im, ax=ax)

    ax.set_xlim((q1_min, q1_max))
    ax.set_ylim((q2_min, q2_max))

    ax.set_title(
        f"3D FM S=1 J=-5, E = {en} meV"
        + f"\n3D Convolution for {len(q1) * len(q2)} points, "
        + f"completed in {t1 - t0:.3f} s with {num_worker:1d} cores"
    )
    ax.grid(alpha=0.6)
    # plt.tight_layout()
    plt.show()
