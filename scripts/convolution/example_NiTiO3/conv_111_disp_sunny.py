from functools import partial
from time import time

import matplotlib.pyplot as plt
import numpy as np
from juliacall import Main as jl  # noqa: N813
from mpl_toolkits.axisartist import Axes

from tavi.data.scan import Scan
from tavi.instrument.resolution.convolution import convolution
from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample


def model_sunny():
    jl.seval("using Sunny, LinearAlgebra")
    jl.include("./scripts/convolution/example_NiTiO3//NiTiO3_R_model.jl")
    return jl.model()


def hexagonal_to_rhombohedral(vq1, vq2, vq3):
    """Convert hexagonal indices to rhombohedral indices."""
    T_H2R = np.array(
        [
            [2 / 3, 1 / 3, 1 / 3],
            [-1 / 3, 1 / 3, 1 / 3],
            [-1 / 3, -2 / 3, 1 / 3],
        ],
        dtype=float,
    ).T
    vq = np.column_stack((vq1, vq2, vq3)).astype(float, copy=False)  # shape (N, 3)
    return list(vq @ T_H2R)  # still (N, 3), now in rhombohedral indices


def _model_disp(vq1, vq2, vq3, model):
    """return energy for given Q points"""
    vq = hexagonal_to_rhombohedral(vq1, vq2, vq3)
    return np.array(jl.disp(model, vq), order="C")


def _model_inten(vq1, vq2, vq3, model):
    """return intensity for given Q points"""
    vq = hexagonal_to_rhombohedral(vq1, vq2, vq3)
    return np.array(jl.inten(model, vq), order="C")


if __name__ == "__main__":
    # setup instrument
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    ctax = TAS(fixed_ef=4.8)
    ctax.load_instrument_params_from_json(instrument_config_json_path)
    path_to_spice_folder = "test_data/exp424/"
    sample = Sample.from_scan(Scan.from_spice(path_to_spice_folder, scan_num=42))
    sample.mosaic_h, sample.mosaic_v = 180, 180
    ctax.mount_sample(sample)

    # ----------------------------------------------------
    # points being measured
    # ----------------------------------------------------
    ql_min, ql_max, ql_step = 2.5, 3.9, 0.1
    en_min, en_max, en_step = 0.1, 4.0, 0.1

    ql_list = np.linspace(ql_min, ql_max, round((ql_max - ql_min) / ql_step + 1))
    en_list = np.linspace(en_min, en_max, round((en_max - en_min) / en_step + 1))
    qe_list = np.array([(0, 0, ql, en) for ql in ql_list for en in en_list])

    axes = ((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en")
    # axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1), "en")
    reso_params = [
        (reso.hkl, reso.en, reso.r0, reso.mat) if reso is not None else None
        for reso in ctax.cooper_nathans(hkle=qe_list, axes=axes)
    ]
    model = model_sunny()
    model_disp = partial(_model_disp, model=model)
    model_inten = partial(_model_inten, model=model)
    conv_model = partial(convolution, model_disp=model_disp, model_inten=model_inten)
    t0 = time()
    # ------------------- muliti core ------------------
    # num_worker = 8
    # with ProcessPoolExecutor(max_workers=num_worker) as executor:
    #     results = executor.map(conv_model, reso_params)
    # measurement_inten = np.asarray(list(results))
    # ------------------- single core ------------------
    num_worker = 1
    sz = len(reso_params)
    measurement_inten = np.empty(shape=sz)
    for i in range(sz):
        measurement_inten[i] = conv_model(reso_params[i])
        # measurement_inten[i] = reso_params[i][2]
    # --------------------------------------------------
    print(f"Convolution completed in {(t1 := time()) - t0:.4f} s")

    # ----------------------------------------------------
    # plot 2D contour
    # ----------------------------------------------------
    # calculate and plot resolution
    ql_rez = np.linspace(ql_min, ql_max, int((ql_max - ql_min) / (ql_step * 5)) + 1)
    en_rez = np.linspace(en_min, en_max, int((en_max - en_min) / (en_step * 5)) + 1)
    qe_rez = np.array([(0, 0, ql, en) for ql in ql_rez for en in en_rez])
    rez_list = ctax.cooper_nathans(hkle=qe_rez, axes=((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en"))

    p = Plot2D()
    for rez in rez_list:
        if rez is None:
            continue
        e_co = rez.get_ellipse(axes=(2, 3), PROJECTION=False)
        e_inco = rez.get_ellipse(axes=(2, 3), PROJECTION=True)
        p.add_reso(e_co, c="w", linestyle="solid")
        p.add_reso(e_inco, c="w", linestyle="dashed")

    # create plot
    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes)
    p.plot(ax)

    # overplot contour
    vq1, ven = np.meshgrid(ql_list, en_list, indexing="ij")
    im = ax.pcolormesh(
        vq1,
        ven,
        measurement_inten.reshape(np.shape(vq1)),
        cmap="turbo",
        vmin=0,
        vmax=1e-4,
    )
    fig.colorbar(im, ax=ax)

    # plot dispersion
    disp = model_disp(np.zeros_like(ql_list), np.zeros_like(ql_list), ql_list)
    for i in range(np.shape(disp)[0]):
        ax.plot(ql_list, disp[i], "-w")

    ax.set_xlim((ql_min - 0.1, ql_max + 0.1))
    ax.set_ylim((en_min - 0.1, en_max))

    ax.set_title(
        f"3D Convolution for {len(ql_list) * len(en_list)} points, "
        + f"completed in {t1 - t0:.3f} s"
        + f" with {num_worker:1d} cores"
    )
    ax.grid(alpha=0.6)
    # plt.tight_layout()
    plt.show()
