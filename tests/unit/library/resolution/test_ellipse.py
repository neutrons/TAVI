from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.testing.compare import compare_images
from mpl_toolkits.axisartist import Axes

from tavi.library.Instrument.instrument import Instrument
from tavi.library.component.collimators import Collimators
from tavi.library.component.crystal import Crystal
from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.experiment import Experiment
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
from tavi.library.plot.plot_ellipse import Plot, grid_helper
from tavi.library.resolution.resolution import Resolution

matplotlib.use("Agg")  # headless backend so image comparison is reproducible


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[4]

@pytest.fixture
def component():
        # Set up sample
        ol = OrientedLattice()
        ol.a = 5.034785
        ol.b = 5.034785
        ol.c = 13.812004
        ol.alpha, ol.beta, ol.gamma = 90, 90, 90
        ol.UB = np.array([
            [-0.016965, -0.026212, -0.071913],
            [-0.201388, -0.193307, 0.007769],
            [-0.108415, 0.120600,-0.003178]
        ])
        ol.in_plane_ref = np.array([-0.993257, 0.107299, -0.043892])
        ol.plane_normal = np.array([-0.04032, 0.035237, 0.998565])
        sample = Sample(ol=ol)

        # Collimators
        collimators = Collimators()
        collimators.pre_mono_h = 40
        collimators.pre_mono_v = 120

        collimators.pre_sample_h = 100
        collimators.pre_sample_v = 120

        collimators.post_sample_h = 80
        collimators.post_sample_v = 120

        collimators.post_ana_h = 120
        collimators.post_ana_v = 120

        # monochromator
        mono = Crystal()
        mono.mosaic_h = 30
        mono.mosaic_v = 30

        # analyzer
        ana = Crystal()
        ana.mosaic_h = 90
        ana.mosaic_v = 90

        # Instrument
        ins = Instrument(monochromator=mono, 
                         analyzer=ana, 
                         collimators=collimators,
                         goniometer=Goniometer(s2_sense = "+"))
        exp = Experiment()
        return sample, ins, exp

def test_ellipse_in_local_q(component):
    sample, instrument, experiment = component
    hkl = (0, 0, 3)
    res = Resolution("Cooper-Nathans", instrument, sample, experiment, axes = None)
    res_mat, r0 = res.get_resolution(hkl, 4.8, 4.8)
    res_2d, axes_angle = res.get_ellipse(res_mat, (0, 3), PROJECTION=False)
    mat_4d = np.array([
        [ 9583.2034, -4671.0378,    -0.0000,   986.5530],
        [-4671.0378, 21359.5229,     0.0000, -4129.1950],
        [    0.0000,     0.0000,  1607.3928,     0.0000],
        [  986.5530, -4129.1950,    -0.0000,   864.3563],
    ])
    mat_2d = np.array([
          [9583.2034, 986.5530],
          [986.5530, 864.3563],
    ])
    assert np.isclose(axes_angle, 90)
    assert np.isclose(r0, 4.585671082367485e-05, 1e-4)
    assert np.allclose(res_mat, mat_4d, atol = 1e-1)
    assert np.allclose(res_2d, mat_2d, atol = 1e-1)

def test_plot(component):
    sample, instrument, experiment = component
    hkl = (0, 0, 3)
    res1 = Resolution("Cooper-Nathans", instrument, sample, experiment)
    res_4d1,_ = res1.get_resolution(hkl, 4.8, 4.8)
    res_2d_co1, axes_angle1 = res1.get_ellipse(res_4d1, ellipse_axes=(0,1), PROJECTION=False)
    res_2d_inc1, _ = res1.get_ellipse(res_4d1, ellipse_axes=(0,1), PROJECTION=True)
    
    plot1 = Plot(axes_angle1)
    plot1.add_ellipse(res_2d_co1, label = "Coherent")
    plot1.add_ellipse(res_2d_inc1, label = "Incohere", ls= "--")


    res2 = Resolution("Cooper-Nathans", instrument, sample, experiment)
    res_4d2,_ = res2.get_resolution(hkl, 4.8, 4.8)
    res_2d_co2, axes_angle2 = res2.get_ellipse(res_mat=res_4d2, ellipse_axes=(0,3), PROJECTION=False)
    res_2d_inc2, _ = res2.get_ellipse(res_mat = res_4d2, ellipse_axes=(0,3), PROJECTION=True)
    
    plot2 = Plot(axes_angle2)
    plot2.add_ellipse(res_2d_co2, label = "Coherent")
    plot2.add_ellipse(res_2d_inc2, label = "Incohere", ls= "--")

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(1, 2, 1, axes_class=Axes, grid_helper=grid_helper(axes_angle1))
    ax2 = fig.add_subplot(1, 2, 2, axes_class=Axes, grid_helper=grid_helper(axes_angle2))
    ax1.grid(True)
    ax2.grid(True)
    ax1.set_title("H–K plane (coh / incoh)")
    ax2.set_title("H–E plane (coh / incoh)")

    plot1.plot(ax=ax1, show=False)
    plot2.plot(ax=ax2, show=False)

    plt.tight_layout()
    plt.show()
    assert np.isclose(axes_angle1, 60)
    assert np.isclose(axes_angle2, 90)