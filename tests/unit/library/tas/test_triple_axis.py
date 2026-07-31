import pytest

from tavi.library.Instrument.instrument import Instrument
from tavi.library.component.collimators import Collimators
from tavi.library.component.crystal import Crystal
from tavi.library.component.goniometer import Goniometer
from tavi.library.experiment.experiment import Experiment
from tavi.library.experiment.peak import DataPoint, MotorAngles
from tavi.library.experiment.utilities import mantid_to_spice
from tavi.library.geometry.oriented_lattice import OrientedLattice
from tavi.library.geometry.sample import Sample
import numpy as np

from tavi.library.plot.plot_ellipse import PlotResolution
from tavi.library.resolution.ellipsoid import ResolutionEllipsoid
from tavi.library.resolution.resolution import Resolution, ResolutionModel
from tavi.library.tas.triple_axis import TAS

@pytest.fixture
def component():
        # Set up sample
        ol = OrientedLattice()
        ol.a = 6.4233
        ol.b = 6.4229
        ol.c = 6.4222
        ol.alpha, ol.beta, ol.gamma = 90, 90, 90
        
        sample = Sample(ol=ol)

        # Collimators
        collimators = Collimators()
        collimators.pre_mono_h = 50
        collimators.pre_sample_h = 40
        collimators.post_sample_h = 40
        collimators.post_ana_h = 80

        collimators.pre_mono_v = 600
        collimators.pre_sample_v = 600
        collimators.post_sample_v = 600
        collimators.post_ana_v = 600

        sample.mosaic_h = 5
        sample.mosaic_v = 5


        # monochromator
        mono = Crystal()
        mono.mosaic_h = 60
        mono.mosaic_v = 30
        mono.sense = "+"

        # analyzer
        ana = Crystal()
        ana.mosaic_h = 60
        ana.mosaic_v = 30
        ana.sense = "+"

        # Instrument
        ins = Instrument(monochromator=mono, 
                         analyzer=ana, 
                         collimators=collimators,
                         goniometer=Goniometer(type= "Y,Z,Y,bisect", s2_sense = "-"))
        exp = Experiment()
        return sample, ins, exp

@pytest.fixture
def generate_peaks():
    ei = 14.4503
    ef = 14.4503
    peaks = (
        DataPoint(hkl=(-3.0, -2.0, -1.001), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.99, "chi": -42.02, "phi": 99.66})),
        DataPoint(hkl=(-3.0, -2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.81, "omega": -41.99, "chi": -44.89, "phi": 78.59})),
        DataPoint(hkl=(-3.0, -2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.94, "chi": -43.69, "phi": 56.93})),
        DataPoint(hkl=(-3.0, -1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.76, "omega": -43.91, "chi": -45.34, "phi": 129.6})),
        DataPoint(hkl=(-3.0, -1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.95, "chi": -48.95, "phi": 26.07})),
        DataPoint(hkl=(-2.0, -3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.96, "chi": -20.79, "phi": 93.85})),
        DataPoint(hkl=(-2.0, -3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.82, "omega": -41.97, "chi": -22.3, "phi": 77.25})),
        DataPoint(hkl=(-2.0, -3.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.94, "chi": -22.11, "phi": 60.49})),
        DataPoint(hkl=(-2.0, -2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.85, "omega": -39.98, "chi": -25.38, "phi": 117.5})),
        DataPoint(hkl=(-2.0, -2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.8, "chi": -30.55, "phi": 100.6})),
        DataPoint(hkl=(-2.0, -2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.68, "chi": -33.6, "phi": 77.83})),
        DataPoint(hkl=(-2.0, -2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.86, "chi": -32.35, "phi": 54.62})),
        DataPoint(hkl=(-2.0, -2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.86, "omega": -39.99, "chi": -28.36, "phi": 36.88})),
        DataPoint(hkl=(-2.0, -1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.95, "chi": -26.02, "phi": 142.2})),
        DataPoint(hkl=(-2.0, -1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.84, "chi": -34.1, "phi": 132.7})),
        DataPoint(hkl=(-2.0, -1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.97, "omega": -27.07, "chi": -44.67, "phi": 114.2})),
        DataPoint(hkl=(-2.0, -1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -27.04, "chi": -47.37, "phi": 42.26})),
        DataPoint(hkl=(-2.0, -1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.78, "chi": -37.88, "phi": 21.79})),
        DataPoint(hkl=(-2.0, -1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.92, "chi": -30.2, "phi": 11.46})),
        DataPoint(hkl=(-2.0, 0.0, -2.999), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.81, "omega": -41.93, "chi": -30.66, "phi": 159.0})),
        DataPoint(hkl=(-2.0, 0.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.63, "chi": -41.62, "phi": 155.6})),
        DataPoint(hkl=(-2.0, 0.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.62, "chi": -46.14, "phi": -2.091})),
        DataPoint(hkl=(-2.0, 0.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.94, "chi": -35.22, "phi": -6.001})),
        DataPoint(hkl=(-2.0, 1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.93, "chi": -32.96, "phi": 177.1})),
        DataPoint(hkl=(-2.0, 1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.52, "omega": -33.71, "chi": -43.78, "phi": -178.0})),
        DataPoint(hkl=(-2.0, 1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.8, "chi": -48.19, "phi": -30.73})),
        DataPoint(hkl=(-2.0, 1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.89, "chi": -37.47, "phi": -25.07})),
        DataPoint(hkl=(-2.0, 2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.85, "omega": -39.89, "chi": -41.0, "phi": -156.9})),
        DataPoint(hkl=(-2.0, 2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.85, "omega": -39.89, "chi": -44.62, "phi": -53.13})),
        DataPoint(hkl=(-2.0, 3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.85, "chi": -42.12, "phi": -127.1})),
        DataPoint(hkl=(-2.0, 3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.81, "omega": -41.82, "chi": -45.01, "phi": -106.0})),
        DataPoint(hkl=(-2.0, 3.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.85, "chi": -43.81, "phi": -84.29})),
        DataPoint(hkl=(-1.0, -3.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.87, "chi": -4.728, "phi": 109.0})),
        DataPoint(hkl=(-1.0, -3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.82, "omega": -37.95, "chi": -6.033, "phi": 94.22})),
        DataPoint(hkl=(-1.0, -3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.72, "omega": -35.9, "chi": -7.06, "phi": 76.59})),
        DataPoint(hkl=(-1.0, -3.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.97, "chi": -7.428, "phi": 58.9})),
        DataPoint(hkl=(-1.0, -3.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.97, "chi": -7.198, "phi": 44.02})),
        DataPoint(hkl=(-1.0, -2.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.97, "chi": -7.141, "phi": 130.8})),
        DataPoint(hkl=(-1.0, -2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.82, "chi": -9.698, "phi": 119.4})),
        DataPoint(hkl=(-1.0, -2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.97, "omega": -27.06, "chi": -12.87, "phi": 101.7})),
        DataPoint(hkl=(-1.0, -2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.94, "omega": -24.53, "chi": -15.18, "phi": 76.93})),
        DataPoint(hkl=(-1.0, -2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.99, "omega": -27.0, "chi": -14.8, "phi": 51.97})),
        DataPoint(hkl=(-1.0, -2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.81, "chi": -12.82, "phi": 33.84})),
        DataPoint(hkl=(-1.0, -2.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.86, "chi": -10.87, "phi": 22.27})),
        DataPoint(hkl=(-1.0, -1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.97, "chi": -11.52, "phi": 145.1})),
        DataPoint(hkl=(-1.0, -1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -26.98, "chi": -16.66, "phi": 136.2})),
        DataPoint(hkl=(-1.0, -1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.43, "omega": -18.78, "chi": -25.38, "phi": 117.5})),
        DataPoint(hkl=(-1.0, -1.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.28, "chi": -33.6, "phi": 77.83})),
        DataPoint(hkl=(-1.0, -1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.44, "omega": -18.79, "chi": -28.36, "phi": 36.88})),
        DataPoint(hkl=(-1.0, -1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.99, "omega": -27.08, "chi": -20.62, "phi": 17.21})),
        DataPoint(hkl=(-1.0, -1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.96, "chi": -15.8, "phi": 7.957})),
        DataPoint(hkl=(-1.0, 0.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.88, "chi": -15.76, "phi": 162.6})),
        DataPoint(hkl=(-1.0, 0.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.52, "chi": -23.71, "phi": 160.8})),
        DataPoint(hkl=(-1.0, 0.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.21, "chi": -41.62, "phi": 155.6})),
        DataPoint(hkl=(-1.0, 0.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.24, "chi": -46.14, "phi": -2.09})),
        DataPoint(hkl=(-1.0, 0.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.5, "chi": -28.29, "phi": -7.949})),
        DataPoint(hkl=(-1.0, 0.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.89, "chi": -20.35, "phi": -9.881})),
        DataPoint(hkl=(-1.0, 1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.84, "omega": -37.95, "chi": -18.56, "phi": -179.3})),
        DataPoint(hkl=(-1.0, 1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -27.04, "chi": -26.58, "phi": -172.7})),
        DataPoint(hkl=(-1.0, 1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.43, "omega": -18.67, "chi": -41.0, "phi": -156.9})),
        DataPoint(hkl=(-1.0, 1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.43, "omega": -18.69, "chi": -44.62, "phi": -53.13})),
        DataPoint(hkl=(-1.0, 1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -26.94, "chi": -30.86, "phi": -35.56})),
        DataPoint(hkl=(-1.0, 1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.94, "chi": -23.01, "phi": -28.55})),
        DataPoint(hkl=(-1.0, 2.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.79, "omega": -43.84, "chi": -19.56, "phi": -163.7})),
        DataPoint(hkl=(-1.0, 2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.79, "chi": -25.53, "phi": -153.0})),
        DataPoint(hkl=(-1.0, 2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.97, "omega": -26.89, "chi": -32.98, "phi": -134.6})),
        DataPoint(hkl=(-1.0, 2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.43, "chi": -37.9, "phi": -105.5})),
        DataPoint(hkl=(-1.0, 2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -26.99, "chi": -35.24, "phi": -75.53})),
        DataPoint(hkl=(-1.0, 2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.76, "chi": -28.98, "phi": -55.91})),
        DataPoint(hkl=(-1.0, 2.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.95, "chi": -23.53, "phi": -44.62})),
        DataPoint(hkl=(-1.0, 3.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.93, "chi": -23.47, "phi": -140.6})),
        DataPoint(hkl=(-1.0, 3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.9, "chi": -27.48, "phi": -124.9})),
        DataPoint(hkl=(-1.0, 3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.72, "omega": -35.82, "chi": -29.78, "phi": -105.0})),
        DataPoint(hkl=(-1.0, 3.0, 1.001), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.89, "chi": -29.05, "phi": -84.85})),
        DataPoint(hkl=(-1.0, 3.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.81, "chi": -26.18, "phi": -68.49})),
        DataPoint(hkl=(0.0, -3.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.82, "omega": -41.97, "chi": 10.73, "phi": 110.2})),
        DataPoint(hkl=(0.0, -3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.88, "chi": 11.51, "phi": 94.65})),
        DataPoint(hkl=(0.0, -3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.79, "chi": 11.36, "phi": 75.84})),
        DataPoint(hkl=(0.0, -3.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.74, "omega": -35.89, "chi": 10.03, "phi": 57.12})),
        DataPoint(hkl=(0.0, -3.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.86, "chi": 8.14, "phi": 41.79})),
        DataPoint(hkl=(0.0, -2.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.88, "chi": 8.201, "phi": 133.0})),
        DataPoint(hkl=(0.0, -2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.62, "chi": 9.651, "phi": 121.6})),
        DataPoint(hkl=(0.0, -2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.94, "omega": -24.43, "chi": 11.19, "phi": 102.9})),
        DataPoint(hkl=(0.0, -2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -43.5, "omega": -21.76, "chi": 11.36, "phi": 75.84})),
        DataPoint(hkl=(0.0, -2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.54, "chi": 9.102, "phi": 48.93})),
        DataPoint(hkl=(0.0, -2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.62, "chi": 6.366, "phi": 30.53})),
        DataPoint(hkl=(0.0, -2.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.84, "omega": -41.99, "chi": 4.35, "phi": 19.35})),
        DataPoint(hkl=(0.0, -1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.91, "chi": 5.76, "phi": 148.2})),
        DataPoint(hkl=(0.0, -1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.94, "omega": -24.55, "chi": 7.122, "phi": 140.1})),
        DataPoint(hkl=(0.0, -1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.37, "omega": -15.23, "chi": 9.651, "phi": 121.6})),
        DataPoint(hkl=(0.0, -1.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -21.35, "omega": -10.68, "chi": 11.36, "phi": 75.84})),
        DataPoint(hkl=(0.0, -1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.23, "chi": 6.366, "phi": 30.53})),
        DataPoint(hkl=(0.0, -1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.96, "omega": -24.46, "chi": 2.99, "phi": 12.35})),
        DataPoint(hkl=(0.0, -1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.74, "omega": -35.91, "chi": 1.388, "phi": 4.371})),
        DataPoint(hkl=(0.0, 0.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.79, "chi": 2.3, "phi": 166.3})),
        DataPoint(hkl=(0.0, 0.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -43.5, "omega": -21.79, "chi": 2.3, "phi": 166.3})),
        DataPoint(hkl=(0.0, 0.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -21.35, "omega": -10.67, "chi": 2.3, "phi": 166.3})),
        DataPoint(hkl=(0.0, 0.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -21.35, "omega": -10.8, "chi": -2.3, "phi": -13.69})),
        DataPoint(hkl=(0.0, 0.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -43.51, "omega": -21.79, "chi": -2.298, "phi": -13.69})),
        DataPoint(hkl=(0.0, 0.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.55, "omega": -33.83, "chi": -2.298, "phi": -13.69})),
        DataPoint(hkl=(0.0, 1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.74, "omega": -35.93, "chi": -1.386, "phi": -175.6})),
        DataPoint(hkl=(0.0, 1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.47, "chi": -2.992, "phi": -167.6})),
        DataPoint(hkl=(0.0, 1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.24, "chi": -6.365, "phi": -149.5})),
        DataPoint(hkl=(0.0, 1.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -21.35, "omega": -10.67, "chi": -11.36, "phi": -104.2})),
        DataPoint(hkl=(0.0, 1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.23, "chi": -9.649, "phi": -58.37})),
        DataPoint(hkl=(0.0, 1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.56, "chi": -7.12, "phi": -39.91})),
        DataPoint(hkl=(0.0, 1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.91, "chi": -5.758, "phi": -31.84})),
        DataPoint(hkl=(0.0, 2.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.84, "omega": -42.0, "chi": -4.349, "phi": -160.6})),
        DataPoint(hkl=(0.0, 2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.64, "chi": -6.365, "phi": -149.5})),
        DataPoint(hkl=(0.0, 2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.96, "omega": -24.54, "chi": -9.103, "phi": -131.1})),
        DataPoint(hkl=(0.0, 2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -43.5, "omega": -21.76, "chi": -11.36, "phi": -104.2})),
        DataPoint(hkl=(0.0, 2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.44, "chi": -11.19, "phi": -77.05})),
        DataPoint(hkl=(0.0, 2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.62, "chi": -9.649, "phi": -58.37})),
        DataPoint(hkl=(0.0, 2.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.88, "chi": -8.199, "phi": -47.01})),
        DataPoint(hkl=(0.0, 3.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.88, "chi": -8.139, "phi": -138.2})),
        DataPoint(hkl=(0.0, 3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.89, "chi": -10.03, "phi": -122.9})),
        DataPoint(hkl=(0.0, 3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.78, "chi": -11.36, "phi": -104.2})),
        DataPoint(hkl=(0.0, 3.0, 1.001), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.74, "omega": -35.89, "chi": -11.51, "phi": -85.34})),
        DataPoint(hkl=(0.0, 3.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.98, "chi": -10.73, "phi": -69.81})),
        DataPoint(hkl=(1.0, -3.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.82, "chi": 26.18, "phi": 111.5})),
        DataPoint(hkl=(1.0, -3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.89, "chi": 29.05, "phi": 95.14})),
        DataPoint(hkl=(1.0, -3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.72, "omega": -35.82, "chi": 29.78, "phi": 74.99})),
        DataPoint(hkl=(1.0, -3.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.91, "chi": 27.48, "phi": 55.14})),
        DataPoint(hkl=(1.0, -3.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.92, "chi": 23.47, "phi": 39.39})),
        DataPoint(hkl=(1.0, -2.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.95, "chi": 23.53, "phi": 135.4})),
        DataPoint(hkl=(1.0, -2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.77, "chi": 28.98, "phi": 124.1})),
        DataPoint(hkl=(1.0, -2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.97, "omega": -26.99, "chi": 35.25, "phi": 104.5})),
        DataPoint(hkl=(1.0, -2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.94, "omega": -24.44, "chi": 37.9, "phi": 74.51})),
        DataPoint(hkl=(1.0, -2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -26.89, "chi": 32.98, "phi": 45.43})),
        DataPoint(hkl=(1.0, -2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.78, "chi": 25.53, "phi": 26.96})),
        DataPoint(hkl=(1.0, -2.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.83, "chi": 19.56, "phi": 16.32})),
        DataPoint(hkl=(1.0, -1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.95, "chi": 23.02, "phi": 151.4})),
        DataPoint(hkl=(1.0, -1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -26.95, "chi": 30.87, "phi": 144.4})),
        DataPoint(hkl=(1.0, -1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.44, "omega": -18.73, "chi": 44.62, "phi": 126.9})),
        DataPoint(hkl=(1.0, -1.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.12, "chi": 56.3, "phi": 72.85})),
        DataPoint(hkl=(1.0, -1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.44, "omega": -18.7, "chi": 41.0, "phi": 23.12})),
        DataPoint(hkl=(1.0, -1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -27.03, "chi": 26.58, "phi": 7.256})),
        DataPoint(hkl=(1.0, -1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.84, "omega": -37.94, "chi": 18.57, "phi": 0.732})),
        DataPoint(hkl=(1.0, 0.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.91, "chi": 20.35, "phi": 170.1})),
        DataPoint(hkl=(1.0, 0.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.52, "chi": 28.29, "phi": 172.1})),
        DataPoint(hkl=(1.0, 0.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.38, "omega": -15.26, "chi": 46.14, "phi": 177.9})),
        DataPoint(hkl=(1.0, 0.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.37, "omega": -15.23, "chi": 41.62, "phi": -24.44})),
        DataPoint(hkl=(1.0, 0.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.54, "chi": 23.71, "phi": -19.21})),
        DataPoint(hkl=(1.0, 0.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.73, "omega": -35.9, "chi": 15.76, "phi": -17.4})),
        DataPoint(hkl=(1.0, 1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.84, "omega": -37.98, "chi": 15.8, "phi": -172.0})),
        DataPoint(hkl=(1.0, 1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -27.1, "chi": 20.62, "phi": -162.8})),
        DataPoint(hkl=(1.0, 1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.44, "omega": -18.81, "chi": 28.36, "phi": -143.1})),
        DataPoint(hkl=(1.0, 1.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -30.37, "omega": -15.29, "chi": 33.6, "phi": -102.2})),
        DataPoint(hkl=(1.0, 1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -37.44, "omega": -18.81, "chi": 25.38, "phi": -62.49})),
        DataPoint(hkl=(1.0, 1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -26.99, "chi": 16.66, "phi": -43.8})),
        DataPoint(hkl=(1.0, 1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.98, "chi": 11.52, "phi": -34.93})),
        DataPoint(hkl=(1.0, 2.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.79, "omega": -43.89, "chi": 10.87, "phi": -157.7})),
        DataPoint(hkl=(1.0, 2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.54, "omega": -33.84, "chi": 12.82, "phi": -146.2})),
        DataPoint(hkl=(1.0, 2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.97, "omega": -27.01, "chi": 14.8, "phi": -128.0})),
        DataPoint(hkl=(1.0, 2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.95, "omega": -24.55, "chi": 15.18, "phi": -103.1})),
        DataPoint(hkl=(1.0, 2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -27.08, "chi": 12.87, "phi": -78.33})),
        DataPoint(hkl=(1.0, 2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.84, "chi": 9.698, "phi": -60.55})),
        DataPoint(hkl=(1.0, 2.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.99, "chi": 7.143, "phi": -49.22})),
        DataPoint(hkl=(1.0, 3.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.99, "chi": 7.201, "phi": -136.0})),
        DataPoint(hkl=(1.0, 3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.97, "chi": 7.429, "phi": -121.1})),
        DataPoint(hkl=(1.0, 3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -71.72, "omega": -35.91, "chi": 7.06, "phi": -103.4})),
        DataPoint(hkl=(1.0, 3.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.83, "omega": -37.97, "chi": 6.033, "phi": -85.77})),
        DataPoint(hkl=(1.0, 3.0, 2.001), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.88, "chi": 4.729, "phi": -71.0})),
        DataPoint(hkl=(2.0, -3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.76, "omega": -43.85, "chi": 43.81, "phi": 95.71})),
        DataPoint(hkl=(2.0, -3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.81, "omega": -41.83, "chi": 45.01, "phi": 74.0})),
        DataPoint(hkl=(2.0, -3.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.85, "chi": 42.13, "phi": 52.91})),
        DataPoint(hkl=(2.0, -2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.85, "omega": -39.91, "chi": 44.62, "phi": 126.9})),
        DataPoint(hkl=(2.0, -2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.67, "chi": 52.91, "phi": 106.3})),
        DataPoint(hkl=(2.0, -2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.2, "omega": -31.52, "chi": 56.3, "phi": 72.85})),
        DataPoint(hkl=(2.0, -2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.72, "chi": 50.44, "phi": 41.38})),
        DataPoint(hkl=(2.0, -2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.86, "omega": -39.89, "chi": 41.0, "phi": 23.12})),
        DataPoint(hkl=(2.0, -1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.91, "chi": 37.47, "phi": 154.9})),
        DataPoint(hkl=(2.0, -1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.82, "chi": 48.2, "phi": 149.3})),
        DataPoint(hkl=(2.0, -1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.97, "omega": -26.96, "chi": 59.75, "phi": 14.67})),
        DataPoint(hkl=(2.0, -1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.71, "chi": 43.78, "phi": 2.008})),
        DataPoint(hkl=(2.0, -1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.93, "chi": 32.96, "phi": -2.934})),
        DataPoint(hkl=(2.0, 0.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.97, "chi": 35.22, "phi": 174.0})),
        DataPoint(hkl=(2.0, 0.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.21, "omega": -31.65, "chi": 46.14, "phi": 177.9})),
        DataPoint(hkl=(2.0, 0.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.94, "omega": -24.51, "chi": 59.12, "phi": -33.78})),
        DataPoint(hkl=(2.0, 0.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.2, "omega": -31.64, "chi": 41.62, "phi": -24.43})),
        DataPoint(hkl=(2.0, 0.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.95, "chi": 30.66, "phi": -20.99})),
        DataPoint(hkl=(2.0, 1.0, -3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.94, "chi": 30.2, "phi": -168.5})),
        DataPoint(hkl=(2.0, 1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.81, "chi": 37.89, "phi": -158.2})),
        DataPoint(hkl=(2.0, 1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -27.06, "chi": 47.37, "phi": -137.7})),
        DataPoint(hkl=(2.0, 1.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -48.94, "omega": -24.58, "chi": 52.01, "phi": -100.7})),
        DataPoint(hkl=(2.0, 1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -53.98, "omega": -27.09, "chi": 44.67, "phi": -65.79})),
        DataPoint(hkl=(2.0, 1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.86, "chi": 34.1, "phi": -47.28})),
        DataPoint(hkl=(2.0, 1.0, 3.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.96, "chi": 26.02, "phi": -37.82})),
        DataPoint(hkl=(2.0, 2.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.86, "omega": -40.01, "chi": 28.36, "phi": -143.1})),
        DataPoint(hkl=(2.0, 2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.88, "chi": 32.35, "phi": -125.4})),
        DataPoint(hkl=(2.0, 2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -63.2, "omega": -31.7, "chi": 33.6, "phi": -102.2})),
        DataPoint(hkl=(2.0, 2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -67.53, "omega": -33.81, "chi": 30.55, "phi": -79.41})),
        DataPoint(hkl=(2.0, 2.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -79.86, "omega": -40.0, "chi": 25.38, "phi": -62.49})),
        DataPoint(hkl=(2.0, 3.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.95, "chi": 22.11, "phi": -119.5})),
        DataPoint(hkl=(2.0, 3.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.99, "chi": 22.3, "phi": -102.7})),
        DataPoint(hkl=(2.0, 3.0, 1.001), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.97, "chi": 20.79, "phi": -86.15})),
        DataPoint(hkl=(3.0, -1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.87, "chi": 59.26, "phi": 154.6})),
        DataPoint(hkl=(3.0, -1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.89, "chi": 54.74, "phi": -3.316})),
        DataPoint(hkl=(3.0, 0.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.95, "chi": 56.85, "phi": -176.2})),
        DataPoint(hkl=(3.0, 0.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.83, "omega": -41.98, "chi": 52.45, "phi": -29.3})),
        DataPoint(hkl=(3.0, 1.0, -2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.78, "omega": -43.99, "chi": 48.95, "phi": -153.9})),
        DataPoint(hkl=(3.0, 1.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.82, "omega": -38.0, "chi": 57.01, "phi": -133.2})),
        DataPoint(hkl=(3.0, 1.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -75.82, "omega": -37.99, "chi": 54.55, "phi": -68.47})),
        DataPoint(hkl=(3.0, 1.0, 2.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.91, "chi": 45.34, "phi": -50.39})),
        DataPoint(hkl=(3.0, 2.0, -1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.95, "chi": 43.7, "phi": -123.1})),
        DataPoint(hkl=(3.0, 2.0, 0.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -83.82, "omega": -41.99, "chi": 44.89, "phi": -101.4})),
        DataPoint(hkl=(3.0, 2.0, 1.0), ei=ei, ef=ef, angles=MotorAngles(angles_dict={"two_theta": -87.77, "omega": -43.99, "chi": 42.02, "phi": -80.34})),
    )
    return peaks

def test_tas_ub(component, generate_peaks):

    sample, instrument, experiment = component
    peaks = generate_peaks
    tas = TAS(instrument=instrument, sample = sample, experiment = experiment)
    ei = 14.4643
    ef = 14.4643
    ub_spice = mantid_to_spice(tas.ub(peaks=(peaks[0], peaks[3])), version="old")
    assert np.allclose(ub_spice, np.array([
        [-0.00139387, -0.03772994,  0.15106212],
        [ 0.03124861,  0.14790957,  0.03723904],
        [-0.15250851,  0.03065114,  0.00624953]]))
    
    resolution = Resolution(
        model=ResolutionModel.CooperNathans,
        instrument=instrument,
        sample=sample,
        experiment=experiment,
        resolution_frame="local",
    )
    res_4d, r0 = resolution.get_resolution(hkl = (-2, -3, 0), ei=ei, ef=ef,rot_mat=None)
    ellipsoid = ResolutionEllipsoid(res_mat=res_4d, r0 = r0, resolution_frame = "local")

    assert np.isclose(ellipsoid.coh_fwhm(projection_axis=1), 0.021729, atol = 1e-6)

    # ellipse, axes_angle = tas.resolution.get_ellipse(res_mat=res_4d, ellipse_axes=(0,1))

    # plot = PlotResolution(axes_angle=axes_angle)
    # plot.add_ellipse(mat = ellipse)
    # plot.plot()