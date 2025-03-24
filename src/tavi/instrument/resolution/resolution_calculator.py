from typing import Union

from tavi.instrument.resolution.ellipsoid import ResoEllipsoid


class ResolutionCalculator:
    """
    Resoluatoin Calculator base class

    Attributes:
        instrument (TAS):

    Methods:
        generate_hkle_list
        validate_instrument_parameters

    """

    def __init__(self, instrument):
        self.instrument = instrument

    def generate_hkle_list(
        self,
        hkl: Union[tuple, list[tuple]],
        en: Union[float, list[float]],
    ) -> tuple[tuple[tuple[float, float, float], float, float], ...]:
        """Generate a list containing tuple ((h, k, l), ei, ef)

        Arguments:
            hkl (tuple | list(tuple)): (h,k,l)
            en (float | list(float)): energy trnasfer

        Return:
            hkle_list (tuple): list of (((h, k, l), ei, ef), ...)
        """
        hkle_list = []
        hkl_list = [hkl] if not isinstance(hkl, list) else hkl
        en_list = [en] if not isinstance(en, list) else en

        for en in en_list:
            ei, ef = self.instrument._get_ei_ef(en=en)
            for hkl in hkl_list:
                hkle_list.append((hkl, ei, ef))
        return tuple(hkle_list)

    def validate_instrument_parameters(self):  # noqa: C901
        """Check if enough instrument parameters are provided for Cooper-Nathans mehtod"""

        try:  # monochromator
            mono = self.instrument.monochromator
        except AttributeError:
            print("Monochromator info are missing.")

        if None in (mono_mosaic := (mono.mosaic_h, mono.mosaic_v)):
            raise ValueError("Mosaic of monochromator is missing.")
        elif not all(val > 0 for val in mono_mosaic):
            raise ValueError("Mosaic of monochromator cannot be negative.")

        try:  # analyzer
            ana = self.instrument.analyzer
        except AttributeError:
            print("Analyzer info are missing.")

        if None in (ana_mosaic := (ana.mosaic_h, ana.mosaic_v)):
            raise ValueError("Mosaic of analyzer is missing.")
        elif not all(val > 0 for val in ana_mosaic):
            raise ValueError("Mosaic of analyzer cannot be negative.")

        # collimators
        if (coll := self.instrument.collimators) is None:
            raise ValueError("Collimators info are missing.")
        elif not all(val > 0 for val in coll.horizontal_divergence):
            raise ValueError("Horizontal divergence of collimators cannot be negative.")
        elif not all(val > 0 for val in coll.vertical_divergence):
            raise ValueError("Vertical divergence of collimators cannot be negative.")

        # sample
        if self.instrument.sample is None:
            raise ValueError("Sample info are missing.")

    def test(self):
        return ResoEllipsoid()
