"""Collimators."""

import numpy as np


class Collimators:
    """Collimators that holds horizontal and vertical divergence for pre_mono, pre_sample, post_sample, post_ana in minute arc."""

    def __init__(
        self,
        pre_mono_h: float = 1.0,
        pre_mono_v: float = 1.0,
        pre_sample_h: float = 1.0,
        pre_sample_v: float = 1.0,
        post_sample_h: float = 1.0,
        post_sample_v: float = 1.0,
        post_ana_h: float = 1.0,
        post_ana_v: float = 1.0,
    ) -> None:
        """Init."""
        self._pre_mono_h = pre_mono_h
        self._pre_mono_v = pre_mono_v
        self._pre_sample_h = pre_sample_h
        self._pre_sample_v = pre_sample_v
        self._post_sample_h = post_sample_h
        self._post_sample_v = post_sample_v
        self._post_ana_h = post_ana_h
        self._post_ana_v = post_ana_v

    @property
    def pre_mono_h(self) -> float:
        """Getter."""
        return np.radians(self._pre_mono_h / 60)

    @pre_mono_h.setter
    def pre_mono_h(self, val: float) -> None:
        """Setter."""
        self._pre_mono_h = val

    @property
    def pre_mono_v(self) -> float:
        """Getter."""
        return np.radians(self._pre_mono_v / 60)

    @pre_mono_v.setter
    def pre_mono_v(self, val: float) -> None:
        """Setter."""
        self._pre_mono_v = val

    @property
    def pre_sample_h(self) -> float:
        """Getter."""
        return np.radians(self._pre_sample_h / 60)

    @pre_sample_h.setter
    def pre_sample_h(self, val: float) -> None:
        """Setter."""
        self._pre_sample_h = val

    @property
    def pre_sample_v(self) -> float:
        """Getter."""
        return np.radians(self._pre_sample_v / 60)

    @pre_sample_v.setter
    def pre_sample_v(self, val: float) -> None:
        """Setter."""
        self._pre_sample_v = val

    @property
    def post_sample_h(self) -> float:
        """Getter."""
        return np.radians(self._post_sample_h / 60)

    @post_sample_h.setter
    def post_sample_h(self, val: float) -> None:
        """Setter."""
        self._post_sample_h = val

    @property
    def post_sample_v(self) -> float:
        """Getter."""
        return np.radians(self._post_sample_v / 60)

    @post_sample_v.setter
    def post_sample_v(self, val: float) -> None:
        """Setter."""
        self._post_sample_v = val

    @property
    def post_ana_h(self) -> float:
        """Getter."""
        return np.radians(self._post_ana_h / 60)

    @post_ana_h.setter
    def post_ana_h(self, val: float) -> None:
        """Setter."""
        self._post_ana_h = val

    @property
    def post_ana_v(self) -> float:
        """Getter."""
        return np.radians(self._post_ana_v / 60)

    @post_ana_v.setter
    def post_ana_v(self, val: float) -> None:
        """Setter."""
        self._post_ana_v = val
