# -*- coding: utf-8 -*-
import json
from pathlib import Path
from typing import Optional, Union

from tavi.data.scan import Scan
from tavi.instrument.components.collimators import Collimators
from tavi.instrument.components.component_base import Distances
from tavi.instrument.components.detector import Detector
from tavi.instrument.components.goni import Goniometer
from tavi.instrument.components.guide import Guide
from tavi.instrument.components.monitor import Monitor
from tavi.instrument.components.mono_ana import MonoAna
from tavi.instrument.components.source import Source
from tavi.sample import Sample


class TASBase(object):
    """
    Base class for the triple-axis spectrometer. Manage instrument congigutarion parameters.

    Attributes:
        source (Source): neutron source
        collimators (Collimators): 4 collimators
        guide (Guide): neutron guide
        monochromator (MonoAna): monochromator
        goniometer (Goniometer): goinometer
        analyzer (MonoAna): analyzer
        detector (Detector): detector
        arms (Distances): distances between components
        sample (Sample): sample being measured

    Methods:
        load_instrument_params_from_json(path_to_json)
    """

    def __init__(self) -> None:
        """Load instrument configuration from json if provided, otherwise leave as None"""
        self.monochromator: MonoAna
        self.analyzer: MonoAna
        self.goniometer: Goniometer
        self.collimators: Collimators
        self.sample: Sample

        self.source: Optional[Source] = None

        self.guide: Optional[Guide] = None
        self.monitor: Optional[Monitor] = None

        self.detector: Optional[Detector] = None
        self.arms: Optional[Distances] = None

        self.json_path: Optional[str] = None  # if loaded from json
        self.reference_scan: Optional[Scan] = None  # if loaded from scan

    def _load_instrument_parameters(self, config_params: dict):
        components = {
            "source": Source,
            "collimators": Collimators,
            "guide": Guide,
            "monochromator": MonoAna,
            "monitor": Monitor,
            "goniometer": Goniometer,
            "analyzer": MonoAna,
            "detector": Detector,
            "distances": Distances,
        }

        for component_name, component_class in components.items():
            param = config_params.get(component_name)
            if param is not None:
                setattr(self, component_name, component_class(param, component_name))

    def load_instrument_params_from_json(
        self,
        path_to_json: str,
    ):
        """Load instrument configuration from a json file"""

        json_file = Path(path_to_json)
        if json_file.is_file():
            with open(json_file, "r", encoding="utf-8") as file:
                config_params = json.load(file)

            self._load_instrument_parameters(config_params)
            self.json_path = path_to_json
        else:
            print("Invalid path for instrument configuration json file.")

    # TODO
    def load_insrument_params_from_scan(
        self,
        scan: Scan,
    ):
        """Load instrument configuration from a scan"""
        pass

    # TODO
    def save_instrument_params_to_json(
        self,
        path_to_json: str,
    ):
        """Save configuration into a dictionary"""
        # convert python dictionary to json file
        # with open("./src/tavi/instrument/instrument_params/takin_test.json", "w") as file:
        #     json.dump(instrument_config, file)
        pass

    def save_instrument_params_to_scan(
        self,
        scan: Scan,
        param: str,
        val: Union[str, float],
    ):
        """Save ONE instrument parameter to scan"""
        pass

    def mount_sample(self, sample: Sample) -> None:
        """Add sample info to the triple-axis spectrometer"""
        self.sample = sample
