import logging
import os
from typing import Iterable, Optional

from tavi.tavi_model.FileSystem.Facilities.load_ansto import LoadANSTO
from tavi.tavi_model.FileSystem.Facilities.load_frmii import LoadFRMII
from tavi.tavi_model.FileSystem.Facilities.load_hzb import LoadHZB
from tavi.tavi_model.FileSystem.Facilities.load_nist import LoadNIST
from tavi.tavi_model.FileSystem.Facilities.load_ornl import LoadORNL
from tavi.tavi_model.FileSystem.Facilities.load_psi import LoadPSI
from tavi.tavi_model.FileSystem.tavi_class_factory import Scan

logger = logging.getLogger("TAVI")


class LoadManager:
    """
    TAVI file management system.

    - Normalizes user input (folder and files) and resolves paths.
    - Ranks likely facility using facility-provided scoring functions.
    - Loads data through the chosen facility loader.

    Parameters
    ----------
    data_folder : str | os.PathLike | None
        Folder path containing data files (optional if `data_file` is provided).
    data_file : str | os.PathLike | Iterable[str | os.PathLike] | None
        A single file or an iterable of files.
    facility : str | None
        If provided, bypass ranking and force this facility (must be supported).
    """

    # currently supported facilities
    supported_facilities = ["ORNL", "ANSTO", "ILL", "PSI", "NIST", "HZB", "FRMII"]

    def __init__(
        self,
        data_folder: Optional[os.PathLike | str] = None,
        data_files: Optional[os.PathLike | str | Iterable[os.PathLike | str]] = None,
        ub_dir: Optional[os.PathLike] = None,
        facility: Optional[str] = None,
    ) -> None:
        self.data_folder = data_folder
        self.data_files = data_files
        self.ub_dir = ub_dir
        self.facility = facility

    # TO DO
    def rank_facility(self) -> None:
        """
        Determine and set the most likely facility using plugin score functions.

        Returns
        -------
        str
            The selected facility name.
        """
        facility_score = dict()
        for facility in self.supported_facilities:
            match facility:
                case "ORNL":
                    facility_score[facility] = LoadORNL(
                        data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir
                    ).score()
                case "ANSTO":
                    facility_score[facility] = LoadANSTO(
                        data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir
                    ).score()
                case "ILL":
                    facility_score[facility] = LoadANSTO(
                        data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir
                    ).score()
                case "PSI":
                    facility_score[facility] = LoadPSI(
                        data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir
                    ).score()
                case "NIST":
                    facility_score[facility] = LoadNIST(
                        data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir
                    ).score()
                case "HZB":
                    facility_score[facility] = LoadHZB(
                        data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir
                    ).score()
                case "FRMII":
                    facility_score[facility] = LoadFRMII(
                        data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir
                    ).score()
        self.facility = max(facility_score, key=facility_score.get)

    def load(self) -> Scan:
        """
        load the either a folder or a single file or a list of files and returns the TaviProject class
        """
        # if self.facility is None, call the ranking system
        if not self.facility:
            self.rank_facility()
        match self.facility:
            case "ORNL":
                return LoadORNL(data_folder=self.data_folder, data_files=self.data_files, ub_dir=self.ub_dir).load()

            # TO DO: extend to other facilities
            # case "ILL":
            # return load_ill.load(self.data_folder, self.data_files)
            # case "ANSTO":
            # return load_ansto.load(self.data_folder, self.data_files)
            # case "PSI":
            # return load_psi.load(self.data_folder, self.data_files)
            # case "NIST":
            # return load_nist.load(self.data_folder, self.data_files)
            # case "MLZ":
            # return load_ill.load(self.data_folder, self.data_files)
            # case "FRMII":
            # return load_frmii.load(self.data_folder, self.data_files)
