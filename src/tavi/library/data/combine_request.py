"""What a user asked a combination to be, before anything is parsed or run."""

from pydantic import BaseModel, ConfigDict

from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.scan import UUID


class CombineRequest(BaseModel):
    """
    What the user asked a combination to be, as they expressed it.

    The numbers are the text they were typed as, and the model that runs this parses them, so a
    half-typed field is reported as the error it is rather than refused at the request's door.
    ``rebin_mode`` names which of them the combination actually reads: a tolerance, or a step and
    the grid edges, or a 2D grid over ``x_axis`` and ``z_axis`` for a histogram, or - for a plain
    append - none of them.
    """

    source_scan_uuids: list[UUID]
    x_axis: str
    y_axis: str
    # The column the rebins weigh by - the monitor or time the points were counted against.
    # Blank is no normalization, which only an append may ask for.
    normalization: str = ""
    multiplication: str = "1"
    rebin_mode: RebinMode = RebinMode.NONE
    tolerance: str = ""
    step: str = ""
    # Blank edges leave the grid to the data, which is what ``CombineOps`` reads ``None`` as.
    start: str = ""
    finish: str = ""
    # The histogram's second map axis and its grid, one Start/Stop/Step per map axis; blank
    # edges are read off the data, as for the 1D equal-step grid.
    z_axis: str = ""
    x_start: str = ""
    x_stop: str = ""
    x_step: str = ""
    y_start: str = ""
    y_stop: str = ""
    y_step: str = ""

    model_config = ConfigDict(arbitrary_types_allowed=True)
