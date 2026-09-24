"""TaviData module."""

from pydantic import BaseModel, ConfigDict, Field

from tavi.library.data.fit_entry import FitEntry
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, ProcessedScan, RawScan, Scan


class PurgeResult(BaseModel):
    """What a ``TaviData.purge`` actually deleted, per store, for callers that announce removals."""

    raw_scans: list[UUID] = Field(default_factory=list)
    plots: list[UUID] = Field(default_factory=list)
    fits: list[UUID] = Field(default_factory=list)

    model_config = ConfigDict(arbitrary_types_allowed=True)


class TaviData(BaseModel):
    """High level data object tracking loaded state."""

    raw_scans: dict[UUID, RawScan] = Field(default_factory=dict)
    plots: dict[UUID, Plot] = Field(default_factory=dict)
    processed_scans: dict[UUID, ProcessedScan] = Field(default_factory=dict)
    fits: dict[UUID, FitEntry] = Field(default_factory=dict)

    def all_scans(self) -> dict[UUID, Scan]:
        """Return every scan, raw and derived, in one lookup - what ``ProcessOps`` resolves uuids against."""
        return {**self.raw_scans, **self.processed_scans}

    def fetch_by_uuid(self, uuid: UUID) -> RawScan | ProcessedScan | Plot | FitEntry:
        """Return the scan, plot, or fit matching uuid, or raise KeyError."""
        for store in (self.raw_scans, self.processed_scans, self.plots, self.fits):
            if uuid in store:
                return store[uuid]

        raise KeyError(f"No such UUID {uuid} available in TaviData for any type.")

    def purge(self, uuids: list[UUID]) -> PurgeResult:
        """
        Delete the named items and everything left dangling by them, and report what went.

        A plot or fit holds no data of its own, only a ``source_scan_uuid`` per series, so a
        series whose scan is going can no longer be resolved and goes with it. A plot survives as
        long as it has a series left - removing one run should not destroy the rest of a fused,
        multi-series plot - whereas a ``FitEntry`` is bound to exactly one series and so has
        nothing to survive on. Unknown uuids are ignored: the tree may ask twice for the same
        item (e.g. a folder and a scan inside it both selected).

        Every store is mutated in place rather than rebound, because the models hold these same
        dicts by reference.
        """
        requested = list(dict.fromkeys(uuids))
        scans = [uuid for uuid in requested if uuid in self.raw_scans]
        plots = [uuid for uuid in requested if uuid in self.plots]
        fits = [uuid for uuid in requested if uuid in self.fits]

        orphaned = set(scans)
        pruned: dict[UUID, Plot] = {}
        for plot_uuid, plot in self.plots.items():
            if plot_uuid in plots:
                continue
            surviving = plot.without_scans(orphaned)
            if surviving is None:
                plots.append(plot_uuid)
            elif surviving is not plot:
                pruned[plot_uuid] = surviving

        fits.extend(
            fit_uuid
            for fit_uuid, fit in self.fits.items()
            if fit_uuid not in fits and fit.series.source_scan_uuid in orphaned
        )

        self.plots.update(pruned)
        for uuid in scans:
            del self.raw_scans[uuid]
        for uuid in plots:
            del self.plots[uuid]
        for uuid in fits:
            del self.fits[uuid]

        return PurgeResult(raw_scans=scans, plots=plots, fits=fits)
