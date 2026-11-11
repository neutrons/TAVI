"""Presenter for the combine panel."""

from tavi.backend.model.interface.combine_model_interface import CombineModelInterface
from tavi.frontend.presenter.abstract_presenter import AbstractPresenter
from tavi.frontend.view.combine_view import CombineView
from tavi.library.data.combine_request import CombineRequest
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RemoveRawScanEvent
from tavi.meta.event.type.presenter_event import ClearFocusEvent, FocusRawScanEvent


class CombinePresenter(AbstractPresenter):
    """Mediates between the raw scans staged in the combine panel and CombineModel."""

    def __init__(self, model: CombineModelInterface) -> None:
        """Create the view, subscribe to the focus events, and hook up the panel's buttons."""
        super().__init__()
        self._model = model
        # The focused raw scans' names, by uuid - what Add Scans stages. Names, not scans: a
        # presenter holds uuids and what its view has to render, never a domain object.
        self._focused_names: dict[UUID, str] = {}
        # What the list shows, in the order it shows it. The view works in names, the model in
        # uuids, so the pairing is kept here rather than being looked back up from either side.
        self._staged: dict[UUID, str] = {}

        self._event_broker = EventBroker()
        self._event_broker.register(ClearFocusEvent, self.handle_clear_focus)
        self._event_broker.register(FocusRawScanEvent, self.handle_raw_scan_focus)
        self._event_broker.register(RemoveRawScanEvent, self.handle_raw_scan_removed)

        self._view: CombineView
        self._view.add_scans_clicked.connect(self.handle_add_scans_clicked)
        self._view.remove_scans_clicked.connect(self.handle_scans_removed)
        self._view.append_clicked.connect(self.handle_append_clicked)
        self._view.combine_clicked.connect(self.handle_combine_clicked)

    def init_view(self) -> None:
        """Create the combine view."""
        self._view = CombineView()

    def handle_clear_focus(self, _: ClearFocusEvent) -> None:
        """Forget the old selection; the staged scans are the user's own list and outlive it."""
        self._focused_names = {}

    def handle_raw_scan_focus(self, e: FocusRawScanEvent) -> None:
        """Track which raw scans Add Scans would stage, and point the column fields at their columns."""
        self._focused_names.update({scan.uuid: scan.tavimeta.friendly_name for scan in e.scans})
        if not e.scans:
            return
        lead = e.scans[0]
        self._view.set_column_options_signal.emit(list(lead.data.data.keys()))
        # Seeded from the selection, as the plotter's axis fields are: a panel still naming the
        # last experiment's columns would combine into an error nothing in it explains.
        self._view.set_axes_signal.emit(*lead.tavimeta.default_axis)

    def handle_raw_scan_removed(self, e: RemoveRawScanEvent) -> None:
        """Unstage a scan that has left the project, so the list never names data nothing can resolve."""
        self._focused_names.pop(e.uuid, None)
        name = self._staged.pop(e.uuid, None)
        if name is not None:
            self._view.remove_scans([name])

    def handle_add_scans_clicked(self) -> None:
        """Stage every focused raw scan; the view lists the ones it did not already hold."""
        self._staged.update(self._focused_names)
        self._view.add_scans(list(self._focused_names.values()))

    def handle_scans_removed(self, names: list[str]) -> None:
        """Unstage the scans the user removed from the list, which has already dropped them."""
        removed = set(names)
        self._staged = {uuid: name for uuid, name in self._staged.items() if name not in removed}

    def handle_append_clicked(self) -> None:
        """Lay the staged scans end to end, with no rebin."""
        # An append applies no weighting - CombineOps refuses a normalization or a scale it has no
        # way to honour - so the fields that describe one are not sent with it.
        self._combine(RebinMode.NONE, normalization="", multiplication="1")

    def handle_combine_clicked(self) -> None:
        """Combine the staged scans with the rebin the panel's visible tab describes."""
        self._combine(
            self._view.rebin_mode(),
            normalization=self._view.normalization(),
            multiplication=self._view.multiplication(),
        )

    def _combine(self, mode: RebinMode, normalization: str, multiplication: str) -> None:
        """Ask the model for the combination the panel describes; it reports anything unusable itself."""
        self._model.combine(
            CombineRequest(
                source_scan_uuids=list(self._staged),
                x_axis=self._view.x_axis(),
                y_axis=self._view.y_axis(),
                normalization=normalization,
                multiplication=multiplication,
                rebin_mode=mode,
                tolerance=self._view.rebin_tolerance(),
                step=self._view.rebin_step(),
                start=self._view.rebin_start(),
                finish=self._view.rebin_stop(),
                z_axis=self._view.histogram_z_axis(),
                x_start=self._view.histogram_x_start(),
                x_stop=self._view.histogram_x_stop(),
                x_step=self._view.histogram_x_step(),
                y_start=self._view.histogram_y_start(),
                y_stop=self._view.histogram_y_stop(),
                y_step=self._view.histogram_y_step(),
            )
        )
