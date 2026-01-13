"""Define event type here."""

from attr import dataclass


class Event:
    """Docstring for Event."""

    pass


@dataclass
class scan_uuid(Event):
    """scan uuid event."""

    scan_uuid_list: list[str]


@dataclass
class selected_uuid(Event):
    """selected data event."""

    selected_uuid: str


@dataclass
class meta_data(Event):
    """meta data event."""

    meta_data_dict: dict


@dataclass
class autoplot_data(Event):
    """template data event."""

    autoplot_data: list[list[float], list[float]]


@dataclass
class template_data(Event):
    """template data event."""

    template_data: str
