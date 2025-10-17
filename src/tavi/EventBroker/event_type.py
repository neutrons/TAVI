from attr import dataclass


class Event:
    pass


@dataclass
class scan_uuid(Event):
    scan_uuid_list: list[str]


@dataclass
class selected_uuid(Event):
    selected_uuid: str


@dataclass
class meta_data(Event):
    meta_data_dict: dict


@dataclass
class random_data(Event):
    random_data: str
