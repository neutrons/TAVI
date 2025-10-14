from attr import dataclass
from tavi.MessageBroker.message import Message

@dataclass
class scan_uuid(Message):
    scan_uuid_list: list[str]