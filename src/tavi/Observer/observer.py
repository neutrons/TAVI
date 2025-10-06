from abc import ABC, abstractmethod
from typing import Any


class Observer(ABC):
    @abstractmethod
    def update(self, subject: Any) -> None:
        """
        Receive update from subject.
        """
        pass
