from typing import Any
from tavi.Observer.observer import Observer

class LoadObserver(Observer):
    def __init__(self):
        self.loaded_data = None

    def update(self, subject: Any) -> None:
        self.loaded_data = subject.file_list
        
    def get_loaded_data(self):
        return self.loaded_data
        
        