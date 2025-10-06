from tavi.Observer.load_observer import LoadObserver


class LoadPresenter:
    def __init__(self, view: any, model: any):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model
        self._view.connect_load_data(self.handle_load_data)
        self.load_observer = LoadObserver()
        # populate meta data

    def handle_load_data(self, data_dir_or_files):
        self._model.attach(self.load_observer)
        self._model.load(folder=data_dir_or_files[0])
        loaded_files = self.load_observer.get_loaded_data()
        print(loaded_files)
