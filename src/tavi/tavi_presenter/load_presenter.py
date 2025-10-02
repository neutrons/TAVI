class LoadPresenter:
    def __init__(self, view: any, model:any):
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model
        self._view.connect_load_data(self.handle_load_data)

        # populate meta data
    
    def handle_load_data(self, data_dir_or_files):
        self._model.load_scans(data_folder=data_dir_or_files[0])