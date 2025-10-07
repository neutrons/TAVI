from tavi.Observer.observer import Observer

class MetaDataPresenter(Observer):
    def __init__(self, view: any, model: any):
        super().__init__()
        """Constructor
        :view: hppt_view class type
        :model:hppt_model class type
        """
        self._view = view
        self._model = model
        
        # display meta data from clicked scan
        # self.metadata_observer = MetaDataObserver(self._view.metadata_widget.set_values)
        self._model.attach(self)

    def update(self, subject) -> None:
        self.selected_meta_data = subject.selected_metadata
        self._view.metadata_widget.set_values(self.selected_meta_data)