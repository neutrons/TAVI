import ...

class MainPresenter:
    def __init__(self, model_dict):
       
        self.file_menu_presenter = FileMenuPresenter(self.exit, model)

        self.view = MainWindowView(self.file_menu_presenter.view,.....,) # sub qtwidgets init?
        self.view.addMenu(self.file_menu_presenter.view)

    def exit(self):
        # future implement of save dirty model
        pass
