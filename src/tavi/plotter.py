import matplotlib.pyplot as plt


class Plot1DManager(object):
    """Manage a plot"""

    def __init__(self) -> None:
        self.title = None
        self.xlim = None
        self.ylim = None
        self.xlabel = None
        self.ylabel = None


class Plot2DManager(object):
    """Manage a plot"""

    def __init__(self) -> None:
        self.title = None
        self.xlim = None
        self.ylim = None
        self.zlim = None
        self.xlabel = None
        self.ylabel = None
