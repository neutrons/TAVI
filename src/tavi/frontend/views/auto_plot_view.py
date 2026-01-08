"""Widgets for the main window."""
from typing import Optional

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy.QtCore import QObject
from qtpy.QtWidgets import (
    QVBoxLayout,
    QWidget,
)
class AutoPlotWidget(QWidget):
    """Widget that displays the plot."""

    def __init__(self, parent: Optional["QObject"] = None) -> None:
        """Constructor for the plotting widget.

        Args:
            parent (QObject): Optional parent

        """
        super().__init__(parent)
        layoutRight = QVBoxLayout()

        self.figure = Figure(figsize=(6, 4.5))
        self.static_canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.static_canvas, self)

        layoutRight.addWidget(self.static_canvas)
        layoutRight.addWidget(self.toolbar)
        self.setLayout(layoutRight)

        # heatmap initialization
        self.ax = self.static_canvas.figure.subplots()

        # draw the plot
        self.static_canvas.draw()

    def update_plot(self, autoplot_data: list[list[float], list[float]]) -> None:
        """Update the plot."""
        # clear
        self.ax.clear()
        x, y = autoplot_data
        # update heatmap
        self.ax.plot(x, y, '.')

        self.static_canvas.draw()
