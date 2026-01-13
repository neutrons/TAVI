"""Entry point."""

import asyncio
import os
import sys

from qtpy.QtWidgets import QApplication

from tavi import __version__ as tavi_version
from tavi.backend.model.interface.TaviProjectInterface import TaviProjectInterface
from tavi.configuration import Configuration
from tavi.frontend.main import start
from tavi.frontend.presenter.main_presenter import MainPresenter


def _print_text_splash():
    # TODO
    pass

def _createArgparser():
    import argparse

    parser = argparse.ArgumentParser(
        prog="TAVI", description="Triple Axis data Visualization Toolkit (TAVI) ", epilog="https://tavi.readthedocs.io/"
    )
    parser.add_argument("-v", "--version", action="version", version=tavi_version)
    parser.add_argument(
        "--headcheck",
        action="store_true",
        help="start the gui then shut it down after 5 seconds. This is used for testing",
    )
    return parser



def main(args=None):
    parser = _createArgparser()
    options, _ = parser.parse_known_args(args)

    # show the ascii splash screen
    _print_text_splash()
    
    return start(options)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
