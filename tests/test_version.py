from tavi import __version__


def test_version():
    print(__version__)
    assert __version__ == "unknown"
    # assert __version__ == "0.0.1"
