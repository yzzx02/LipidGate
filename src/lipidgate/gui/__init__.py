"""Keep data adapters and worker entry points independent of Qt DLL loading."""

__all__ = ["MainWindow", "main"]


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(name)
    from . import app
    return getattr(app, name)
