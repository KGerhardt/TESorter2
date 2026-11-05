try:
    from importlib.metadata import version as _version
    __version__ = _version("tesorter2")
except Exception:  # running from a source tree that was never installed
    __version__ = "unknown"
