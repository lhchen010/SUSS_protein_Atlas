"""Runtime compatibility for the licensed DeepTMHMM 1.0 distribution.

DeepTMHMM requests a Matplotlib style name that was renamed in Matplotlib 3.6.
Python imports this module automatically when its directory is on PYTHONPATH.
"""

try:
    import matplotlib.style as _style

    _original_use = _style.use

    def _compatible_use(style):
        if style == "seaborn-whitegrid" and style not in _style.available:
            style = "seaborn-v0_8-whitegrid"
        return _original_use(style)

    _style.use = _compatible_use
except ImportError:
    pass
