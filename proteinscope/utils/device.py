"""CUDA / device utilities shared by all PyTorch-dependent modules.

Usage::

    from utils.device import get_device, cuda_available, device_info

    device = get_device()          # "cuda" or "cpu"
    model = model.to(device)

All callers should use `get_device()` rather than repeating the
`torch.cuda.is_available()` check inline, so that the GPU selection
policy is controlled from one place.
"""

from __future__ import annotations

import logging

_log = logging.getLogger(__name__)
_device: str | None = None  # cached after first call


def get_device() -> str:
    """Return "cuda" if a CUDA-capable GPU is available, otherwise "cpu".

    The result is cached so `torch.cuda.is_available()` is called only once
    per process.  Safe to call from any thread after the first call.
    """
    global _device
    if _device is not None:
        return _device

    try:
        import torch
        if torch.cuda.is_available():
            gpu_name = torch.cuda.get_device_name(0)
            vram_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
            _log.info("CUDA available — using GPU: %s (%.1f GB VRAM)", gpu_name, vram_gb)
            _device = "cuda"
        else:
            _log.info("CUDA not available — running on CPU")
            _device = "cpu"
    except ImportError:
        _log.warning("torch not installed — defaulting to CPU")
        _device = "cpu"

    return _device


def cuda_available() -> bool:
    """Convenience wrapper: True if a CUDA GPU is available."""
    return get_device() == "cuda"


def device_info() -> dict:
    """Return a dict with device metadata for health-check endpoints."""
    info: dict = {"device": get_device(), "cuda_available": cuda_available()}
    if cuda_available():
        try:
            import torch
            props = torch.cuda.get_device_properties(0)
            info["gpu_name"] = props.name
            info["gpu_vram_gb"] = round(props.total_memory / 1e9, 2)
            info["cuda_version"] = torch.version.cuda
            info["device_count"] = torch.cuda.device_count()
        except Exception:
            pass
    return info
