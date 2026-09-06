"""
simplenmrium — NMRium/NOMAD client for simpleNMR.

    from simplenmrium import NMRiumData

    data = NMRiumData("exam_CMCse_1_NOMAD-pp.json")
    payload = data.createJsonDict()
"""

from .json_converter import NMRiumData
from .nmrium_reader import NMRDataset

__all__ = ["NMRiumData", "NMRDataset"]

__version__ = "0.1.0"
