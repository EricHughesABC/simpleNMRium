"""
nmrium_reader.py

Parses NMRium/NOMAD export JSON files into a clean object model. This is
Eric's original nmr_analysis.py, kept as the single source of truth for
"how do I read an NMRium file" logic, with two changes from the original:

1. RDKit/atom-info/SMILES/molfile-envelope construction has been split OUT
   of NMRDataset and moved to json_converter.py. This module's job is
   purely "parse the NMRium JSON shape into objects" — building the
   simpleNMR contract's allAtomsInfo/carbonAtomsInfo envelopes is a
   different concern (contract-shape knowledge, not NMRium-shape
   knowledge) and belongs next to the rest of the SimpleNMRBuilder wiring
   instead. NMRDataset here still exposes .molecules / .smiles / .molfile
   as plain parsed values so json_converter.py has something to build
   from.
2. Range1D.all_peaks previously filtered on `signal.kind == 'signal'`;
   Zone2D.all_peaks did not have the equivalent filter despite the
   docstring claiming it did. Both now consistently skip non-'signal'
   kind entries (NMRium can carry other signal "kind" values, e.g.
   reference peaks) so 1D and 2D peak counts are derived the same way.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Union, Tuple
from pathlib import Path
from enum import Enum
import json


class Dimension(Enum):
    """Spectrum dimensionality"""
    ONE_D = 1
    TWO_D = 2


# ============================================================================
# 1D Spectrum Classes
# ============================================================================

@dataclass
class Peak1D:
    """Represents a single peak in a 1D spectrum"""

    id: str
    x: float  # Chemical shift (ppm)
    original_x: float
    y: float  # Intensity/height
    width: float
    fwhm: Optional[float] = None
    shape_kind: Optional[str] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Peak1D":
        shape = data.get("shape", {})
        return cls(
            id=data.get("id", ""),
            x=data.get("x", 0.0),
            original_x=data.get("originalX", data.get("x", 0.0)),
            y=data.get("y", 0.0),
            width=data.get("width", 0.0),
            fwhm=shape.get("fwhm"),
            shape_kind=shape.get("kind"),
        )

    def __repr__(self) -> str:
        return f"Peak1D(δ={self.x:.3f} ppm, height={self.y:.1f})"


@dataclass
class Signal1D:
    """Represents a signal (multiplet) in a 1D spectrum"""

    id: str
    delta: float
    original_delta: float
    multiplicity: Optional[str] = None
    peaks: List[Peak1D] = field(default_factory=list)
    coupling_constants: List[Dict[str, Any]] = field(default_factory=list)
    assignment: Optional[str] = None
    kind: str = "signal"

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Signal1D":
        peaks = [Peak1D.from_dict(p) for p in data.get("peaks", [])]
        return cls(
            id=data.get("id", ""),
            delta=data.get("delta", 0.0),
            original_delta=data.get("originalDelta", data.get("delta", 0.0)),
            multiplicity=data.get("multiplicity"),
            peaks=peaks,
            coupling_constants=data.get("js", []),
            assignment=data.get("assignment"),
            kind=data.get("kind", "signal"),
        )

    @property
    def n_peaks(self) -> int:
        return len(self.peaks)

    @property
    def has_peaks(self) -> bool:
        return bool(self.peaks)

    @property
    def max_intensity(self) -> float:
        return max((p.y for p in self.peaks), default=0.0)

    def __repr__(self) -> str:
        mult = f", {self.multiplicity}" if self.multiplicity else ""
        return f"Signal1D(δ={self.delta:.3f} ppm{mult}, {self.n_peaks} peaks)"


@dataclass
class Range1D:
    """Represents a range (integrated region) in a 1D spectrum"""

    id: str
    from_ppm: float
    to_ppm: float
    original_from: float
    original_to: float
    integration: float
    absolute: float
    signals: List[Signal1D] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Range1D":
        signals = [Signal1D.from_dict(s) for s in data.get("signals", [])]
        return cls(
            id=data.get("id", ""),
            from_ppm=data.get("from", 0.0),
            to_ppm=data.get("to", 0.0),
            original_from=data.get("originalFrom", data.get("from", 0.0)),
            original_to=data.get("originalTo", data.get("to", 0.0)),
            integration=data.get("integration", 0.0),
            absolute=data.get("absolute", 0.0),
            signals=signals,
        )

    @property
    def width(self) -> float:
        return abs(self.to_ppm - self.from_ppm)

    @property
    def center(self) -> float:
        return (self.from_ppm + self.to_ppm) / 2

    @property
    def n_signals(self) -> int:
        return len(self.signals)

    @property
    def all_peaks(self) -> List[Peak1D]:
        return [peak for signal in self.signals if signal.kind == "signal" for peak in signal.peaks]

    def __repr__(self) -> str:
        return f"Range1D({self.from_ppm:.3f}-{self.to_ppm:.3f} ppm, integral={self.integration:.2f}, {self.n_signals} signals)"


# ============================================================================
# 2D Spectrum Classes
# ============================================================================

@dataclass
class Peak2D:
    """Represents a single peak in a 2D spectrum"""

    id: str
    x: float
    y: float
    z: float
    original_x: float
    original_y: float
    min_x: Optional[float] = None
    max_x: Optional[float] = None
    min_y: Optional[float] = None
    max_y: Optional[float] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Peak2D":
        return cls(
            id=data.get("id", ""),
            x=data.get("x", 0.0),
            y=data.get("y", 0.0),
            z=data.get("z", 0.0),
            original_x=data.get("originalX", data.get("x", 0.0)),
            original_y=data.get("originalY", data.get("y", 0.0)),
            min_x=data.get("minX"),
            max_x=data.get("maxX"),
            min_y=data.get("minY"),
            max_y=data.get("maxY"),
        )

    def __repr__(self) -> str:
        return f"Peak2D(F2={self.x:.3f}, F1={self.y:.3f} ppm, intensity={self.z:.1e})"


@dataclass
class Axis2D:
    """Axis information for a 2D signal"""

    delta: float
    original_delta: float
    nucleus: str
    resolution: Optional[float] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Axis2D":
        return cls(
            delta=data.get("delta", 0.0),
            original_delta=data.get("originalDelta", data.get("delta", 0.0)),
            nucleus=data.get("nucleus", ""),
            resolution=data.get("resolution"),
        )


@dataclass
class Signal2D:
    """Represents a correlation signal in a 2D spectrum"""

    id: str
    x_axis: Axis2D
    y_axis: Axis2D
    peaks: List[Peak2D] = field(default_factory=list)
    kind: str = "signal"

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Signal2D":
        peaks = [Peak2D.from_dict(p) for p in data.get("peaks", [])]
        x_data = data.get("x", {})
        y_data = data.get("y", {})
        return cls(
            id=data.get("id", ""),
            x_axis=Axis2D.from_dict(x_data),
            y_axis=Axis2D.from_dict(y_data),
            peaks=peaks,
            kind=data.get("kind", "signal"),
        )

    @property
    def n_peaks(self) -> int:
        return len(self.peaks)

    @property
    def max_intensity(self) -> float:
        return max((p.z for p in self.peaks), default=0.0)

    @property
    def x_shift(self) -> float:
        return self.x_axis.delta

    @property
    def y_shift(self) -> float:
        return self.y_axis.delta

    @property
    def average_x(self) -> float:
        if not self.peaks:
            return 0.0
        return sum(p.x for p in self.peaks) / len(self.peaks)

    @property
    def average_y(self) -> float:
        if not self.peaks:
            return 0.0
        return sum(p.y for p in self.peaks) / len(self.peaks)

    @property
    def average_intensity(self) -> float:
        if not self.peaks:
            return 0.0
        return sum(p.z for p in self.peaks) / len(self.peaks)

    @property
    def total_intensity(self) -> float:
        return sum(p.z for p in self.peaks)

    def __repr__(self) -> str:
        return f"Signal2D(F2={self.x_shift:.3f}, F1={self.y_shift:.3f} ppm, {self.n_peaks} peaks)"


@dataclass
class Zone2D:
    """Represents a zone (region) in a 2D spectrum"""

    id: str
    x_from: float
    x_to: float
    y_from: float
    y_to: float
    signals: List[Signal2D] = field(default_factory=list)
    kind: str = "signal"

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Zone2D":
        signals = [Signal2D.from_dict(s) for s in data.get("signals", [])]
        x_data = data.get("x", {})
        y_data = data.get("y", {})
        return cls(
            id=data.get("id", ""),
            x_from=x_data.get("from", 0.0),
            x_to=x_data.get("to", 0.0),
            y_from=y_data.get("from", 0.0),
            y_to=y_data.get("to", 0.0),
            signals=signals,
            kind=data.get("kind", "signal"),
        )

    @property
    def x_width(self) -> float:
        return abs(self.x_to - self.x_from)

    @property
    def y_width(self) -> float:
        return abs(self.y_to - self.y_from)

    @property
    def x_center(self) -> float:
        return (self.x_from + self.x_to) / 2

    @property
    def y_center(self) -> float:
        return (self.y_from + self.y_to) / 2

    @property
    def n_signals(self) -> int:
        return len(self.signals)

    @property
    def all_peaks(self) -> List[Peak2D]:
        return [peak for signal in self.signals if signal.kind == "signal" for peak in signal.peaks]

    @property
    def average_x(self) -> float:
        peaks = self.all_peaks
        if not peaks:
            return 0.0
        return sum(p.x for p in peaks) / len(peaks)

    @property
    def average_y(self) -> float:
        peaks = self.all_peaks
        if not peaks:
            return 0.0
        return sum(p.y for p in peaks) / len(peaks)

    @property
    def average_intensity(self) -> float:
        peaks = self.all_peaks
        if not peaks:
            return 0.0
        return sum(p.z for p in peaks) / len(peaks)

    @property
    def total_intensity(self) -> float:
        return sum(p.z for p in self.all_peaks)

    @property
    def max_intensity(self) -> float:
        peaks = self.all_peaks
        return max((p.z for p in peaks), default=0.0)

    def __repr__(self) -> str:
        return f"Zone2D(F2: {self.x_from:.3f}-{self.x_to:.3f}, F1: {self.y_from:.3f}-{self.y_to:.3f} ppm, {self.n_signals} signals)"


# ============================================================================
# Spectrum
# ============================================================================

@dataclass
class Spectrum:
    """Represents a single NMR spectrum (1D or 2D)"""

    nucleus: List[str]
    dimension: int
    experiment: str
    name: str
    exp_id: str
    experiment_number: int

    solvent: str
    temperature: List[float]
    number_of_scans: List[int]
    pulse_sequence: List[str]

    origin_frequency: List[float]
    base_frequency: List[float]
    spectral_width: List[float]
    spectrum_size: List[int]
    number_of_points: List[int]

    phc0: List[float]
    phc1: List[float]

    field_strength: float
    acquisition_time: float
    relaxation_time: float

    ranges: List[Range1D] = field(default_factory=list)
    zones: List[Zone2D] = field(default_factory=list)

    owner: Optional[str] = None
    title: Optional[str] = None
    creator: Optional[str] = None
    probe_name: Optional[List[str]] = None
    noise: Optional[Dict[str, float]] = None

    raw_info: Optional[Dict[str, Any]] = None

    @classmethod
    def from_dict(cls, spectrum_data: Dict[str, Any]) -> "Spectrum":
        info = spectrum_data.get("info", {})

        spectrum = cls(
            nucleus=info.get("nucleus", []) if isinstance(info.get("nucleus"), list) else [info.get("nucleus", "")],
            dimension=info.get("dimension", 1),
            experiment=info.get("experiment", "unknown"),
            name=info.get("name", "undefined"),
            exp_id=info.get("expId", ""),
            experiment_number=info.get("experimentNumber", 0),
            solvent=info.get("solvent", ""),
            temperature=info.get("temperature", []) if isinstance(info.get("temperature"), list) else [info.get("temperature", 0.0)],
            number_of_scans=info.get("numberOfScans", []) if isinstance(info.get("numberOfScans"), list) else [info.get("numberOfScans", 0)],
            pulse_sequence=info.get("pulseSequence", []),
            origin_frequency=info.get("originFrequency", []),
            base_frequency=info.get("baseFrequency", []),
            spectral_width=info.get("spectralWidth", []),
            spectrum_size=info.get("spectrumSize", []),
            number_of_points=info.get("numberOfPoints", []),
            phc0=info.get("phc0", []),
            phc1=info.get("phc1", []),
            field_strength=info.get("fieldStrength", 0.0),
            acquisition_time=info.get("acquisitionTime", 0.0),
            relaxation_time=info.get("relaxationTime", 0.0),
            owner=info.get("owner"),
            title=info.get("title"),
            creator=info.get("creator"),
            probe_name=info.get("probeName"),
            noise=info.get("noise"),
            raw_info=info,
        )

        if "ranges" in spectrum_data and spectrum_data["ranges"]:
            ranges_data = spectrum_data["ranges"].get("values", [])
            spectrum.ranges = [Range1D.from_dict(r) for r in ranges_data]

        if "zones" in spectrum_data and spectrum_data["zones"]:
            zones_data = spectrum_data["zones"].get("values", [])
            spectrum.zones = [Zone2D.from_dict(z) for z in zones_data]

        return spectrum

    @property
    def is_1d(self) -> bool:
        return self.dimension == 1

    @property
    def is_2d(self) -> bool:
        return self.dimension == 2

    @property
    def primary_nucleus(self) -> str:
        return self.nucleus[0] if self.nucleus else "unknown"

    @property
    def secondary_nucleus(self) -> Optional[str]:
        return self.nucleus[1] if len(self.nucleus) > 1 else None

    @property
    def n_ranges(self) -> int:
        return len(self.ranges)

    @property
    def n_zones(self) -> int:
        return len(self.zones)

    @property
    def total_peaks_1d(self) -> int:
        return sum(len(r.all_peaks) for r in self.ranges)

    @property
    def total_peaks_2d(self) -> int:
        return sum(len(z.all_peaks) for z in self.zones)

    def get_all_peaks_1d(self) -> List[Tuple[float, float]]:
        peaks = []
        for range_obj in self.ranges:
            for peak in range_obj.all_peaks:
                peaks.append((peak.x, peak.y))
        return peaks

    def get_all_peaks_2d(self) -> List[Tuple[float, float, float]]:
        peaks = []
        for zone in self.zones:
            for peak in zone.all_peaks:
                peaks.append((peak.x, peak.y, peak.z))
        return peaks

    def get_zone_averages(self) -> List[Tuple[float, float, float, float]]:
        """Returns [(avg_x, avg_y, avg_intensity, total_intensity), ...] per zone."""
        averages = []
        for zone in self.zones:
            averages.append((zone.average_x, zone.average_y, zone.average_intensity, zone.total_intensity))
        return averages

    def get_all_integrals(self) -> List[Tuple[float, float, float]]:
        integrals = []
        for range_obj in self.ranges:
            integrals.append((range_obj.center, range_obj.width, range_obj.integration))
        return integrals

    def summary(self) -> str:
        lines = [
            f"Spectrum: {self.experiment.upper()}",
            f"  Dimension: {self.dimension}D",
            f"  Nucleus: {' / '.join(self.nucleus)}",
            f"  Experiment #{self.experiment_number}",
            f"  Pulse sequence: {self.pulse_sequence}",
            f"  Solvent: {self.solvent}",
            f"  Field: {self.field_strength:.2f} T",
        ]
        if self.is_1d:
            lines.append(f"  Ranges: {self.n_ranges}")
            lines.append(f"  Total peaks: {self.total_peaks_1d}")
        else:
            lines.append(f"  Zones: {self.n_zones}")
            lines.append(f"  Total peaks: {self.total_peaks_2d}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"Spectrum(exp={self.experiment}, dim={self.dimension}D, nucleus={self.nucleus})"


# ============================================================================
# Dataset
# ============================================================================

class NMRDataset:
    """
    Container for all spectra + the molecule(s) in an NMRium/NOMAD file.

    Deliberately does NOT touch RDKit or build simpleNMR contract envelopes
    (allAtomsInfo, carbonAtomsInfo, smiles/molfile scalar wrapping) — that's
    json_converter.py's job. This class exposes the raw parsed values
    (.molecules, .smiles, .correlations) so the converter has something to
    build from without re-parsing the file itself.
    """

    def __init__(self, nomad_dict: Dict[str, Any]):
        self.raw_data = nomad_dict
        self.spectra: List[Spectrum] = []
        self.source: Optional[str] = None

        data = nomad_dict["nmriumData"]["data"]
        self.molecules: List[Dict[str, Any]] = data.get("molecules", [])
        self.smiles: List[str] = nomad_dict.get("smiles", [])
        self.correlations: Dict[str, Any] = data.get("correlations", {})
        self.title: Optional[str] = nomad_dict.get("title")

        self._parse_spectra()

    @classmethod
    def from_file(cls, filename: Union[str, Path]) -> "NMRDataset":
        filepath = Path(filename)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        if not filepath.is_file():
            raise ValueError(f"Path is not a file: {filepath}")
        with open(filepath, "r", encoding="utf-8") as f:
            nomad_dict = json.load(f)
        instance = cls(nomad_dict)
        instance.source = str(filepath.absolute())
        return instance

    @classmethod
    def from_fileobj(cls, fileobj) -> "NMRDataset":
        try:
            nomad_dict = json.load(fileobj)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in file object: {e}")
        instance = cls(nomad_dict)
        if hasattr(fileobj, "name"):
            instance.source = fileobj.name
        return instance

    @classmethod
    def from_json_string(cls, json_string: str) -> "NMRDataset":
        try:
            nomad_dict = json.loads(json_string)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON string: {e}")
        return cls(nomad_dict)

    def _parse_spectra(self):
        try:
            spectra_list = self.raw_data["nmriumData"]["data"]["spectra"]
            for spectrum_data in spectra_list:
                self.spectra.append(Spectrum.from_dict(spectrum_data))
        except KeyError as e:
            raise ValueError(f"Invalid NMRium structure: missing key {e}")
        except TypeError as e:
            raise ValueError(f"Invalid NMRium structure: {e}")

    def get_1d_spectra(self) -> List[Spectrum]:
        return [s for s in self.spectra if s.is_1d]

    def get_2d_spectra(self) -> List[Spectrum]:
        return [s for s in self.spectra if s.is_2d]

    def get_by_experiment(self, experiment_type: str) -> List[Spectrum]:
        return [s for s in self.spectra if s.experiment.lower() == experiment_type.lower()]

    def get_by_nucleus(self, nucleus: str) -> List[Spectrum]:
        return [s for s in self.spectra if nucleus in s.nucleus]

    def get_spectrum(self, index: int) -> Spectrum:
        return self.spectra[index]

    @property
    def n_spectra(self) -> int:
        return len(self.spectra)

    @property
    def primary_molfile(self) -> Optional[str]:
        """The first molecule's raw V3000/V2000 molblock text, or None."""
        if not self.molecules:
            return None
        return self.molecules[0].get("molfile")

    @property
    def primary_smiles(self) -> Optional[str]:
        return self.smiles[0] if self.smiles else None

    def summary(self) -> str:
        lines = []
        if self.title:
            lines.append(f"NMR Dataset: {self.title}")
        elif self.source:
            lines.append(f"NMR Dataset from: {Path(self.source).name}")
        else:
            lines.append("NMR Dataset")

        lines.append(f"Total spectra: {self.n_spectra}")
        lines.append("")

        for i, spectrum in enumerate(self.spectra):
            lines.append(f"[{i}] {spectrum.experiment.upper()} ({spectrum.dimension}D)")
            lines.append(f"    Nucleus: {' / '.join(spectrum.nucleus)}")
            if spectrum.is_1d:
                lines.append(f"    Ranges: {spectrum.n_ranges}, Peaks: {spectrum.total_peaks_1d}")
            else:
                lines.append(f"    Zones: {spectrum.n_zones}, Peaks: {spectrum.total_peaks_2d}")

        lines.append("")
        lines.append(f"smiles: {self.smiles}")

        return "\n".join(lines)

    def detailed_summary(self) -> str:
        lines = [self.summary(), "", "=" * 60, ""]
        for i, spectrum in enumerate(self.spectra):
            lines.append(f"\nSpectrum [{i}]: {spectrum.experiment.upper()}")
            lines.append("-" * 60)
            lines.append(spectrum.summary())

            if spectrum.is_1d and spectrum.ranges:
                lines.append("\n  Ranges:")
                for j, range_obj in enumerate(spectrum.ranges[:5]):
                    lines.append(f"    [{j}] {range_obj}")
                if len(spectrum.ranges) > 5:
                    lines.append(f"    ... and {len(spectrum.ranges) - 5} more")

            elif spectrum.is_2d and spectrum.zones:
                lines.append("\n  Zones:")
                for j, zone in enumerate(spectrum.zones[:5]):
                    lines.append(f"    [{j}] {zone}")
                if len(spectrum.zones) > 5:
                    lines.append(f"    ... and {len(spectrum.zones) - 5} more")

        return "\n".join(lines)

    def to_file(self, filename: Union[str, Path], indent: int = 2):
        filepath = Path(filename)
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(self.raw_data, f, indent=indent)

    def __repr__(self) -> str:
        source_str = f" from {Path(self.source).name}" if self.source else ""
        return f"NMRDataset({self.n_spectra} spectra{source_str})"

    def __len__(self) -> int:
        return self.n_spectra

    def __getitem__(self, index: int) -> Spectrum:
        return self.spectra[index]
