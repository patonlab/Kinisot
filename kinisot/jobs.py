"""JSON job files: everything flags cannot carry (series, channels, several isotopologues).

    kinisot --job job.json

A job file holds the structures, the labels and the settings; the command
line then gives only the outputs (-o, --json, --csv, -q, --overwrite).
File names are relative to the job file. docs/job_files.md describes the
format with examples.
"""

import json
import os

from .api import compute_kie
from .ensemble import Conformers
from .exceptions import KinisotInputError
from .pathways import Series, channels

__all__ = ["load_job", "run_job", "SETTINGS"]

# job-file key: compute_kie keyword
SETTINGS = {
    "scale": "scale",
    "scale_type": "scale_type",
    "imag_cutoff": "imag_cutoff",
    "tunneling": "tunneling",
    "barrier": "barrier",
    "project": "project",
    "calc": "calculator",
    "delta": "delta",
    "weights": "weights",
    "weight_uncertainty": "weight_uncertainty",
}
STRUCTURES = ("transition_structure", "product", "series", "channels")
TOP_KEYS = set(SETTINGS) | set(STRUCTURES) | {
    "temperature", "reactants", "iso", "reference", "isotopologues", "shares", "barriers", "amounts", "energy_unit",
    "description",
}  # fmt: skip
CHANNEL_KEYS = {"name", "reactants", "transition_structure", "series", "amount"}
SERIES_KEYS = {"steps", "free_energies", "commitment", "energy_unit"}
SPECIES_KEYS = {"files", "free_energies", "degeneracy", "energy_unit"}


def _check_keys(where, data, allowed):
    if not isinstance(data, dict):
        raise KinisotInputError("%s must be a JSON object" % where)
    unknown = sorted(set(data) - set(allowed))
    if unknown:
        raise KinisotInputError(
            "%s: unknown key(s) %s (allowed: %s)" % (where, ", ".join(unknown), ", ".join(sorted(allowed)))
        )


def _path(base, name):
    if not isinstance(name, str):
        raise KinisotInputError("expected a file name, got %r" % (name,))
    return name if os.path.isabs(name) else os.path.normpath(os.path.join(base, name))


def _species(base, value, where):
    """A file name, a list of conformer files, or {"files": [...], "free_energies": [...], "degeneracy": [...]}."""
    if isinstance(value, str):
        return _path(base, value)
    if isinstance(value, list):
        return [_path(base, f) for f in value]
    _check_keys(where, value, SPECIES_KEYS)
    if "files" not in value:
        raise KinisotInputError('%s needs "files"' % where)
    return Conformers(
        [_path(base, f) for f in value["files"]],
        free_energies=value.get("free_energies"),
        degeneracy=value.get("degeneracy"),
        energy_unit=value.get("energy_unit", "kcal/mol"),
    )


def _side(base, value, where):
    if not isinstance(value, list) or not value:
        raise KinisotInputError("%s must be a non-empty list of species" % where)
    return [_species(base, v, "%s[%d]" % (where, i)) for i, v in enumerate(value)]


def _series(base, value, where):
    _check_keys(where, value, SERIES_KEYS)
    if "steps" not in value:
        raise KinisotInputError('%s needs "steps"' % where)
    steps = [_species(base, v, "%s.steps[%d]" % (where, i)) for i, v in enumerate(value["steps"])]
    return Series(
        steps,
        free_energies=value.get("free_energies"),
        commitment=value.get("commitment"),
        energy_unit=value.get("energy_unit", "kcal/mol"),
    )


def load_job(path):
    """Read and check a job file. Returns a dictionary run_job() takes."""
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except OSError as err:
        raise KinisotInputError("cannot read the job file %s: %s" % (path, err.strerror or err)) from None
    except json.JSONDecodeError as err:
        raise KinisotInputError("the job file %s is not valid JSON: %s" % (path, err)) from None
    _check_keys("the job file", data, TOP_KEYS)
    base = os.path.dirname(os.path.abspath(path))
    given = [k for k in STRUCTURES if k in data]
    if len(given) != 1:
        raise KinisotInputError("a job file needs exactly one of %s" % ", ".join('"%s"' % k for k in STRUCTURES))
    kind = given[0]
    job = {"kind": kind, "settings": {SETTINGS[k]: data[k] for k in SETTINGS if k in data}}
    job["temperatures"] = _temperatures(data.get("temperature", 298.15))
    if kind == "channels":
        if "reactants" in data:
            raise KinisotInputError('with "channels", each channel gives its own "reactants"')
        job["channels"] = []
        for c, channel in enumerate(data["channels"]):
            where = "channels[%d]" % c
            _check_keys(where, channel, CHANNEL_KEYS)
            if ("transition_structure" in channel) == ("series" in channel) or "reactants" not in channel:
                raise KinisotInputError('%s needs "reactants" and one of "transition_structure" or "series"' % where)
            ts = (
                _series(base, channel["series"], where + ".series")
                if "series" in channel
                else _side(base, channel["transition_structure"], where + ".transition_structure")
            )
            job["channels"].append(
                {
                    "name": channel.get("name", "channel %d" % (c + 1)),
                    "rct": _side(base, channel["reactants"], where + ".reactants"),
                    "ts": ts,
                    "amount": channel.get("amount", 1.0),
                }
            )
        for key in ("shares", "barriers", "amounts", "energy_unit"):
            if key in data:
                job[key] = data[key]
    else:
        for key in ("shares", "barriers", "amounts"):
            if key in data:
                raise KinisotInputError('"%s" belongs to a job with "channels"' % key)
        if "reactants" not in data:
            raise KinisotInputError('the job file needs "reactants"')
        job["rct"] = _side(base, data["reactants"], "reactants")
        if kind == "series":
            job["ts"] = _series(base, data["series"], "series")
        elif kind == "transition_structure":
            job["ts"] = _side(base, data["transition_structure"], "transition_structure")
        else:
            job["prd"] = _side(base, data["product"], "product")
    if ("isotopologues" in data) == ("iso" in data):
        raise KinisotInputError('give "iso" (one isotopologue) or "isotopologues" (a list), one of them')
    if "isotopologues" in data:
        entries = data["isotopologues"]
        if not isinstance(entries, list) or not entries:
            raise KinisotInputError('"isotopologues" must be a non-empty list')
    else:
        entries = [{"iso": data["iso"], "reference": data.get("reference")}]
    if "isotopologues" in data and "reference" in data:
        raise KinisotInputError('with "isotopologues", give each one its own "reference"')
    job["isotopologues"] = []
    for i, entry in enumerate(entries):
        _check_keys("isotopologues[%d]" % i, entry, {"name", "iso", "reference"})
        if entry.get("iso") is None:
            raise KinisotInputError('isotopologues[%d] needs "iso"' % i)
        if kind == "channels":
            for key in ("iso", "reference"):
                value = entry.get(key)
                if value is not None and (not isinstance(value, list) or len(value) != len(job["channels"])):
                    raise KinisotInputError(
                        'with channels, "%s" is a list with one entry (that channel\'s labels) per channel' % key
                    )
        name = entry.get("name") or (
            " | ".join(_label_text(x) for x in entry["iso"]) if kind == "channels" else _label_text(entry["iso"])
        )
        job["isotopologues"].append({"name": str(name), "iso": entry["iso"], "reference": entry.get("reference")})
    return job


def _temperatures(value):
    from .cli import parse_temperatures  # imported here: the command line imports this module

    try:
        values = [float(t) for t in value] if isinstance(value, list) else parse_temperatures(str(value))
    except (TypeError, ValueError) as err:
        raise KinisotInputError("the job file's temperature: %s" % err) from None
    if min(values) <= 0:
        raise KinisotInputError("the temperature must be positive")
    return values


def _label_text(labels):
    return labels if isinstance(labels, str) else " / ".join(str(x) for x in labels)


def run_job(job, temperature):
    """Every isotopologue of a job at one temperature: a list of (name, result)."""
    results = []
    for entry in job["isotopologues"]:
        options = dict(job["settings"], temperature=temperature)
        if job["kind"] == "channels":
            jobs = []
            for c, channel in enumerate(job["channels"]):
                jobs.append(
                    {
                        "name": channel["name"],
                        "rct": channel["rct"],
                        "ts": channel["ts"],
                        "iso": entry["iso"][c],
                        "reference": entry["reference"][c] if entry["reference"] is not None else None,
                        "amount": channel["amount"],
                    }
                )
            result = channels(
                jobs, shares=job.get("shares"), barriers=job.get("barriers"), amounts=job.get("amounts"),
                energy_unit=job.get("energy_unit", "kcal/mol"), **options,
            )  # fmt: skip
        else:
            result = compute_kie(
                rct=job["rct"], ts=job.get("ts"), prd=job.get("prd"), iso=entry["iso"], reference=entry["reference"],
                **options,
            )  # fmt: skip
        results.append((entry["name"], result))
    return results
