"""ASE backend: Hessians from any ASE calculator (machine-learned potentials included).

Two entry points:

* ``parse_ase_json(path)`` reads a Hessian saved as the dictionary of
  :class:`ase.vibrations.VibrationsData` (``todict()`` encoded with
  ``ase.io.jsonio``), which is also what ``save_hessian_json`` writes.
* ``hessian_from_calculator(atoms, calculator)`` computes the Hessian with
  central finite differences (``ase.vibrations.Vibrations``) or, when the
  calculator provides ``get_hessian``, analytically.

``hessian_for_geometry`` combines both for the command line: it computes the
Hessian of a geometry file with ``--calc`` and caches it next to the file as
``<name>.hessian.json`` so that a scan over substitutions and temperatures
computes it once.

ASE is imported lazily; install it with ``pip install kinisot[ase]``. Units
are converted from ASE's eV and Angstrom to Hartree and Bohr; the masses
ASE carries (standard atomic weights) only identify elements, both
isotopologues come from Kinisot's isotope table. Because finite-difference
Hessians leave sizeable translational/rotational residuals, the API projects
the external modes out by default for inputs from this backend.
"""

import importlib
import os

import numpy as np

from ..exceptions import KinisotInputError, KinisotParseError
from ..hessian import HessianInput, linear_from_geometry
from ..isotopes import element_symbol, light_mass

__all__ = [
    "parse_ase_json",
    "is_ase_json",
    "hessian_from_calculator",
    "hessian_for_geometry",
    "save_hessian_json",
    "build_calculator",
    "CALCULATORS",
]

# --calc names -> (module, callable, default kwargs, pip hint). Any other value is taken as "module:callable".
CALCULATORS = {
    "emt": ("ase.calculators.emt", "EMT", {}, "ase (built in; a test potential, not for chemistry)"),
    "mace_mp": ("mace.calculators", "mace_mp", {"default_dtype": "float64"}, "mace-torch"),
    "mace_off": ("mace.calculators", "mace_off", {"default_dtype": "float64"}, "mace-torch"),
    "mace_omol": ("mace.calculators", "mace_omol", {"default_dtype": "float64"}, "mace-torch"),
    "orb": (
        "orb_models.forcefield.pretrained",
        "orb_v3_conservative_inf_omat",
        {},
        "orb-models (wrap with ORBCalculator)",
    ),
    "sevennet": ("sevenn.calculator", "SevenNetCalculator", {}, "sevenn"),
    "aimnet2": ("aimnet2calc", "AIMNet2ASE", {}, "aimnet2calc"),
}


def _require_ase():
    try:
        import ase  # noqa: F401
        import ase.units  # noqa: F401
    except ImportError:
        raise KinisotInputError(
            "this input needs the Atomic Simulation Environment: pip install ase (or pip install kinisot[ase])"
        ) from None


def build_calculator(spec):
    """Instantiate an ASE calculator from a --calc specification.

    ``spec`` is a name from CALCULATORS, optionally followed by ``:model``
    (e.g. ``mace_mp:medium``), or ``module.path:callable`` for anything else.
    The callable is called with no positional arguments (plus ``model=`` when given).
    """
    _require_ase()
    name, _, model = spec.partition(":")
    kwargs = {}
    if name in CALCULATORS:
        module_name, attribute, kwargs, hint = CALCULATORS[name]
        kwargs = dict(kwargs)
    elif "." in name and model:
        module_name, attribute, hint = name, model, name
        model = ""
    else:
        raise KinisotInputError(
            "unknown calculator %r: use one of %s, or module.path:callable" % (spec, ", ".join(CALCULATORS))
        )
    try:
        module = importlib.import_module(module_name)
    except ImportError as err:
        raise KinisotInputError(
            "cannot import %s for --calc %s (%s); install %s" % (module_name, spec, err, hint)
        ) from None
    try:
        factory = getattr(module, attribute)
    except AttributeError:
        raise KinisotInputError("%s has no %s (--calc %s)" % (module_name, attribute, spec)) from None
    if model:
        kwargs["model"] = model
    try:
        return factory(**kwargs)
    except TypeError:
        return factory(model) if model else factory()


def _to_hessian_input(atoms, hessian_ev_a2, source, program_note="ase", energy_ev=None, frequencies=None):
    from ase import units

    hessian = np.asarray(hessian_ev_a2, dtype=float).reshape(3 * len(atoms), 3 * len(atoms))
    hessian = 0.5 * (hessian + hessian.T) * units.Bohr**2 / units.Hartree
    atomic_numbers = tuple(int(z) for z in atoms.get_atomic_numbers())
    positions = np.asarray(atoms.get_positions(), dtype=float) / units.Bohr
    masses = tuple(light_mass(element_symbol(z)) for z in atomic_numbers)
    return HessianInput(
        hessian=hessian,
        masses=masses,
        atomic_numbers=atomic_numbers,
        source=source,
        program=program_note,
        level_of_theory=atoms.info.get("level_of_theory"),
        linear=linear_from_geometry(positions, masses),
        positions=positions,
        program_frequencies=tuple(frequencies) if frequencies is not None else None,
        energy=energy_ev / units.Hartree if energy_ev is not None else None,
    )


def _frequencies_with_masses(vibrations_data, masses):
    """ASE's frequencies (cm-1, negative for imaginary) with Kinisot's light masses."""
    freqs = vibrations_data.with_new_masses(list(masses)).get_frequencies()
    real = np.where(np.iscomplex(freqs), -np.abs(freqs), np.real(freqs)).astype(float)
    order = np.argsort(np.abs(real))
    n_external = 5 if linear_from_geometry(vibrations_data.get_atoms().get_positions(), masses) else 6
    return sorted(real[order[n_external:]])


def is_ase_json(path):
    """Whether ``path`` is a JSON file holding a VibrationsData dictionary."""
    if os.path.splitext(path)[1].lower() != ".json":
        return False
    try:
        with open(path, encoding="utf-8") as handle:
            head = handle.read(4096)
    except OSError:
        return False
    return '"hessian"' in head and '"atoms"' in head


def parse_ase_json(path):
    """Read a Hessian saved as VibrationsData JSON (``save_hessian_json`` or ``VibrationsData.todict``)."""
    _require_ase()
    from ase.io.jsonio import decode
    from ase.vibrations import VibrationsData

    try:
        with open(path, encoding="utf-8") as handle:
            data = decode(handle.read())
    except (OSError, ValueError, KeyError) as err:
        raise KinisotParseError("%s: cannot read as VibrationsData JSON (%s)" % (path, err)) from None
    try:
        vibrations = VibrationsData.fromdict(data)
    except (KeyError, ValueError, TypeError, AssertionError) as err:
        raise KinisotParseError("%s: not a VibrationsData dictionary (%s)" % (path, err)) from None
    atoms = vibrations.get_atoms()
    if len(vibrations.get_indices()) != len(atoms):
        raise KinisotParseError(
            "%s: the Hessian covers %d of %d atoms; Kinisot needs all atoms free"
            % (path, len(vibrations.get_indices()), len(atoms))
        )
    energy = atoms.info.get("energy")
    masses = tuple(light_mass(element_symbol(int(z))) for z in atoms.get_atomic_numbers())
    return _to_hessian_input(
        atoms,
        vibrations.get_hessian_2d(),
        source=path,
        program_note=atoms.info.get("program", "ase"),
        energy_ev=energy,
        frequencies=_frequencies_with_masses(vibrations, masses),
    )


def hessian_from_calculator(atoms, calculator, delta=0.01, nfree=2, analytic=True, source=None, name=None):
    """Compute the Hessian of ``atoms`` with an ASE ``calculator`` and return a HessianInput.

    Uses the calculator's analytic ``get_hessian`` when it has one and
    ``analytic`` is true, otherwise central finite differences of the forces
    (``ase.vibrations.Vibrations``: ``delta`` in Angstrom, ``nfree`` 2 or 4)
    in a temporary directory (or ``name`` as cache directory). The energy is
    stored for the Skodje-Truhlar barrier. The geometry must be a stationary
    point of the same calculator.
    """
    _require_ase()
    import tempfile

    from ase.vibrations import Vibrations, VibrationsData

    atoms = atoms.copy()
    atoms.calc = calculator
    energy = float(atoms.get_potential_energy())
    hessian = None
    if analytic and hasattr(calculator, "get_hessian"):
        try:
            hessian = np.asarray(calculator.get_hessian(atoms), dtype=float)
        except Exception:  # noqa: BLE001 -- fall back to finite differences on any calculator error
            hessian = None
    if hessian is None:
        if name is None:
            with tempfile.TemporaryDirectory() as directory:
                vib = Vibrations(atoms, name=os.path.join(directory, "vib"), delta=delta, nfree=nfree)
                vib.run()
                hessian = vib.get_vibrations().get_hessian_2d()
        else:
            vib = Vibrations(atoms, name=name, delta=delta, nfree=nfree)
            vib.run()
            hessian = vib.get_vibrations().get_hessian_2d()
    vibrations = VibrationsData.from_2d(atoms, np.asarray(hessian).reshape(3 * len(atoms), 3 * len(atoms)))
    masses = tuple(light_mass(element_symbol(int(z))) for z in atoms.get_atomic_numbers())
    atoms.info["energy"] = energy
    atoms.info.setdefault("program", "ase:%s" % type(calculator).__name__)
    return _to_hessian_input(
        atoms,
        vibrations.get_hessian_2d(),
        source=source or "%s (%s)" % (atoms.get_chemical_formula(), type(calculator).__name__),
        program_note=atoms.info["program"],
        energy_ev=energy,
        frequencies=_frequencies_with_masses(vibrations, masses),
    )


def save_hessian_json(data, path, atoms=None, calc_spec=None):
    """Write a HessianInput as VibrationsData JSON so it can be reused as an input file.

    ``atoms`` (an ASE Atoms) supplies the geometry when given; otherwise it is
    rebuilt from the HessianInput's atomic numbers and positions.
    """
    _require_ase()
    from ase import Atoms, units
    from ase.io.jsonio import encode
    from ase.vibrations import VibrationsData

    if atoms is None:
        if data.positions is None:
            raise KinisotInputError("cannot save %s as ASE JSON: it has no geometry" % data.source)
        atoms = Atoms(numbers=list(data.atomic_numbers), positions=np.asarray(data.positions) * units.Bohr)
    atoms = atoms.copy()
    atoms.calc = None
    if data.energy is not None:
        atoms.info["energy"] = data.energy * units.Hartree
    if data.level_of_theory:
        atoms.info["level_of_theory"] = data.level_of_theory
    atoms.info["program"] = data.program or "ase"
    if calc_spec:
        atoms.info["kinisot_calc"] = calc_spec
    hessian_ev = np.asarray(data.hessian) * units.Hartree / units.Bohr**2
    vibrations = VibrationsData.from_2d(atoms, hessian_ev)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(encode(vibrations.todict()))
        handle.write("\n")
    return path


def hessian_for_geometry(path, calc_spec, delta=0.01, nfree=2, cache=True, recompute=False):
    """Hessian of the geometry file ``path`` with the calculator ``calc_spec``, cached as ``<stub>.hessian.json``."""
    _require_ase()
    from ase.io import read
    from ase.io.jsonio import decode

    stub = os.path.splitext(path)[0]
    cached = stub + ".hessian.json"
    if cache and not recompute and os.path.exists(cached) and os.path.getmtime(cached) >= os.path.getmtime(path):
        with open(cached, encoding="utf-8") as handle:
            info = decode(handle.read())["atoms"].info
        if info.get("kinisot_calc") == calc_spec:
            return parse_ase_json(cached)
    try:
        atoms = read(path)
    except Exception as err:  # noqa: BLE001 -- ase.io raises a variety of types
        raise KinisotParseError("%s: cannot read the geometry with ASE (%s)" % (path, err)) from None
    if isinstance(atoms, list):
        atoms = atoms[-1]
    data = hessian_from_calculator(atoms, build_calculator(calc_spec), delta=delta, nfree=nfree, source=path)
    if cache:
        save_hessian_json(data, cached, atoms=atoms, calc_spec=calc_spec)
        data = parse_ase_json(cached)
    return data
