import os
import pyxmolpp2
import numpy as np
from tqdm import tqdm
from typing import Tuple, List, Union, Callable, Optional, Iterable
from bionmr_utils.md import LatticeVectors, Degrees, Atom, Frame, VectorXYZ, BestShiftFinder, Trajectory


def extract_time_per_file_ns(path_to_trajectory: str) -> float:
    """
    Extract time step from trajectory

    :param path_to_trajectory: path to directory
    :return time between trajectory frames in nanoseconds

    """
    path_to_firt_run_in = os.path.join(path_to_trajectory, "5_run", "run00001.in")
    with open(path_to_firt_run_in) as first_run_in:
        for line in first_run_in:
            row = line.strip().split()
            if row[0] == 'ntwx':
                ntwx = int(row[2].strip(","))
            if row[0] == 'dt':
                dt = float(row[2].strip(","))
    time_step_ns = dt / 1000 * ntwx
    return time_step_ns


def extract_lattice_vectors_rst7(lattice_path: str) -> LatticeVectors:
    """
    Extract lattice vectors from Amber input coordinate (rst7) file

    :param lattice_path: path to input coordinate (rst7) file

    """
    with open(lattice_path, "r") as in_file:
        last_line = in_file.readlines()[-1].split()
        vectors = [float(coordinate) for coordinate in last_line[0:3]]
        angles = [Degrees(float(coordinate)) for coordinate in last_line[3:]]
    return LatticeVectors(vectors[0], vectors[1], vectors[2],
                          angles[0], angles[1], angles[2])


def _atom_name_mass_map(atom: Atom) -> float:
    """
    Mass map for H, C, O, N, P, S atoms

    :param atom: an atom
    :return atomic weight

    """
    atomics_weight = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.973762, 'S': 32.06}
    return atomics_weight[atom.name.str[0]]


def extract_mass_center(traj: Iterable[Frame],
                        time_per_file: float,
                        lattice_vectors: LatticeVectors,
                        volume: np.array = None,
                        atom_selector: Callable[[Atom], bool] = lambda _: True,
                        atom_mass_map: Callable[[Atom], float] = _atom_name_mass_map
                        ) -> Tuple[np.array, VectorXYZ]:
    """
    Extract mass center

    :param traj: trajectory or trajectory slice
    :param time_per_file: time between trajectory frames in nanoseconds
    :param lattice_vectors: translational vectors for periodic boundary conditions
    :param volume: cell volume along the trajectory. Must match *whole* trajectory in length.
                   If None no lattice_vectors rescale is performed
    :param atom_selector: selector for atom
    :param atom_mass_map: function to determinate mass of atom
    :returns: (np.array of time points in nanoseconds, VectorXYZ of corresponding XYZ)

    """
    time = []
    mass_centers = VectorXYZ()
    atoms = None
    prev_cm = None
    for frame in traj:
        if atoms is None:
            atoms = frame.asAtoms.filter(atom_selector)
            masses = [atom_mass_map(atom) for atom in atoms]
            bsf = BestShiftFinder(lattice_vectors)
            reference_volume = lattice_vectors[0].dot(lattice_vectors[1].cross(lattice_vectors[2]))

        if volume is not None:
            factor_scale = (volume[frame.index] / reference_volume) ** (1 / 3)
            bsf.scale_lattice_by(factor_scale)

        current_cm = atoms.mass_center(masses)

        if prev_cm:
            current_cm += bsf.find_best_shift(prev_cm, current_cm)[1]

        time.append(frame.index * time_per_file)
        mass_centers.append(current_cm)

        prev_cm = current_cm

        if volume is not None:
            bsf.scale_lattice_by(1 / factor_scale)
    return np.array(time), mass_centers


def extract_vectors(trajectory: Union[Trajectory, pyxmolpp2.trajectory.TrajectorySlice],
                    get_vectors: Callable[[Frame], List[Tuple[Atom, Atom]]],
                    alignment_selector: Optional[Callable[[Atom], bool]] = None
                    ):
    """
    Get vectors from trajectory

    :param trajectory:
    :param get_vectors: function to determinate tracked vector
    :param alignment_selector: Optional atom selector for frame alignment
    :return dict of (rid, aname): auto-correlation
    """
    if alignment_selector is not None:
        ref = trajectory[0]
        ref_alignment_atoms = ref.asAtoms.filter(alignment_selector)
        alignment_atoms = None
    pair_vectors = None
    vectors = None
    for frame in tqdm(trajectory):
        if alignment_selector is not None:
            if alignment_atoms is None:
                alignment_atoms = frame.asAtoms.filter(alignment_selector)
                frame_atoms = frame.asAtoms
            frame_atoms.transform(alignment_atoms.alignment_to(ref_alignment_atoms))
        if pair_vectors is None:
            pair_vectors = {}
            vectors = {}
            for atom1, atom2 in get_vectors(frame):
                pair_vectors[atom1.rId, atom2.aName] = (atom1, atom2)
                vectors[atom1.rId, atom2.aName] = VectorXYZ()

        for (rid, aname), (atom1, atom2) in pair_vectors.items():
            vectors[rid, aname].append(atom1.r - atom2.r)

    return vectors
