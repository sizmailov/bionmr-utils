import os
import numpy as np
from tqdm import tqdm
from typing import Tuple, List, Union, Callable, Optional, Iterable, Iterator
from pyxmolpp2 import Trajectory, AtomSelection, Frame, Atom, XYZ, aName
from pyxmolpp2.pipe import Align


def extract_time_per_file_ns(path_to_trajectory: str,
                             subdir="5_run",
                             runin="run00001.in") -> float:
    """
    Extract time step from trajectory

    :param path_to_trajectory: path to directory
    :return time between trajectory frames in nanoseconds

    """
    path_to_firt_run_in = os.path.join(path_to_trajectory, subdir, runin)
    with open(path_to_firt_run_in) as first_run_in:
        for line in first_run_in:
            row = line.strip().split()
            if row[0] == 'ntwx':
                ntwx = int(row[2].strip(","))
            if row[0] == 'dt':
                dt = float(row[2].strip(","))
    time_step_ns = dt / 1000 * ntwx
    return time_step_ns


def extract_rotation_matrices(trajectory: Union[Trajectory, Trajectory.Slice],
                              atom_selector: Callable[[Atom], bool] = (aName == "CA")
                              ) -> np.array:
    """
    Extract rotation matrices

    :param trajectory: trajectory or trajectory slice
    :param atom_selector: selector for atom
    :return: rotation matrices as numpy array of shape (N, 3, 3), where N is length of trajectory

    """
    matrix3d_all_rotations = []
    trajectory = trajectory | Align(by=atom_selector)

    ref = trajectory[0]
    ref_ca = ref.atoms.filter(atom_selector)
    frame_ca = None

    for frame in tqdm(trajectory, desc="extract vectors"):
        if frame_ca is None:
            frame_ca = frame.atoms.filter(atom_selector)
        alignment_transformation = ref_ca.alignment_to(frame_ca)
        matrix_one_rotation = alignment_transformation.matrix3d()
        matrix3d_all_rotations.append(matrix_one_rotation)

    return np.array(matrix3d_all_rotations)


def extract_mass_center(trajectory: Union[Trajectory, Trajectory.Slice],
                        time_per_file: float,
                        atom_selector: Callable[[Atom], bool] = lambda _: True,
                        ) -> Tuple[np.array, np.array]:
    """
    Extract mass center

    :param trajectory: trajectory or trajectory slice
    :param time_per_file: time between trajectory frames in nanoseconds
    :param atom_selector: selector for atom
    :returns: (np.array of time points in nanoseconds, np.array of corresponding XYZ coordinates mass center)

    """
    time = []
    mass_centers = []
    prev_cm = None
    atoms = None
    for frame in trajectory:
        if atoms is None:
            atoms = frame.atoms.filter(atom_selector)

        current_cm = atoms.mean(weighted=True).values
        time.append(frame.index * time_per_file)
        mass_centers.append(current_cm)

        if prev_cm is not None:
            closest = frame.cell.closest_image_to(XYZ(*prev_cm), XYZ(*current_cm))
            current_cm += closest.shift.values
        prev_cm = current_cm

    return np.array(time), np.array(mass_centers)


def extract_vectors(trajectory: Union[Trajectory, Trajectory.Slice, List[Frame]],
                    get_selection: Callable[[Frame], Tuple[AtomSelection, AtomSelection]],
                    alignment_selector: Optional[Callable[[Atom], bool]] = None,
                    ) -> Tuple[List[Tuple[int, str]], Iterable[np.array]]:
    """
    Get vectors from trajectory

    :param trajectory:
    :param get_selection: function to determinate tracked vector
    :param alignment_selector: Optional atom selector for frame alignment
    :return: residue ids and atom names corresponding to atom in the end of chosen vector,
             generator of vectors - list of (x, y, z) coordinates of chosen vectors per frame in trajectory
    """

    if alignment_selector:
        trajectory = trajectory | Align(by=alignment_selector)

    iterator = iter(trajectory)
    frame = next(iterator)

    atoms_selection_1, atoms_selection_2 = get_selection(frame)
    atoms_selection_1_coords, atoms_selection_2_coords = atoms_selection_1.coords, atoms_selection_2.coords
    rids_anames_pairs = [(atom.residue.id.serial, atom.name) for atom in atoms_selection_2]

    def generator_of_vectors() -> Iterator[np.array]:
        # yield vector from first iteration
        yield atoms_selection_1_coords.values - atoms_selection_2_coords.values
        for _ in iterator:
            yield atoms_selection_1_coords.values - atoms_selection_2_coords.values

    return rids_anames_pairs, generator_of_vectors()
