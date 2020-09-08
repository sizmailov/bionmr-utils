from typing import Union, Tuple
import pyxmolpp2


def traj_from_dir(path: str,
                  reference_pdb: Union[str, None] = None,
                  stride: int = 1,
                  first: int = 1,
                  last: int = 1,
                  pattern: str = "run%05d",
                  subdir="5_run",
                  filetype: str = "dat") \
        -> Tuple[pyxmolpp2.Trajectory, pyxmolpp2.Frame]:
    import os

    from pyxmolpp2 import Trajectory, TrjtoolDatFile, AmberNetCDF, PdbFile, GromacsXtcFile, TrajectoryInputFile

    if reference_pdb is None:
        for pdb in [
            os.path.join(path, "5_run", "run00001.pdb"),
            os.path.join(path, "5_run", "step-1.pdb"),
            os.path.join(path, "ref", "ref.pdb")
        ]:
            if os.path.isfile(pdb):
                reference_pdb = pdb
                break
        if reference_pdb is None:
            raise RuntimeError("Reference pdb is not specified")

    if not os.path.isfile(reference_pdb):
        if not os.path.isfile(os.path.join(path, reference_pdb)):
            raise RuntimeError("Reference pdb is not found: `%s`" % reference_pdb)
        reference_pdb = os.path.join(path, reference_pdb)

    ref_frame = PdbFile(reference_pdb).frames()[0]

    coordinate_files = [os.path.join(path, subdir, f"{pattern % i}.{filetype}")
                        for i in range(first, last + 1, stride)
                        ]

    traj = Trajectory(ref_frame)

    portion_type: type

    if filetype == "dat":
        portion_type = TrjtoolDatFile
    elif filetype == "pdb":
        portion_type = PdbFile
    elif filetype == "nc":
        portion_type = AmberNetCDF
    elif filetype == "xtc":
        portion_type = GromacsXtcFile
    else:
        raise RuntimeError("Unknown trajectory coordinate file type `%s`" % filetype)

    for coordinate_file in coordinate_files:
        if not os.access(coordinate_file, os.O_RDONLY):
            raise RuntimeError("Can't access file `%s`" % coordinate_file)
        traj.extend(portion_type(coordinate_file))

    return traj, ref_frame
