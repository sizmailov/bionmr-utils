from typing import Union, Tuple, Callable
import pyxmolpp2


def traj_from_dir(path: str,
                  reference_pdb: Union[str, None] = None,
                  stride: int = 1,
                  first: int = 1,
                  last: int = 1,
                  pattern: str = "run%05d",
                  subdir="5_run",
                  filetype: str = "dat") \
        -> Tuple[pyxmolpp2.trajectory.Trajectory, pyxmolpp2.polymer.Frame]:
    import os
    from tqdm import tqdm

    from pyxmolpp2.trajectory import Trajectory
    from pyxmolpp2.trjtool import DatFile
    from pyxmolpp2.amber import NetCDFTrajectoryFile
    from pyxmolpp2.pdb import PdbFile, \
        AlteredPdbRecords, \
        StandardPdbRecords, \
        RecordName, \
        FieldName

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

    altered_records = AlteredPdbRecords(StandardPdbRecords.instance())
    altered_records.alter_record(RecordName("ATOM"), FieldName("serial"), [7, 12])

    ref_frame = PdbFile(reference_pdb, altered_records).get_frame()

    coordinate_files = [os.path.join(path, subdir, "%s.%s" % ((pattern % i), filetype)) for i in
                        range(first, last + 1, stride)]

    traj = Trajectory(ref_frame, check_portions_to_match_reference=True)

    if filetype == "dat":
        portion_type = DatFile  # type: Callable[[str], pyxmolpp2.trajectory.TrajectoryPortion]
    elif filetype == "pdb":
        portion_type = lambda filename: PdbFile(filename, altered_records)  # noqa: E731
    elif filetype == "nc":
        portion_type = NetCDFTrajectoryFile
    else:
        raise RuntimeError("Unknown trajectory coordinate file type `%s`" % filetype)

    for coordinate_file in tqdm(coordinate_files, leave=False, desc="checking input files"):
        if not os.access(coordinate_file, os.O_RDONLY):
            raise RuntimeError("Can't access file `%s`" % coordinate_file)
        traj.push_trajectory_portion(portion_type(coordinate_file))

    return traj, ref_frame
