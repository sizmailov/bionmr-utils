import bionmr_utils.data
from bionmr_utils.md import (Frame)


def rename_inplace_charmm_to_amber(frame: Frame):
    """

    :param frame:
    :return: None
    """
    import pandas as pd
    from pkg_resources import resource_stream

    def to_atom_name(s):
        assert len(s) <= 4
        assert s.encode('ascii')
        return s

    def to_residue_name(s):
        assert len(s) <= 4
        assert s.encode('ascii')
        return s

    anames = pd.concat([
        pd.read_csv(resource_stream(bionmr_utils.data.__name__, "rename_tables/anames_%s.csv" % (suffix)))
        for suffix in ["protein", "DNA"]
    ]
    )

    termini_anames = pd.concat([
        pd.read_csv(resource_stream(bionmr_utils.data.__name__, "rename_tables/anames_termini_%s.csv" % (suffix)))
        for suffix in ["protein", "DNA"]
    ]
    )

    rnames = pd.concat([
        pd.read_csv(resource_stream(bionmr_utils.data.__name__, "rename_tables/rnames_%s.csv" % (suffix)))
        for suffix in ["protein", "DNA"]
    ]
    )

    rename_aname_map = {
        (to_residue_name(amber_rName), to_atom_name(row.charmm_aName)): to_atom_name(row.amber_aName)
        for index, row in anames.iterrows()
        for amber_rName in row.amber_rNames.split("|")
    }

    rename_termini_aname_map = {
        (
            row.residue_index, to_residue_name(amber_rName), to_atom_name(row.charmm_aName)
        ): to_atom_name(row.amber_aName)
        for index, row in termini_anames.iterrows()
        for amber_rName in row.amber_rNames.split("|")
    }

    rename_rname_map = {
        (to_residue_name(row.charmm_rName), row.residue_index): to_residue_name(row.amber_rName)
        for index, row in rnames.iterrows()
    }

    # rename atoms in all residues
    def rname_charmm_to_amber(charmm_rName, residue_index=1):
        key = (charmm_rName, residue_index)
        if key in rename_rname_map:
            return rename_rname_map[key]
        return charmm_rName

    # rename atoms in termini residues
    for chain in frame.molecules:
        residues = chain.residues
        for residue_index in [0, -1]:
            for atom in residues[residue_index].atoms:
                termini_key = (residue_index, rname_charmm_to_amber(atom.residue.name, residue_index), atom.name)
                if termini_key in rename_termini_aname_map:
                    print("%s -> %s" % (atom.name, rename_termini_aname_map[termini_key].str))
                    atom.name = rename_termini_aname_map[termini_key]

    # rename atoms in all residues
    for atom in frame.atoms:
        key = (rname_charmm_to_amber(atom.residue.name, residue_index=1), atom.name)
        if key in rename_aname_map:
            atom.name = rename_aname_map[key]

    # rename intermediate residues
    for residue in frame.residues[1:-1]:
        residue.name = rname_charmm_to_amber(residue.name, 1)

    # rename termini residues
    for residue, residue_index in zip(frame.residues[:1] + frame.residues[-1:], [0, -1]):  # type: ignore
        residue.name = rname_charmm_to_amber(residue.name, residue_index)
