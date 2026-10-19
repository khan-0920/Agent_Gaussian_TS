"""Chemical conversion helper run in the profile's RDKit-capable prepare environment."""
import sys
from .io import read_json, write_json


def placer_residue(path, atom_map):
    from rdkit import Chem
    molecule = Chem.SDMolSupplier(str(path), removeHs=False)[0]
    if molecule is None:
        raise ValueError('invalid custom ligand SDF')
    molecule = Chem.AddHs(molecule, addCoords=True)
    byindex = {index: name for name, index in atom_map.items()}
    hydrogen, names = 1, []
    for index, atom in enumerate(molecule.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            while 'H%d' % hydrogen in byindex.values():
                hydrogen += 1
            names.append('H%d' % hydrogen)
            hydrogen += 1
        else:
            if index not in byindex:
                raise ValueError('PLACER ligand preparation changed heavy-atom indexing')
            names.append(byindex[index])
    sdf = Chem.MolToMolBlock(molecule) + '\n$$$$\n'
    return {'sdf': sdf, 'atom_id': names, 'leaving': [False]*len(names), 'pdbx_align': [0]*len(names)}


def main():
    from rdkit import rdBase
    request = read_json(sys.argv[1])
    result = {code: placer_residue(item['sdf'], item['atom_map']) for code, item in request.items()}
    write_json(sys.argv[2], {'rdkit_version': rdBase.rdkitVersion, 'residues': result})


if __name__ == '__main__':
    main()
