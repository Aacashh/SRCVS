from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import openmm.app as app
from rdkit import Chem

_orig_assign = Molecule.assign_partial_charges
def _fast_assign(self, partial_charge_method=None, **kwargs):
    _orig_assign(self, partial_charge_method="gasteiger", **kwargs)
Molecule.assign_partial_charges = _fast_assign

rdmol = next(iter(Chem.SDMolSupplier('out/05_select/final_leads.sdf', removeHs=False)))
off_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
smirnoff = SMIRNOFFTemplateGenerator(molecules=[off_mol])

pdb = app.PDBFile('out/01_target/receptor_prepared.pdb')
forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
forcefield.registerTemplateGenerator(smirnoff.generator)

modeller = app.Modeller(pdb.topology, pdb.positions)
lig_top = off_mol.to_topology().to_openmm()
lig_pos = off_mol.conformers[0].to_openmm()
modeller.add(lig_top, lig_pos)

system = forcefield.createSystem(modeller.topology)
print("System created successfully!")
