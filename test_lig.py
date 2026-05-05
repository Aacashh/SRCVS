from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import openmm.app as app
from rdkit import Chem

_orig_assign = Molecule.assign_partial_charges
def _fast_assign(self, partial_charge_method=None, **kwargs):
    print(f"Bypassing {partial_charge_method} and using gasteiger")
    _orig_assign(self, partial_charge_method="gasteiger", **kwargs)
Molecule.assign_partial_charges = _fast_assign

rdmol = next(iter(Chem.SDMolSupplier('out/05_select/final_leads.sdf', removeHs=False)))
off_mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
print("Molecule loaded")
smirnoff = SMIRNOFFTemplateGenerator(molecules=[off_mol])
ff = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
ff.registerTemplateGenerator(smirnoff.generator)
modeller = app.Modeller(off_mol.to_topology().to_openmm(), off_mol.conformers[0].to_openmm())
sys = ff.createSystem(modeller.topology)
print("System created")
