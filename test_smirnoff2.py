from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import openmm.app as app
import time
mol = Molecule.from_smiles('c1ccccc1-c2ccccc2-c3ccccc3')
mol.generate_conformers(n_conformers=1)
mol.assign_partial_charges('gasteiger')
t0 = time.time()
smirnoff = SMIRNOFFTemplateGenerator(molecules=[mol])
ff = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
ff.registerTemplateGenerator(smirnoff.generator)
modeller = app.Modeller(mol.to_topology().to_openmm(), mol.conformers[0].to_openmm())
sys = ff.createSystem(modeller.topology)
print('System created in', time.time() - t0, 's')
