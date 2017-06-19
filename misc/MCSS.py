from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import rdkit.Chem.MCS as FMCS

import os

class MCS (object):
	def __init__(self, mol, numatoms, numbonds):
		self.mol = mol
		self.numatoms = int(numatoms)
		self.numbonds = int(numbonds)
		self.total_similar = self.numatoms + self.numbonds
		return

os.system('./PMG C6H6O -o phenol.dat')

# Run MCS search
suppl = Chem.SDMolSupplier('phenol.dat')
ref_mol = Chem.MolFromSmiles('c1ccccc1') # benzene

MCS_results = []

for i, mol in enumerate(suppl):
		if i % 100 == 0: print 'processing entry: ', i
		try:
			res = FMCS.FindMCS([ref_mol, mol]) 
		
			MCS_results.append(
				MCS(mol, res.numAtoms, res.numBonds)
				)
		except Exception, e:
			print str(e)

MCS_results.sort(key = lambda x: x.total_similar, reverse = True)

for index, mol in enumerate(MCS_results):
	if index > 5: break
	AllChem.Compute2DCoords(mol.mol)
	Draw.MolToFile(mol.mol, 'MCS_mol_%s.png' %index)


