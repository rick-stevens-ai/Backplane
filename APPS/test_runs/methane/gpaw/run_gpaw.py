from ase.build import molecule
from gpaw import GPAW, FermiDirac
atoms = molecule('CH4')
atoms.center(vacuum=4.0)
calc = GPAW(mode='lcao', basis='dzp', xc='PBE', occupations=FermiDirac(0.1), txt='gpaw.log')
atoms.calc = calc
energy = atoms.get_potential_energy()
print(f'TOTAL_ENERGY={energy:.8f} eV')