import os
os.chdir('..')
load('fpMKIII.sage')
Z=BinaryQF(81,44,6)
with open('data/Sp4Z.20_Ups.json', 'r') as f:
	L = json.load(f)
E = SPMF(L)
F=twistForms(E,3)
print(F)