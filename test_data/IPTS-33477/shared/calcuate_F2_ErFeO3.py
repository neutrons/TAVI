# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter

def angular_distance(setting1, setting2):
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(setting1, setting2)))

def greedy_sort(settings):
    n = len(settings)
    visited = np.zeros(n, dtype=bool)
    path = []
    current_index = 0
    path.append(current_index)
    visited[current_index] = True

    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            distance_matrix[i, j] = angular_distance(settings[i], settings[j])

    for _ in range(n - 1):
        next_index = np.argmin([distance_matrix[current_index, j] if not visited[j] else np.inf for j in range(n)])
        visited[next_index] = True
        path.append(next_index)
        current_index = next_index
    return [i for i in path]

CreateSingleValuedWorkspace(OutputWorkspace='sample')

LoadCIF(Workspace='sample', InputFile='/HFIR/HB1A/IPTS-33477/shared/ErFeO3/ErFeO3.cif')

generator = ReflectionGenerator(mtd['sample'].sample().getCrystalStructure())

hkls = generator.getHKLsUsingFilter(1.68, 45, ReflectionConditionFilter.StructureFactor)

ds = np.array(generator.getDValues(hkls))
F2s = np.array(generator.getFsSquared(hkls))

hkls = np.array(hkls)

UB = np.loadtxt('/HFIR/HB1A/IPTS-33477/shared/ErFeO3/UBmatrix.dat')
#UB = np.loadtxt('/HFIR/HB1A/IPTS-33477/shared/hexaferrite/UBmatrix.dat.hexaferrite')

*lattice, lamda = np.loadtxt('/HFIR/HB1A/IPTS-33477/shared/ErFeO3/lattice.dat')

SetUB(Workspace='sample', UB=UB)

B = mtd['sample'].sample().getOrientedLattice().getB()

Qc = np.einsum('ij,kj->ik', B, hkls)

q = np.sqrt(np.sum(Qc**2, axis=0))

theta = np.arcsin(lamda*q/2)
theta2 = theta*2*180/np.pi

omega = theta2/2

hphi = np.einsum('ij,kj->ik', UB, hkls)

newphi = np.zeros_like(ds)
newchi = np.zeros_like(ds)

for i in range(len(ds)):
    newphi[i] = np.arctan(hphi[1,i]/hphi[0,i])*180/np.pi
    if newphi[i] > 0 and hphi[0,i] < 0:
        newphi[i] -= 180
    elif newphi[i] <= 0 and hphi[0,i] < 0:
        newphi[i] += 180

    newchi[i] = np.arctan(hphi[2,i]/np.sqrt(hphi[0,i]**2+hphi[1,i]**2))*180/np.pi

mask = (theta2 >= 10) & (theta2 < 90) & (newchi > -70) & (newchi < 40) & (np.abs(theta2-61) > 3) & (np.abs(theta2-72) > 3)

hkls = hkls[mask]
ds = ds[mask]
F2s = F2s[mask]
theta2 = theta2[mask]
omega = omega[mask]
newphi = newphi[mask]
newchi = newchi[mask]

weights = [45/50, 42/40, 100/76]

sort = greedy_sort(np.column_stack([weights[0]*theta2, weights[1]*newchi, weights[2]*newphi]))

with open('/HFIR/HB1A/IPTS-33477/shared/ErFeO3/scanlist_Xtal.dat', 'w') as f:
    for i in sort:
        f.write('{:4.0f}{:4.0f}{:4.0f}{:12.2f}{:12.4f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}\n'.format(*hkls[i], F2s[i], ds[i], theta2[i], omega[i], newchi[i], newphi[i]))
