# Written by Ding Hairui.
# E-mail: ding_hairui@foxmail.com

import numpy as np
import sys
import os
# to visualize the progress bar
# you can comment out all the 'tqdm' in this code if you don't need
from tqdm import tqdm

conf = sys.argv[1]


# read .conf file
print('reading .conf file ...\n')
with open(conf)as fp:
	conf = [i for i in fp.readlines()]

idx_data = 0
idx_out = 0
idx_lat = 0
idx_kp = 0
idx_kpath = 0
idx_freq = 0
idx_dt = 0
idx_mass = 0
idx_atoms = 0
for i in range(len(conf)):
	if 'data_file' in conf[i]:
		idx_data = i
	if 'out_file' in conf[i]:
		idx_out = i
	if 'lattice' in conf[i]:
		idx_lat = i
	if 'k_points' in conf[i]:
		idx_kp = i
	if 'k_path' in conf[i]:
		idx_kpath = i
	if 'freq' in conf[i]:
		idx_freq = i
	if 'sample_dt' in conf[i]:
		idx_dt = i
	if 'mass' in conf[i]:
		idx_mass = i
	if 'atoms' in conf[i]:
		idx_atoms = i
		break

data_file = conf[idx_data].strip().split()[-1]
out_file = conf[idx_out].strip().split()[-1]
out_file = out_file
print('data file:', data_file)
print('out file:', out_file, '\n')

lat_a = float(conf[idx_lat+1].strip())
lat = np.array([[float(j) for j in i.strip().split()]for i in conf[idx_lat+2 : idx_lat+5]])
lat *= lat_a
k_lat = np.linalg.inv(lat) * 2 * np.pi
print('lattice =')
print(lat, '\n')
print('k lattice =')
print(k_lat, '\n')

n_ksymb = int(conf[idx_kp+1].strip())
k_symb = [i.split()[0] for i in conf[idx_kp+2 : idx_kp+2+n_ksymb]]
k_coor = [[float(j) for j in i.strip().split()[1:4]] for i in conf[idx_kp+2 : idx_kp+2+n_ksymb]]
kpoints = {k_symb[i]: k_coor[i] for i in range(n_ksymb)}
kpoints = {i: np.matmul(k_lat, np.array(kpoints[i])) for i in kpoints}

n_kpath = int(conf[idx_kpath+1].strip().split()[0])
dk = float(conf[idx_kpath+1].strip().split()[1])
kpath_s = [i.strip().split() for i in conf[idx_kpath+2 : idx_kpath+2+n_kpath]]
print('k path generated with dk =', dk)

nkps = [int(np.ceil(np.linalg.norm(np.array(kpoints[k[0]]) - np.array(kpoints[k[1]])) / dk)) for k in kpath_s]
nkp = sum(nkps)
for i in range(len(kpath_s)):
	print('   %s-%s %s' % (kpath_s[i][0], kpath_s[i][1], nkps[i]))
print('total number of k points:', nkp, '\n')
kpath = [np.linspace(kpoints[kpath_s[i][0]], kpoints[kpath_s[i][1]], nkps[i]).astype(np.float32) for i in range(len(kpath_s))]
kpath_all = np.vstack(kpath)

freq_info = conf[idx_freq+1].strip().split()
freq_start = float(freq_info[0])
freq_end = float(freq_info[1])
n_freq = int(freq_info[2])
freq = np.linspace(freq_start, freq_end, n_freq).astype(np.float32)
omega = freq * 2 * np.pi
print('frequency range from %s to %s THz' % (freq_start, freq_end))
print('total number of frequency =', n_freq, '\n')

dt = float(conf[idx_dt+1].strip().split()[0])

mass = [int(i) for i in conf[idx_mass+1].strip().split()]
mass = np.array(mass).reshape([len(mass), 1])
print('total number of atoms in the unit cell =', mass.shape[0], '\n')

n_atoms = int(conf[idx_atoms+1].strip())
atoms_info = np.array([[int(j) for j in i.split()] for i in conf[idx_atoms+2 : idx_atoms+2+n_atoms]])
atoms_map = atoms_info[:, 1]
atoms_map = [np.where(atoms_map==i)[0] for i in range(mass.shape[0])]
n_sc = atoms_map[0].shape[0]
sc_idx = atoms_info[:, 2:5]
rl = np.array([lat[0]*i[0] + lat[1]*i[1] + lat[2]*i[2] for i in sc_idx]).astype(np.float32)
print('total number of supercells =', n_sc)
print('total number of atoms =', n_atoms, '\n')

print('read .conf file done\n')
print('loading data file ...')


# SED in meV*fs
unit = 1.660539069 / 1.602176634 * 1e5

# load atomic velocity file
v = np.load(data_file)
total_frames = v.shape[0]
t = np.array([i * dt for i in range(total_frames)]).astype(np.float32)
if v.shape[1] != n_atoms:
	raise RuntimeError('number of atoms in data file not equal to the .conf file')
print('load data done')
print('total number of frames =', total_frames, ', data sample dt =', dt, '\n')


# window function to alleviate frequency leaks
def window(seq):
	res = 0.5 * (1 - np.cos(2 * np.pi * np.array(seq) / max(seq)))
	return res.astype(np.float32)
window_f = window(range(total_frames))
window_f = np.reshape(window_f, [total_frames, 1, 1])
window_x = window([i[0] for i in rl])
window_x = np.reshape(window_x, [1, n_atoms, 1])
window_y = window([i[1] for i in rl])
window_y = np.reshape(window_y, [1, n_atoms, 1])
window_z = window([i[2] for i in rl])
window_z = np.reshape(window_z, [1, n_atoms, 1])
v[:, :, 0:1] *= window_x
v[:, :, 1:2] *= window_y
v[:, :, 2:3] *= window_z
v *= window_f


# main calculation of SED

sed = np.zeros([n_freq, nkp])
expt = np.exp(-1j * (omega.reshape([n_freq, 1, 1, 1]) * t.reshape([1, total_frames, 1, 1])))
for k in tqdm(range(nkp)):
	vk = []
	expr = np.exp(1j * (kpath_all[k:k+1] * rl).sum(axis=1)).reshape([1, n_atoms, 1])
	vtmp = v * expr
	va = 1j * np.zeros([total_frames, mass.shape[0], 3])
	for j in range(mass.shape[0]):
		va[:, j, :] = vtmp[:, atoms_map[j], :].sum(axis=1)
	vk.append(va)
	vk = np.vstack(vk).reshape([1, total_frames, mass.shape[0], 3])
	vk = np.abs((vk * expt * dt).sum(axis=1)) ** 2
	vk = (vk * mass).sum(axis=1).sum(axis=1) / 4 / np.pi / total_frames / n_sc
	sed[:, k] = vk

sed *= unit
np.save(out_file, sed)

