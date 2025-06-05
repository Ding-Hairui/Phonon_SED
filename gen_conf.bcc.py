# !!!!! you should modify this script according to how you create your supercell
# In my calculation, I create the supercell based on a two-oxygen-atom cubic unit cell and run MD
# But finally I use the BCC Bravais lattice to calculate phonon SED
# So in each loop, two oxygen atoms are appended

# supercell number
la = lb = lc = 21
filename = 'bcc.conf'

with open(filename)as fp:
	lines = [i for i in fp.readlines()]

idx = 0
for i in range(len(lines)):
	if 'atoms' in lines[i]:
		idx = i+1
		break

new_lines = lines[:idx]

atom_idx = []
atom_unit_cell = []
la_idx = []
lb_idx = []
lc_idx = []
count = 0
for i1 in range(la):
	for i2 in range(lb):
		for i3 in range(lc):
			atom_idx.append(count)
			atom_unit_cell.append(0)
			la_idx.append(i2+i3)
			lb_idx.append(i1+i3)
			lc_idx.append(i1+i2)
			count += 1

			atom_idx.append(count)
			atom_unit_cell.append(0)
			la_idx.append(i2+i3+1)
			lb_idx.append(i1+i3+1)
			lc_idx.append(i1+i2+1)
			count += 1

new_lines.append('%s\n' % (count))
for i in range(count):
	new_lines.append('%s %s %s %s %s\n' % (atom_idx[i], atom_unit_cell[i], la_idx[i], lb_idx[i], lc_idx[i]))

with open(filename, 'w')as fp:
	fp.writelines(new_lines)

