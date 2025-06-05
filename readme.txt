Python3 scripts to calculate phonon SED.

For small systems (< 10 000 atoms), CPU version is fast enough.
For large systems (> 10 000 atoms), if you know how to allocate GPU memory, GPU version is more recommanded.
Becuase the atomic velocity trajectory may be too large to load into GPU at once, you may need to modify some parameters in GPU version script.

To run SED.py, you should provide two files:
1. atomic velocity trajectory file, in the format of numpy binary file, with shape of [n_frames, n_atoms, n_dimensions], in Angstrom/femtosecond
2. .conf file, including following keys:
  - data_file: file name of atomic velocity file
  - out_file:  file name of output SED file
  - lattice:   lattice parameters in angstrom
  - k_points:  k or q points to generate path in the reciprocal space
                 - number of k_points
		 - name, fractional coordinate of k_point
  - k_path:    k or q path in the reciprocal space
                 - number of k_path, spacing of discrete k points on the path
		 - starting k points, ending k points
  - frequency: range of frequency in THz
                 - starting f, ending f, number of slices
  - sample_dt: time interval between frames in your atomic velocity in picosecond
  - mass:      atomic mass in a unit cell in dalton, the sequence should match atomic velocity file
                 - m1, m2, m3, ... 
  - atoms:     information of atoms in the supercell, all indexed start from 0
                 - number of all atoms
		 - atom index in the supercell, atom index in the unit cell, supercell index along the directions of lattice vector a, b, c

For convenience, you can use gen_conf.py to create 'atoms' block in .conf file.
Note that the sequence of atoms in each unit cell should be consistent,
and the sequence of all atoms in 'atoms' block should match the atomic velocity file

Finally, run 'python SED.py your_file.conf'

The resulted SED is in the format of numpy binary file, with shape of [n_frequency, n_k_points]

Here, we provide an example of the oxygen sublattice of ice X with 21*21*21 supercells, including:
  - .conf file:            bcc.conf
  - calculated SED file:   sed.cpu/gpu.npy
The atomic velocity file is too large to upload. You can e-mail to the author or run GPUMD with provided NEP potential.
