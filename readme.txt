Python scripts to calculate phonon SED.

For small systems (< 10 000 atoms), CPU version is fast enough.
For large systems (> 10 000 atoms), GPU version is more recommanded.

To run SED.py, you should provide two files:
1. atomic velocity, in the format of numpy binary file, with a shape of [n_frames, n_atoms, n_dimensions], in Angstrom/femtosecond
2. .conf file, including following keys:
  - data_file: path of atomic velocity file
  - out_file:  path of output SED file
  - lattice:   lattice parameters in angstrom
  - k_points:  k or q points to generate path in the reciprocal space
  - k_path:    k or q path in the reciprocal space
      - number of k_path, spacing of discrete k points on the path
      - starting k points, ending k points
  - frequency: range of frequency in THz
      - starting f, ending f, number of interpolation
  - sample_dt: time interval between frames in your atomic velocity in picosecond
  - mass:      atomic mass in a unit cell, the sequence should match atomic velocity file
  - atoms:     information of atoms in the supercell, all indexed start from 0
      - number of all atoms
      - atom index in the supercell, atom index in the unit cell, supercell index along the directions of lattice vector a, b, c

You can use gen_conf.py to create 'atoms' block in .conf file.
Note that the sequence of atoms in each unit cell should be consistent,
and the sequence of all atoms in 'atoms' block should match the atomic velocity file

Finally, run 'python sed.py your_file.conf'

The resulted SED is saved to numpy binary file, with a shape of [n_frequency, n_k_points]

Due to the atomic velocity file is too large to upload,
we just provide an example of a very very small ice X.

To get a high-quality SED of ice,
you can use the parameters we provided in our paper,
and run MD with nep.txt here.
