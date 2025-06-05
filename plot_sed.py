import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

data = np.load('sed.cpu.npy')
# flip the frequency axis
data = np.flipud(data)

# truncate large values to balance the color
# or customize a new color bar
# truncation is much easier :D
norm = mpl.colors.LogNorm(vmin=data.min(), vmax=1e-2)

plt.imshow(data, cmap='rainbow', aspect='auto', norm=norm)
plt.xlabel('q path')
plt.ylabel('Frequency (THz)')
plt.xticks([i-0.5 for i in [0, 36, 87, 131, 175]], ['N', 'G', 'H', 'P', 'G'])
plt.yticks([i-0.5 for i in [0, 1000, 2000, 3000, 4000]], [40, 30, 20, 10, 0])
plt.colorbar()
plt.show()

