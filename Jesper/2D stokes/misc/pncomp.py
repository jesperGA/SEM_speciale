import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

# Reading the data from CSV files
Pn = pd.read_csv('misc/PnPn_error.csv', header=None).values
Pn2 = pd.read_csv('misc/PnPn-2_error.csv', header=None).values

# Setting up the plot's appearance and using LaTeX for text rendering
plt.figure(figsize=[10, 5])
plt.rcParams.update({
    'font.size': 16,
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts}'
})

# Getting the default color cycle for plotting
# DC = plt.rcParams['axes.prop_cycle'].by_key()['color']

colors_df = pd.read_csv('misc/colors.csv', header=None)
DC = colors_df.values / 1  # Normalize RGB values to 0-1

# Plotting the data
plt.semilogy(Pn[:,0], Pn[:,1], '-o', label=r'$u \in \mathbb{P}_N$', linewidth=2, markersize=10, color=DC[1])
plt.semilogy(Pn[:,0], Pn[:,2], '-o', label=r'$p \in \mathbb{P}_N$', linewidth=2, markersize=10, color=DC[3])
plt.semilogy(Pn2[:,0], Pn2[:,1], '-o', label=r'$u \in \mathbb{P}_{N-2}$', linewidth=2, markersize=10, color=DC[0])
plt.semilogy(Pn2[:,0], Pn2[:,2], '-o', label=r'$p \in \mathbb{P}_{N-2}$', linewidth=2, markersize=10, color=DC[2])

# Adding labels and legend with LaTeX formatting
plt.xlabel('$n_{dof}$', fontsize=18, labelpad=20)
plt.ylabel('$\\| Error \\|_{\infty}$', fontsize=18, labelpad=20)
plt.legend(loc='upper center', ncol=2, fontsize=18,bbox_to_anchor=(0.5, 1.2))

# Set x-ticks to only display integers
ax = plt.gca()  # Get current axes
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_yscale('log')

# Turn on minor ticks and minor grid
ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth='0.5', color='gray')
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='gray', which='minor')

# Additional plot customizations
plt.grid(True)
plt.show()