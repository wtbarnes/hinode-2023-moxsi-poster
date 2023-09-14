"""
Example script for creating a figure in PGF format
"""
import pathlib

import matplotlib.pyplot as plt
import numpy as np


figure_dir = pathlib.Path(__file__).parent.parent / 'src' / 'figures'
angle = np.linspace(-2*np.pi, 2*np.pi, 1000)

y_1 = np.cos(angle)
y_2 = np.sin(angle)

# This is the width of a column in the poster
# as defined by the \colwidth variable. Adjust this
# depending on the column width you use and then
# scale all figure sizes to this variable
column_width = 13.5028346256
golden_ratio = (1 + np.sqrt(5))/2

fig = plt.figure(figsize=(0.65*column_width, 0.65*column_width/golden_ratio), layout='constrained')
ax = fig.add_subplot()
ax.plot(angle, y_1, label=r'$\cos{(x)}$')
ax.plot(angle, y_2, label=r'$\sin{(x)}$')
ax.set_xlabel(r'$\phi$ [rad]')
ax.set_ylabel(r'$y$')
ax.legend(loc=3)
ax.set_title('A very interesting title')

fig.savefig(figure_dir / 'example_plot.pgf')
