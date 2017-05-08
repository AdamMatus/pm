import numpy as np
import matplotlib as mpl
mpl.use('pgf')

def figsize(scale):
    fig_width_pt = 455.24408                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "text.fontsize": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{polski}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

# I make my own newfig and savefig functions
def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))


import scipy.io as sio
# vowel 'e'
fig, ax = newfig(0.5)

mat = sio.loadmat('./e.mat')
e = mat['e']

ax.plot(e)
ax.plot([106,183],[0.5,0.5])
ax.text(110,0.55,'77 probek')
plt.ylim([-0.4,0.7])

savefig('e_vowel')

# cep d
fig, ax = newfig(0.5)

mat = sio.loadmat('./d.mat')
d = mat['d']

ax.plot(d)

plt.annotate('pierwsza \n"rahmoniczna"', xy=(77, 2), xytext=(115, 4),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

plt.annotate('druga \n"rahmoniczna"', xy=(154, -2), xytext=(180, -5.5),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
plt.ylim([-7,7])
savefig('e_cepstrum')
