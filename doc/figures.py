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

# melfb

fig, ax = newfig(1)

mat = sio.loadmat('./melfb.mat');
t = mat['t']
y = mat['y']

ax.plot(t,y)
ax.set_xlabel('Czestotliwosc [Hz]');
savefig('melfb')

import numpy as np
#mel

fig, ax = newfig(1)

mat = sio.loadmat('./mel.mat');
ty = mat['t']
yy = mat['y']

ax.plot(ty[0],yy[0],'b')
ax.set_xlabel('Czestotliwosc [Hz]');
ax.set_ylabel('Skala mela [Mel]');
savefig('mel')
#spectrum

fig, ax = newfig(1)

mat = sio.loadmat('./spectrum.mat');
spec = mat['spectrum'];
t = np.linspace(0,6250, num=128)

plt.annotate('formant', xy=(600, 2), xytext=(1200, 2.2),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

plt.annotate('formant', xy=(2600, 1.2), xytext=(3000, 2),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

plt.annotate('formant', xy=(4200, -0.2), xytext=(4800, 0.5),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

plt.annotate('formant', xy=(4850, -1.3), xytext=(5300, -0.5),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

ax.plot(t,spec);
ax.axes.get_yaxis().set_ticks([])
ax.set_xlabel('Czestotliwosc [Hz]');
savefig('spectrum')

# DTW
#
x = np.arange(0, 6.28, 0.1)
y = np.cos(x)

fig, ax = newfig(1)

ax0 = plt.subplot2grid((3,3), (0, 0), rowspan=2)
ax0.plot(y, x)
ax0.set_ylim([0, 6.28])
ax0.axes.get_xaxis().set_ticks([])
ax0.set_title('Sygnal wzorcowy')

ax1 = plt.subplot2grid((3,3), (0, 1), rowspan=2, colspan=2)
ax1.plot([0,0.33*3.14],[0,3.14],'g', label='fun 1')
ax1.plot([0.33*3.14,3.14],[3.14,6.28],'g')
ax1.plot([0,(3.0/2)*3.14],[0,6.28],'r', label='fun 2')
ax1.plot([0,6.28],[0,6.28],'k', label='fun 3')
ax1.plot([3.14, 6.28],[0,6.28],'b', label='fun 4')
ax1.set_ylim([0, 6.28])
ax1.set_xlim([0, 6.28])
ax1.set_title('Trajektorie najmniejszej zagregowanej odleglosci')
ax1.set_xlabel('t[s] - sygnal prownywany')
ax1.set_ylabel('t[s] - sygnal wzorcowy')

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc=2)

ax2 = plt.subplot2grid((3,3), (2, 1), colspan=2)
x = np.arange(0, 1.5*3.14, 0.1)
y = -np.cos(1.3*x)
ax2.plot(x,y,'r');
x = np.arange(3.14, 6.28, 0.1)
y = -np.cos(2*x)
ax2.plot(x,y,'b');
x = np.arange(0, 6.28, 0.1)
y = -np.cos(x)
ax2.plot(x,y,'k');
x = np.arange(0, (1.0/6)*6.28, 0.1)
y = -np.cos(3*x)
ax2.plot(x,y,'g');
x = np.arange((1.0/6)*6.28, 4*3.14/3, 0.1)
y = -np.cos((x)+(2.0/3)*3.14)
ax2.plot(x,y,'g');
ax2.set_xlim([0, 6.28])
ax2.set_title('Sygnaly porownywane')

plt.tight_layout()
savefig('dtw')

#lbz

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111)



ax0 = plt.subplot2grid((3,2), (0,0))
mat = sio.loadmat('./D1.mat');
D = mat['D']
mat = sio.loadmat('./C1.mat');
C = mat['C']
ax0.plot(D[0],D[1],',')
ax0.plot(C[0],C[1],'.r')

ax0.set_title('Etap podwajania')

ax0 = plt.subplot2grid((3,2), (0,1))
mat = sio.loadmat('./D2.mat');
D = mat['D']
mat = sio.loadmat('./C2.mat');
C = mat['C']
ax0.plot(D[0],D[1],',')
ax0.plot(C[0],C[1],'.r')

ax0.set_title('Etap podwajania')

ax0 = plt.subplot2grid((3,2), (1,0))
mat = sio.loadmat('./D3.mat');
D = mat['D']
mat = sio.loadmat('./C3.mat');
C = mat['C']
ax0.plot(D[0],D[1],',')
ax0.plot(C[0],C[1],'.r')


ax0 = plt.subplot2grid((3,2), (1,1))
mat = sio.loadmat('./D4.mat');
D = mat['D']
mat = sio.loadmat('./C4.mat');
C = mat['C']
ax0.plot(D[0],D[1],',')
ax0.plot(C[0],C[1],'.r')


ax0 = plt.subplot2grid((3,2), (2,0))
mat = sio.loadmat('./D5.mat');
D = mat['D']
mat = sio.loadmat('./C5.mat');
C = mat['C']
ax0.plot(D[0],D[1],',')
ax0.plot(C[0],C[1],'.r')

ax0 = plt.subplot2grid((3,2), (2,1))
mat = sio.loadmat('./D6.mat');
D = mat['D']
mat = sio.loadmat('./C6.mat');
C = mat['C']
ax0.plot(D[0],D[1],',')
ax0.plot(C[0],C[1],'.r')

savefig('lbz')
