import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import glob 

"""
# Snapshots for different k
srcdir = './snaps_alpha=0.5/'
files = glob.glob(srcdir+'k=*.out')
k_list = [float(f.split('/')[-1][2:-4]) for f in files]
data_list = [np.loadtxt(f) for f in files]
snaps = data_list[0].shape[1]-2
xx = data_list[0][:, 0]

target_k = [900, 50, 2, 0.15625]
k_data_list = []

for k, data in zip(k_list, data_list):
	if k in target_k: 
		k_data_list.append((k, data))

colors = pl.cm.jet(np.linspace(0.3,0.7,len(k_data_list)))

for i_snap in range(1, snaps+1): 
	plt.figure()
	for i, (k, data) in enumerate(k_data_list):
		plt.plot(xx, data[:, i_snap+1], color=colors[i], label=r'$\kappa^0=%s$' % k)
	plt.ylim(-0.05, 1.05)
	plt.xlim(xx.min(), xx.max())
	plt.xlabel(r'$X$')
	plt.ylabel(r'$C_O$')
	plt.title(str(i_snap/snaps)+r'$T_{max}$')
	plt.legend(loc='lower right')
	plt.savefig(srcdir+'snap'+str(i_snap)+'.png', dpi=300)
"""


# Snapshots for different alpha 
srcdir = './snaps_k=50/'
files = glob.glob(srcdir+'alpha=*.out')
alpha_list = [float(f.split('/')[-1][6:-4]) for f in files]
data_list = [np.loadtxt(f) for f in files]
snaps = data_list[0].shape[1]-2
xx = data_list[0][:, 0]

colors = pl.cm.jet(np.linspace(0,1,len(alpha_list)))

for i_snap in range(1, snaps+1): 
	plt.figure()
	for i, (alpha, data) in enumerate(zip(alpha_list, data_list)):
		plt.plot(xx, data[:, i_snap+1], color=colors[i], label=r'$\alpha=%s$' % alpha)
	plt.ylim(-0.05, 1.05)
	plt.xlim(xx.min(), xx.max())
	plt.xlabel(r'$X$')
	plt.ylabel(r'$C_O$')
	plt.title(r'$%s T_{max}$' % (i_snap/snaps))
	plt.legend(loc='lower right')
	plt.savefig(srcdir+'snap'+str(i_snap)+'.png', dpi=300)
# plt.show()
