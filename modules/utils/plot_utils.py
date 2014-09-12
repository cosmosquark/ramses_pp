#Methods for general plotting
import numpy as np

#Okamoto et al. 2008 gas fraction fitting formula (Hoeft et al. 2006 use same equation)
def okamoto(ax, mhalo, Mc, z, f_bar, alpha, label=None):
	if label == None:
		label = r'Ok+08 $M_{c} = %1.2e M_{\odot} h^{-1}$ $\alpha = %1.2f$'%(Mc, alpha)
	#Okamoto 2008: f_B(M,z) = <f_b> {1 + (2^{alpha/3} - 1) * (M/M_c(z))^{-alpha}}^{-3/alpha}
	mhalo = np.sort(mhalo, axis=None)
	fb_ok = f_bar * ( 1 + ( 2**(float(alpha)/3.) - 1 ) * ( mhalo / Mc )**-alpha )**(-3./float(alpha))
	ax.plot(np.log10(mhalo), fb_ok, label=label, linewidth=2.0)

def legend(ax):
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.1,
				 box.width, box.height * 0.9])
	#ax.set_position([box.x0 - box.width * 0.05, box.y0,
	#			 box.width * 0.9, box.height])

	# Put a legend below current axis
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
		  fancybox=True, shadow=True, ncol=4, prop={'size':7.0})

def adjustFigAspect(fig,aspect=1):
	'''
	Adjust the subplot parameters so that the figure has the correct
	aspect ratio.
	'''
	xsize,ysize = fig.get_size_inches()
	minsize = min(xsize,ysize)
	xlim = .4*minsize/xsize
	ylim = .4*minsize/ysize
	if aspect < 1:
		xlim *= aspect
	else:
		ylim /= aspect
	fig.subplots_adjust(left=.5-xlim,
						right=.5+xlim,
						bottom=.5-ylim,
						top=.5+ylim)