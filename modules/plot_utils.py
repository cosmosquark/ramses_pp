#Methods for general plotting

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