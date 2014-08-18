import pylab as pl
from nd_binner import proj

def plot_binned_data(ax,edges,data,norm=False,*args,**kwargs):
	#if norm:
	#	data = norm_binned_data(edges,data)
	x = []
	y = [0.]
	for i in range(0,len(edges)):
		x.extend([edges[i],edges[i]])
	for i in range(0,len(data)):
		y.extend([data[i],data[i]])
	y.append(0.)
	ax.plot(x,y,*args,**kwargs)
	length = x[-1]-x[0]
	ax.set_xlim((x[0]-0.1*length,x[-1]+0.1*length))

def integr_binned_data(edges,data):
	assert( len(edges) == len(data)+1)
	integr = 0.
	for i in range(0,len(data)):
		integr+= (edges[i+1]-edges[i])*data[i]
	return integr

def norm_binned_data(edges,data):
	norm = integr_binned_data(edges,data)
	for i in range(0,len(data)):
		data[i] /= norm
	return data
	
def make_3dTo2d_plots(fig,N,edges,cmap=None,binwidthWeighting=True):
	x_edges,y_edges,z_edges = edges[0],edges[1],edges[2]
	ax = fig.add_subplot(221,aspect='equal')
	if binwidthWeighting:
		ax.pcolormesh(x_edges,y_edges,proj(N,z_edges,axis=2).T,cmap=cmap)#,vmin=1e-30)
	else:
		ax.pcolormesh(x_edges,y_edges,np.sum(N,axis=2).T,cmap=cmap)#,vmin=1e-30)
	ax.set_xlabel(r"$P_{12,x}$")
	ax.set_ylabel(r"$P_{12,y}$")
	ax = fig.add_subplot(222,aspect='equal')
	if binwidthWeighting:
		ax.pcolormesh(y_edges,z_edges,proj(N,x_edges,axis=0).T,cmap=cmap)#,vmin=1e-30)
	else:
		ax.pcolormesh(y_edges,z_edges,np.sum(N,axis=0).T,cmap=cmap)#,vmin=1e-30)
	ax.set_xlabel(r"$P_{12,y}$")
	ax.set_ylabel(r"$P_{12,z}$")
	ax = fig.add_subplot(223,aspect='equal')
	if binwidthWeighting:
		CS = ax.pcolormesh(x_edges,z_edges,proj(N,y_edges,axis=1).T,cmap=cmap)#,vmin=1e-30)
	else:
		ax.pcolormesh(x_edges,z_edges,np.sum(N,axis=1).T,cmap=cmap)#,vmin=1e-30)
	pl.colorbar(CS)
	ax.set_xlabel(r"$P_{12,x}$")
	ax.set_ylabel(r"$P_{12,z}$")
	#fig.tight_layout()

def make_3d_shell_plots(edges,shelldata,totaldata):
	fig = pl.figure()
	fig.subplots_adjust(bottom=0.18,left=0.16)
	ax1 = fig.add_subplot(111)
	for s in shelldata:
		plot_binned_data(ax1,edges[0],proj(s,edges=[edges[1],edges[2]],axis=(1,2)),lw=3)
	plot_binned_data(ax1,edges[0],proj(totaldata,edges=[edges[1],edges[2]],axis=(1,2)),lw=3,ls='dashed',color='black')
	
	fig = pl.figure()
	fig.subplots_adjust(bottom=0.18,left=0.16)
	ax2 = fig.add_subplot(111)
	for s in shelldata:
		plot_binned_data(ax2,edges[1],proj(s,edges=[edges[0],edges[2]],axis=(0,2)),lw=3)
	plot_binned_data(ax2,edges[1],proj(totaldata,edges=[edges[0],edges[2]],axis=(0,2)),lw=3,ls='dashed',color='black')
	
	fig = pl.figure()
	fig.subplots_adjust(bottom=0.18,left=0.16)
	ax3 = fig.add_subplot(111)
	for s in shelldata:
		plot_binned_data(ax3,edges[2],proj(s,edges=[edges[0],edges[1]],axis=(0,1)),lw=3)
	plot_binned_data(ax3,edges[2],proj(totaldata,edges=[edges[0],edges[1]],axis=(0,1)),lw=3,ls='dashed',color='black')
	
	return ax1,ax2,ax3
