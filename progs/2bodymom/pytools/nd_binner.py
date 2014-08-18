import numpy as np
from nd_list import nd_list
from nd_discrete_sampler import nd_sampler
import copy


# returns the indices where you would want to put the thing given (single) data element and edges, if point is on an edge, returns the left index
def get_bin_indices(data,edges):
	if not hasattr(data,'__iter__'): # a single array was passed instead of a multidimension array, do some list wrapping
		data = [data]
		edges = [edges]
	
	indices = [None]*len(data)
	for i in range(len(data)):
		indices[i] = np.searchsorted(edges[i],data[i])-1 # determine correct index, here it is crucial that the edges are monotonically increasing! minus one here is to return left edge index
	return indices # can be used for indexing, ( arr[indices])

def bin_events(events,edges,binF): # binF is the function to call of events to get the classification
	if not hasattr(edges,'__iter__'):
		edges = [edges]
	shape = [ len(edge) for edge in edges ]
	H = nd_list(shape,H=[]) # create a nd_list with empty lists to be filled as the items
	for event in events:
		H[ get_bin_indices(binF(event),edges) ].append(event)
	return H

def count_events(H,edges):
	if not hasattr(edges,'__iter__'):
		edges = [edges]
	shape = [ len(edge)-1 for edge in edges ] # determine the shape of the array i am supposed to return
	n = np.zeros(shape,dtype=int)
	for index,value in np.ndenumerate(n):
		n[index] += len(H[index]) #len(get_listElement(H,index))
	return n

def data_binned_events(H,edges,f): # f is the function to calculate the data for the bin given a list of whatever is in H
	n = np.zeros( [len(edge)-1 for edge in edges ],dtype=float)
	for index,_ in np.ndenumerate(n):
		n[index] = f(H[index])
	return n

def get_bin_volumes(edges): # calculate the bins from n-dimensional rectangular array given the edges
	shape = tuple( [ len(edge)-1 for edge in edges ] ) # shape of the volume array to return
	V = np.zeros(shape)
	for index,value in np.ndenumerate(V): # iterate over all elements of the V array
		#print "lower bounds of current bin with indices ", index , " = " , [ edges[i][index[i]] for i in range(len(index)) ]
		#print "upper bounds of current bin with indices ", index , " = " , [ edges[i][index[i]+1] for i in range(len(index)) ]
		#print "lenght of sides of current bin ", [ edges[i][index[i]+1] - edges[i][index[i]] for i in range(len(index)) ]
		#print "volume of bin = ", np.prod([ edges[i][index[i]+1] - edges[i][index[i]] for i in range(len(index)) ])
		V[index] = np.prod([ edges[i][index[i]+1] - edges[i][index[i]] for i in range(len(index)) ])
	return V

# P is the probabilty of choosing the bin
# H contains the lists of events for each bin
# nsamples is total number of events to sample
def sample_array_random(P,H,nsamples):
	sampler = nd_sampler(P)
	samples = [None]*nsamples
	for i in range(nsamples):
		ind = sampler.sample()
		bin_event_List = H[ind]
		samples[i] = bin_event_List[ np.random.randint(0,len(bin_event_List)) ]
	return samples

# loop over list instead of choosing random events in the bin list
# this will reduce the number of oversampled events
def sample_array_ordered(P,H,nsamples):
	N_s = np.zeros(np.shape(P),dtype=int) # this will contain the index for the lists in each bin
	sampler = nd_sampler(P)
	samples = [None]*nsamples
	for i in range(nsamples):
		ind = sampler.sample()
		bin_event_List = H[ind]
		#print N_s[ind]
		#print "ind ", ind
		#print get_listElement(H,ind)
		if (len(bin_event_List) > N_s[ind]):
			samples[i] = bin_event_List[N_s[ind]]
		else:
			print "warning oversampling @ ", ind
			N_s[ind] = 0
			samples[i] = bin_event_List[N_s[ind]]
		N_s[ind] += 1
	return samples

# sample an array but with oversampling forbidden
# If a bin is out of events it is skipped
# Note that this can take a long time if your
# probability distribution is very small in the bins
# where i have to sample to get all the events you requested
def sample_array_oversample_restricted(P,H,sample_edges,nsamples):
	print "======================================"
	print "========== SAMPLER REPORT ============\n"
	print "sample_array_oversample_restricted reporting in..."
	print "  > %d samples requested, MAX samples is %d " % (nsamples,np.sum(count_events(H,sample_edges)))
	assert ( np.sum(count_events(H,sample_edges)) >= nsamples)
	N_s = np.zeros(np.shape(P),dtype=int) # this will contain the index for the lists in each bin
	sampler = nd_sampler(P)
	samples = [None]*nsamples
	i=0
	rej=0
	while i < nsamples:
		ind = sampler.sample()
		bin_event_List = H[ind]
		if (len(bin_event_List) > N_s[ind]):
			samples[i] = bin_event_List[N_s[ind]]
			N_s[ind] += 1
			i+=1
		else:
			#print "warning oversampling @ ", ind, " Rejecting sample! "
			rej +=1
	
	print "  > Due to the restriction on oversampling %d samples were rejected because the bin was out of events" % (rej)
	print "  > You requested %d samples, so rejection/succes ratio is %.2f " % (nsamples,float(rej)/nsamples)
	print "======================================\n"
	return samples

def proj(d,edges,axis=None,ignoreVolume=False): #project down multiple axis
	#print "projection requested on shape ", np.shape(d), " len of edges is ", len(edges), " axis = ", axis
	if (len(np.shape(d)) == 1 ): # we are projecting onto a scalar!
		return proj_single(d,edges,0,ignoreVolume)
	
	else:
		#assert(len(edges)==len(axis))
		if not (hasattr(axis,'__iter__')): # we are projection down a single axis, make compatible with for loop by making them tuples with size 1
			axis = (axis,)
			edges = (edges,)
		t = copy.deepcopy(d)
		for i in range(len(edges)):
			t = proj_single(t,edges[i],axis[i],ignoreVolume)
			axis = tuple( [ axis[j] if j <= i else axis[j]-1 for j in range(len(axis)) ] ) # once we have projected axis[i] out, all axis[j] > axis[i] have to diminish by 1
		return t

def proj_single(d,edges,axis,ignoreVolume): #project down a single axis
	#print "single projection with shape ", np.shape(d), " down axis", axis
	w = [ edges[i+1]-edges[i] for i in range(len(edges)-1) ] if not ignoreVolume else np.ones(len(edges)-1) # widths of bins, if ignoreVolume is on, all just ones
	#print "widths", w
	temp = np.rollaxis(d,axis) # bring the axis "axis" to front to integrate
	assert(len(w) == len(temp))
	p_shape = np.shape(d)[0:axis] + np.shape(d)[axis+1:] # shape of the new matrix, the axis to integrate over is purged from shape
	if len(p_shape) > 0 : # we are projecting onto an array
		p = np.zeros(p_shape)
		for index,value in np.ndenumerate(p): # iterate over all element of the projected matrix
			for i in range(len(w)): # for all elements of the projected matrix, calculate integral
				p[index] += temp[ (i,) + index ]*w[i]
	else: # we are projecting onto a scalar
		p = 0.
		for i in range(len(w)):
			p += temp[i]*w[i]
	return p

def binned_mean(d,edges):
	assert(len(d)==len(edges)-1)
	ret = 0.
	for i in range(len(d)):
		ret += d[i]*(edges[i+1]**2-edges[i]**2)
	return 0.5*ret/proj(d,edges)

def binned_std(d,edges):
	assert(len(d)==len(edges)-1)
	ret = 0.
	for i in range(len(d)):
		ret += d[i]*(edges[i+1]**3-edges[i]**3)
	ret /= proj(d,edges)
	return np.sqrt( 1./3.*ret - binned_mean(d,edges)**2)

if __name__=="__main__":
	print get_bin_indices([0.5,0.2],[[0.0,0.2,0.4,0.6,0.8,1.0],[0.0,0.2,0.4,0.6,0.8,1.0]])
