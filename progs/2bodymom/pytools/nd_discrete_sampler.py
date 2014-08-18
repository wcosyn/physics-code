##@package nd_discrete_sampler.py
# @author Camille Colle, Ghent University
#
# This script contains a very compact class
# to sample from a given discrete probability
# distribution (given as a array) with arbitrary
# many dimensions (tested for # dimensions 1 to 4)
import numpy as np

## Sample efficiently from a multidimensional array
#  The array should contain positive values, however
#  no explicit checks are done in order to save 
#  computation time.
#  
#  Sampling goes fairly fast, even for python :)
#  1e6 samples from a 40x40x40x40 array takes about 12 seconds
#  on a Intel(R) Core(TM) i3-2120 CPU @ 3.30GHz processor.
class nd_sampler:
	def __init__(self,array):
		self.array = np.array(array,dtype=float)  # make sure the array is a numpy array, note that this copies the array.
		self.ndim  = self.array.ndim  # number of dimensions of the array
		self.P     = [None]*self.ndim # the total size of P is garantueed to be smaller than ndim*size(array), with size array n1*n2*n3*...*nm with ni the array size in the ith dimension obviously can be quite huge for high dimensional arrays
		self.summed= [None]*self.ndim # the total size of self.summed is garantueed to be smaller then (ndim-1)*size(array[:,:,...,:,0]), obviously this can still be quite huge for high dimensional data
		self._construct_marginal_probDistr() #construct the marginal prob distributions needed form sampling
	
	# note the convention is that self.P[index] represents the marginal distribution over dimensions [0,1,...,index], if index == ndim -1 this is not really a marginal distribution anymore, but that is a semantics issue
	def _construct_marginal_probDistr(self):
		self._construct_summed_arrays() #generate summed arrays on beforehand. More memory needed but MUCH less CPU time
		self.P[self.ndim-1] = np.cumsum(self.array,axis=-1) # make a cumsum over last index, still have to normalise this
		temp = np.rollaxis(self.P[self.ndim-1],-1) # bring last index to front, rollaxis does not copy array! It makes a different indexed reference
		for i in range(len(temp)): # normalise the cumsum by iterating over last index (brought to front by rollaxis!)
			temp[i] = np.nan_to_num(temp[i]/temp[-1]) # normalising cumsum is just dividing by last element, if temp[-1] contains zeros "nan"s will be produced, replace 0/0 by 0 trough np.nan_to_num. This will still give an warning though
		# we dont have to roll axis back, self.P[ndim-1] is automatically updated also with above normalisation
		for i in range(self.ndim-2,-1,-1): # iterate in reverse
			self.P[i] = np.cumsum(self.summed[i+1],axis=-1) # get the sum over last index to marginalise the last index (i+1) and make a cumsum to convert into cum prob
			temp = np.rollaxis(self.P[i],-1) # bring last index to front
			for i in range(len(temp)): #normalise the cumsum by iterating over last index (brought to front by rollaxis!)
				temp[i] /= temp[-1] # normalising cumsum is just dividing by last element
	
	# note the convention is that summed[index] represents the sum along all axis in [index,index+1,...,ndim-1]
	def _construct_summed_arrays(self):
		self.summed[self.ndim-1] = np.sum(self.array,axis=-1)
		for i in range(self.ndim-2,-1,-1): #iterate in reverse, compute sums iteratively
			self.summed[i] = np.sum(self.summed[i+1],axis=-1)
	
	def print_marginal_probDistr(self):
		print "original array (shape = ", np.shape(self.array), " ) \n", self.array
		for i in range(self.ndim):
			print " marginal array # %d (shape = " % i , np.shape(self.P[i]), " ), ndim = %d \n " % self.P[i].ndim, self.P[i]
			print ""
	
	def sample(self): # get a sample from our array
		indices = () # create empty tuple
		# now iterate over every dimension/marginal distribution, pulling an index for each dimension, saving and updating the indices
		# note that a (numpy) array accessed [] with a empty tuple just returns the whole array, this is exactly what we want.
		# also remark that as the marginal cum distributions are monotone increasing we can use numpy.searchsorted for log(n) retrieval.
		for i in range(self.ndim): #iterate over all dimensions/marginal distributions
			indices += ( np.searchsorted(self.P[i][indices],np.random.random()), ) # pull a random number, determine what index corresponds to this via numpy.searchsorted and append index to indices tuple
		return indices
	



def test(nsamples):
	# module pylab required
	import pylab as pl
	#! for large arrays you can expect to see large deviations if nsamples is not large enough!
	np.random.seed(18876)
	t = np.random.random((10,10,10))
	n = Nd_sampler(t)
	s = np.zeros(np.shape(t))
	
	for i in range(nsamples): # the actual sampling is this short!
		s[n.sample()] += 1.
	
	while (t.ndim > 2):
		t = np.sum(t,axis=-1) # project to 2D for plotting
	while (s.ndim > 2):
		s = np.sum(s,axis=-1) # project to 2D for plotting
	
	fig = pl.figure()
	fig.suptitle("Left and Right should ideally resemble one another")
	ax = fig.add_subplot(121,aspect='equal')
	ax.pcolormesh(t)
	ax = fig.add_subplot(122,aspect='equal')
	ax.pcolormesh(s)
	
	fig = pl.figure()
	fig.suptitle("Absolute difference, corrected for sample size ")
	ax = fig.add_subplot(111,aspect='equal')
	mfact = np.sum(t)/np.sum(s)
	CS = ax.pcolormesh( (t-s*mfact))
	pl.colorbar(CS)
	
	pl.show()


if __name__=="__main__":
	test(1000000)
