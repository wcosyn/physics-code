## @package nd_list.py
#  A multidimensional list that supports indexing by arrays and tuples
#  
#  Note that there is no method to get the shape from the array as one can
#  modify this at will. There is no garantee that self.dims will be preserved
#  if __setitem__ is called
#
#
#  @author Camille Colle, Ghent University

import copy

class nd_list(object):
	def __init__(self,dims,H=None): # init multidimensional list
		self.dims = dims
		self.H    = self._constr_nd_list(level=len(dims)-1,H=H)
	
	def _constr_nd_list(self,level,H=None):
		#print "constr called with level ", level, " and H ", H, " and dims ", self.dims
		if (level==-1):
			return H
		else:
			Hnew = []
			for i in range(self.dims[level]):
				Hnew.append(copy.deepcopy(H))
			return self._constr_nd_list(level-1,Hnew)
	
	def __getitem__(self,ind):
		if not hasattr(ind,'__iter__'): # if index is a single number instead of tuple/array, wrap it
			ind = [ind]
		temp = self.H
		for i in ind: # unpack the list progressively
			temp = temp[i]
		return temp # return the list element that has indices ind
	
	def __setitem__(self,ind,value):
		if not hasattr(ind,'__iter__'): # if index is a single number instead of tuple/array, wrap it
			ind = [ind]
		temp = self.H
		# unpack the list progressively, have to use a trick here, 
		# stop one earlier, because if list contains immutable types 
		# a copy will be made instead of using the reference and self.H
		# will not be modified, just a local copy in temp
		for i in range(len(ind)-1):
			temp = temp[ind[i]]
		temp[ind[-1]] = value
	
	def __str__(self):
		return self.H.__str__()



if __name__=="__main__":
	l = nd_list((5,5),0)
	print l
	print l[0]
	print l[0][0] # sequential indexing ok
	l[0,0] = 666  # tuple/array indexing ok
	print l[0]
	print l
	l[0,1] = [2,3,4] # you can assign whatever you want to the list items, for example an other nd_list :)
	print l
	print l[0][1]
	# the following will destroy the original shape of our nd_list, so we supply no shape function as one can destroy the original shape as follows
	l[0] = 1
	l[1] = 2
	l[2] = 3
	l[3] = 4
	l[4] = 5
	print l
