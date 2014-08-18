import ctypes


lib = ctypes.cdll.LoadLibrary("/home/ccolle/2bodymom/src/MeanFieldNucleusWrapper.so")
# return types
lib.getA.restype             = ctypes.c_int
lib.getZ.restype             = ctypes.c_int
lib.getPLevels.restype       = ctypes.c_int
lib.getNLevels.restype       = ctypes.c_int
lib.getTotalLevels.restype   = ctypes.c_int
lib.getMassA.restype         = ctypes.c_double
lib.getMassA_min_pp.restype  = ctypes.c_double
lib.getMassA_min_pn.restype  = ctypes.c_double
lib.getMassA_min_nn.restype  = ctypes.c_double
lib.getFinalMProton.restype  = ctypes.c_int
lib.getFinalMNeutron.restype = ctypes.c_int
# argument types
lib.getA.argtype             = ctypes.c_void_p
lib.getZ.argtype             = ctypes.c_void_p
lib.getPLevels.argtype       = ctypes.c_void_p
lib.getNLevels.argtype       = ctypes.c_void_p
lib.getTotalLevels.argtype   = ctypes.c_void_p
lib.getExcitation.argtype    = ctypes.c_void_p
lib.getMassA.argtype         = ctypes.c_void_p
lib.getMassA_min_pp.argtype  = ctypes.c_void_p
lib.getMassA_min_pn.argtype  = ctypes.c_void_p
lib.getMassA_min_nn.argtype  = ctypes.c_void_p

class MeanFieldNucleus:
	s = ["He","C","O","Fe","Pb","Al","Cu","Au"]
	names = { b : a for a,b in enumerate(s) } # acts enum like, now you can construct using MeanFieldNucleus.names["C"]
	
	def __init__(self,n):
		self.nuc = lib.get_MeanFieldNucleus(n) # [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
		self.openPtr = True
	
	## because destructor __del__ can be called after global variable lib has already dissappeared 
	#  make sure you call the __del__ destructor manually before Python's garbage collector comes along.
	# 
	# I think this is more elegant then reloading a seperate lib instance for each
	# class object which would resolve the need to manually calling __del__ at the cost of probably much more memory
	# . (something doing like self.lib = ctypes.cdll.LoadLibrary("..."); self.lib.function.restype = ....;)
	#
	# don't forget to pass the pointer to nuc to EVERY function!
	
	def __del__(self):
		if lib: # library reference still exists
			lib.del_MeanFieldNucleus(self.nuc)
		else:   #oops probably leaked some memory
			print('Warning! The shared library was already removed from memory when the MeanFieldNucleus destructor ' \
				'was called. Memory leaks unavoidable! Make sure you manually delete all instances of MeanFieldNucleus '\
				'in stead of letting the Python garbage collector doing it!')

	def getA(self):
		return lib.getA(self.nuc)
	
	def getZ(self):
		return lib.getZ(self.nuc)

	def getN(self):
		return self.getA()-self.getZ()

	def getPLevels(self):
		return lib.getPLevels(self.nuc)
	
	def getNLevels(self):
		return lib.getNLevels(self.nuc)
	
	def getTotalLevels(self):
		return lib.getTotalLevels(self.nuc)
	
	def getMassA(self):
		return lib.getMassA(self.nuc)
	
	def getMassA_min_pp(self):
		return lib.getMassA_min_pp(self.nuc)
	
	def getMassA_min_pn(self):
		return lib.getMassA_min_pn(self.nuc)
	
	def getMassA_min_nn(self):
		return lib.getMassA_min_nn(self.nuc)

	def getFinalMProton(self):
		return lib.getFinalMProton(self.nuc)

	def getFinalMNeutron(self):
		return lib.getFinalMNeutron(self.nuc)

	def getN_array(self):
		lib.getN_array.restype = ctypes.POINTER(ctypes.c_int * lib.getTotalLevels(self.nuc) )
		res = lib.getN_array(self.nuc)
		return [ i for i in res.contents ]
	
	def getL_array(self):
		lib.getL_array.restype = ctypes.POINTER(ctypes.c_int * lib.getTotalLevels(self.nuc) )
		res = lib.getL_array(self.nuc)
		return [ i for i in res.contents ]
	
	def getJ_array(self):
		lib.getJ_array.restype = ctypes.POINTER(ctypes.c_int * lib.getTotalLevels(self.nuc) )
		res = lib.getJ_array(self.nuc)
		return [ i for i in res.contents ]
	
	def getExcitation(self):
		lib.getExcitation.restype = ctypes.POINTER(ctypes.c_double * lib.getTotalLevels(self.nuc) )
		res = lib.getExcitation(self.nuc)
		return [ i for i in res.contents ]
	"""
	def getShellOccupation(self,shell): # get the number of particles in a certain shell
		if shell < self.getPLevels(): # proton
			n = self.getZ()
			for s=0
		elif self.getPLevels() <= shell < self.getTotalLevels():
			n = self.getN()
		else:
			raise Exception("Invalid shell combo")
	"""
if __name__=="__main__":
	#print MeanFieldNucleus.names["Pb"]
	n = MeanFieldNucleus(MeanFieldNucleus.names["Pb"])
	print("getPLevels()      : {:d}".format(n.getPLevels()))
	print("getTotalLevels()  : {:d}".format(n.getTotalLevels()))
	print("Shell combinations: {:d}".format((n.getPLevels()+1)*n.getPLevels()/2))
	print n.getMassA()
	print n.getMassA_min_nn()
	print n.getExcitation()
	print n.getA()
	print n.getZ()
	print n.getN_array()
	print n.getL_array()
	print n.getJ_array()
	print n.getFinalMProton()
	print n.getFinalMNeutron()
	del n
