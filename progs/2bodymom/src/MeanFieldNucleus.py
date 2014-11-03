import ctypes


lib = ctypes.cdll.LoadLibrary("/home/camille/Code/trunk/progs/2bodymom/src/MeanFieldNucleusWrapper.so")
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
lib.getRange.restype         = ctypes.c_double
lib.getWF_r_step.restype     = ctypes.c_double
lib.getFinalMProton.restype  = ctypes.c_int
lib.getFinalMNeutron.restype = ctypes.c_int
lib.getTotalDensity.restype  = ctypes.c_double
lib.getProtonDensity.restype = ctypes.c_double
lib.getNeutronDensity.restype= ctypes.c_double
lib.getWaveFunction.restype  = None # is C equiv of void function
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
lib.getRange.argtype         = ctypes.c_void_p
lib.getWF_r_step.argtype     = ctypes.c_void_p
lib.getTotalDensity.argtypes  = [ctypes.c_void_p,ctypes.c_double]
lib.getProtonDensity.argtypes = [ctypes.c_void_p,ctypes.c_double]
lib.getNeutronDensity.argtypes= [ctypes.c_void_p,ctypes.c_double]
lib.getWaveFunction.argtypes  = [ctypes.c_void_p,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double] # times four because of 4 component wavefunctions!

class MeanFieldNucleus:
	s = ["He","C","O","Fe","Pb","Al","Cu","Au"]
	names = { b : a for a,b in enumerate(s) } # acts enum like, now you can construct using MeanFieldNucleus.names["C"]
	
	def __init__(self,n):
		self.nuc = lib.get_MeanFieldNucleus(n) # [0-7] =[He,C,O,Fe,Pb,Al,Cu,Au]
	
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

        def getRange(self):
                return lib.getRange(self.nuc)

	def getFinalMProton(self):
		return lib.getFinalMProton(self.nuc)

	def getFinalMNeutron(self):
		return lib.getFinalMNeutron(self.nuc)
        
        def getWF_r_step(self):
                return lib.getWF_r_step(self.nuc)

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
        def getWaveFunction(self,shellindex,m,r,costheta,phi):
                real = (ctypes.c_double * 4 )()
                imag = (ctypes.c_double * 4 )()
                lib.getWaveFunction(self.nuc,real,imag,shellindex,m,r,costheta,phi)
                return map(lambda x: complex(*x),zip(real,imag))
                

class MeanFieldNucleusThick(MeanFieldNucleus):
        s = ["He","C","O","Fe","Pb","Al","Cu","Au"]
	names = { b : a for a,b in enumerate(s) } # acts enum like, now you can construct using MeanFieldNucleus.names["C"]
        
        def __init__(self,n):
            self.nuc = lib.get_MeanFieldNucleusThick(n)

        def __del__(self):
                if lib:
                    lib.del_MeanFieldNucleusThick(self.nuc)
		else:   #oops probably leaked some memory
			print('Warning! The shared library was already removed from memory when the MeanFieldNucleus destructor ' \
				'was called. Memory leaks unavoidable! Make sure you manually delete all instances of MeanFieldNucleus '\
				'in stead of letting the Python garbage collector doing it!')
        def getTotalDensity(self,r):
            return lib.getTotalDensity(self.nuc,r)

    	def getTotalDensity_raw(self,r):
	    if (r < self.getWF_r_step()): # prevent coming too close to zero
		    r = self.getWF_r_step()
	    return self.getTotalDensity(r)/r/r

        def getProtonDensity(self,r):
            return lib.getProtonDensity(self.nuc,r)
        
        def getNeutronDensity(self,r):
            return lib.getNeutronDensity(self.nuc,r)

if __name__=="__main__":
	#print MeanFieldNucleus.names["Pb"]
	n = MeanFieldNucleus(MeanFieldNucleus.names["Fe"])
	print("#getPLevels()      : {:d}".format(n.getPLevels()))
	print("#getTotalLevels()  : {:d}".format(n.getTotalLevels()))
	print("#Shell combinations: {:d}".format((n.getPLevels()+1)*n.getPLevels()/2))
	print '#',n.getMassA()
	print '#',n.getMassA_min_nn()
	print '#',n.getExcitation()
	print '#',n.getA()
	print '#',n.getZ()
        print '#         : ',range(0,len(n.getN_array()))
        print '# n array : ',n.getN_array()
        print '# l array : ',n.getL_array()
        print '# j array : ',n.getJ_array()
	print '#',n.getFinalMProton()
	print '#',n.getFinalMNeutron()
        print n.getWaveFunction(1,1,1.,0.,0.)
	del n
        
        nthick = MeanFieldNucleusThick(MeanFieldNucleus.names["Pb"])
        import numpy as np
        import pylab as pl
        import scipy.integrate
        print "integral of density is ", 4.*np.pi*scipy.integrate.quad(nthick.getTotalDensity,0,nthick.getRange(),limit=250)[0]
        rr = np.linspace(nthick.getWF_r_step(),nthick.getRange(),100)
        pl.plot(rr,[ nthick.getTotalDensity(r)/r/r for r in rr ])
        pl.show()
        del nthick
        
