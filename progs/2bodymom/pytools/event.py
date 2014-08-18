import numpy as np
# just a container class for cleaner & shorter code
# should characterize the kinematics of a 2 nucleon 
# knockout event.
class Event:
	def __init__(self,id=None,xB=None,Q2=None,omega=None,q=None,shellindex1=None,shellindex2=None,k1=None,p1=None,p2=None):
		self.id             = id                # [ scalar ]
		self.xB             = xB		# [ scalar ]
		self.Q2             = Q2		# [ scalar ]
		self.omega          = omega		# [ scalar ]
		self.q              = q			# [ scalar ]
		self.shellindex1    = shellindex1	# [ scalar ] 
		self.shellindex2    = shellindex2	# [ scalar ] 
		self.k1             = np.array(k1)	# [3-vector] 
		self.p1             = np.array(p1)	# [3-vector] 
		self.p2             = np.array(p2)	# [3-vector] 

def parse_column_info(line):
	verts = line.split() # split on whitespace
	fields = []
	i=0
	for v in verts:
		field = v.split("[")[0]
		#assert(field in vars(Event())) # make sure member is present in Event class
		if not (field in vars(Event())):
			raise Exception("Event has no member named {:s}.".format(field))
		n = int(v.split("[")[1].rstrip("]"))
		if n==1:
			fields.append((field,i))
		else:
			fields.append((field,i,i+n))
		i+=n
	return fields

def parseEventsGroupShells(fname):
	f = file(fname,"r")
	events = []
	fields = None
	for line in f:
		line = line.strip() # strip whitespace
		if line.startswith('#::') and line.endswith('::#'): # column info
			fields = parse_column_info(line[3:-3]) # pass line without special identifiers #:: and ::#

		elif not ( line.startswith('#') or line==''): # skip all comment and empty lines
			if (fields==None): # we encountered data line without a column info line, we don't now how to parse!
				raise Exception("Data File Format Error! Data line encountered without column info line!")
			# here we have a data line
			cols = line.strip().split()       # split the line on whitespace
			event = Event()              # make blank event
			
			for field in fields: # loop over the different fields to set (members of Event() to set) 
				if len(field)==2: # scalar
					setattr(event,field[0],float(cols[field[1]]))
				elif len(field)==3: # vector
					setattr(event,field[0],np.array( map(float,cols[field[1]:field[2]]) ) )
				else:
					print("Unsupported field length :{:d}".format(len(field)))
					exit(-1)
			event.Pcm = event.k1 + event.p2 # you can do this, even though Pcm is not a declared member of event! Whoa! right?
			if event.id != events[-1][-1].id  or len(events)==0: # if id's differ of events is empty start new list of shellcombos
					events.append([])
			events[-1].append(event)
	return events

# group the events
# put the events with the same id
# in the same sublist
def groupShells(events):
	gEvents = []
	for ev in events:
		if len(gEvents)==0 or ev.id != gEvents[-1][-1].id: # start a new shell group when len is zero or we have new shellIndep id
			gEvents.append([])
		gEvents[-1].append(ev)
	return gEvents

def parseEvents(fname):
	f = file(fname,"r")
	events = []
	fields = None
	for line in f:
		line = line.strip() # strip whitespace
		if line.startswith('#::') and line.endswith('::#'): # column info
			fields = parse_column_info(line[3:-3]) # pass line without special identifiers #:: and ::#

		elif not ( line.startswith('#') or line==''): # skip all comment and empty lines
			if (fields==None): # we encountered data line without a column info line, we don't now how to parse!
				raise Exception("Data File Format Error! Data line encountered without column info line!")
			# here we have a data line
			cols = line.strip().split()       # split the line on whitespace
			event = Event()              # make blank event
			
			for field in fields: # loop over the different fields to set (members of Event() to set) 
				if len(field)==2: # scalar
					setattr(event,field[0],float(cols[field[1]]))
				elif len(field)==3: # vector
					setattr(event,field[0],np.array( map(float,cols[field[1]:field[2]]) ) )
				else:
					print("Unsupported field length :{:d}".format(len(field)))
					exit(-1)
			event.Pcm = event.k1 + event.p2 # you can do this, even though Pcm is not a declared member of event! Whoa! right?
			events.append(event)
	return events

