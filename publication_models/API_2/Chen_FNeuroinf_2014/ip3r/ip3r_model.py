import steps.interface

# IP3 receptor model

from steps.model import *
from steps.geom import *

###############################################################################

def getModel():
	mdl = Model()
	with mdl:
		
		# chemical species objects
		# Calcium
		# IP3
		
		# receptor state objects
		# IP3 receptor in 'naive' state
		# bound IP3 
		# bound IP3 and Ca (open)
		# 1 bound Ca to inactivation site
		# 2 bound Ca to inactivation sites
		# 3 bound Ca to inactivation sites
		
		# chemical species objects
		Ca, IP3, R, RIP3, Ropen, RCa, R2Ca, R3Ca, R4Ca = Species.Create()				# Calcium
		
		surfsys = SurfaceSystem.Create()
	with surfsys, mdl:
		
		# The 'forward' binding reactions: 
			
		# The 'backward' binding reactions:
		R.s + IP3.o <r['R_bind_IP3_f']> RIP3.s
		RIP3.s + Ca.o <r['RIP3_bind_Ca_f']> Ropen.s
		R.s + Ca.o <r['R_bind_Ca_f']> RCa.s
		RCa.s + Ca.o <r['RCa_bind_Ca_f']> R2Ca.s
		R2Ca.s + Ca.o <r['R2Ca_bind_Ca_f']> R3Ca.s
		R3Ca.s + Ca.o <r['R3Ca_bind_Ca_f']> R4Ca.s
		
		# Ca ions passing through open IP3R channel
		Ca.i + Ropen.s <r['R_Ca_channel_f']> Ropen.s + Ca.o
	
	# The reaction constants
	r['R_bind_IP3_f'].setRates(1000000000.0, 25800.0)
	r['RIP3_bind_Ca_f'].setRates(8000000000.0, 2000.0)
	r['R_bind_Ca_f'].setRates(8889000.0, 5.0)
	r['RCa_bind_Ca_f'].setRates(20000000.0, 10.0)
	r['R2Ca_bind_Ca_f'].setRates(40000000.0, 15.0)
	r['R3Ca_bind_Ca_f'].setRates(60000000.0, 20.0)
	
	# Corresponds to Ca input ~ 20000/ms for open receptor
	r['R_Ca_channel_f'].setRates(8000000.0, 8000000.0)          
	
	return mdl

###############################################################################

def getGeom():
	mesh = Geometry()
	with mesh:
		cyt = Compartment.Create(None, 1e-19)
		
		ER = Compartment.Create(None, 2e-20)					
		
		ERmemb = Patch.Create(ER, cyt, 'surfsys')
	
	return mesh

###############################################################################
