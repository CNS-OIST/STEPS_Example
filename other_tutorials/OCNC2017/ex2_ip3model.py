# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Example 2: Surface reactions: IP3 receptor model 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import steps.model as smodel
import steps.geom as sgeom

###############################################################################

def getModel():

	# Create model container object
	mdl = smodel.Model()
	
	# chemical species objects
	Ca = smodel.Spec('Ca', mdl)				# Calcium
	IP3 = smodel.Spec('IP3', mdl)			# IP3
	
	# receptor state objects
	R = smodel.Spec('R', mdl)				# IP3 receptor in 'naive' state
	RIP3 = smodel.Spec('RIP3', mdl)			# bound IP3 
	Ropen = smodel.Spec('Ropen', mdl)		# bound IP3 and Ca (open)
	RCa = smodel.Spec('RCa', mdl)			# 1 bound Ca to inactivation site
	R2Ca = smodel.Spec('R2Ca', mdl)			# 2 bound Ca to inactivation sites
	R3Ca = smodel.Spec('R3Ca', mdl)			# 3 bound Ca to inactivation sites
	R4Ca = smodel.Spec('R4Ca', mdl)			# 4 bound Ca to inactivation sites
	
	surfsys = smodel.Surfsys('ssys', mdl)
	
	# The 'forward' binding reactions: 
	R_bind_IP3_f = smodel.SReac('R_bind_IP3_f', surfsys, \
		olhs=[IP3], slhs=[R], srhs=[RIP3])
	RIP3_bind_Ca_f = smodel.SReac('RIP3_bind_Ca_f', surfsys, \
		olhs=[Ca], slhs=[RIP3], srhs = [Ropen])
	R_bind_Ca_f = smodel.SReac('R_bind_Ca_f', surfsys, \
		olhs=[Ca], slhs=[R], srhs=[RCa])
	RCa_bind_Ca_f = smodel.SReac('RCa_bind_Ca_f', surfsys, \
		olhs=[Ca], slhs=[RCa],srhs = [R2Ca])
	R2Ca_bind_Ca_f = smodel.SReac('R2Ca_bind_Ca_f', surfsys, \
		olhs=[Ca], slhs= [R2Ca], srhs = [R3Ca])
	R3Ca_bind_Ca_f = smodel.SReac('R3Ca_bind_ca_f', surfsys, \
		olhs=[Ca], slhs=[R3Ca], srhs=[R4Ca])
		
	# The 'backward' binding reactions:
	R_bind_IP3_b = smodel.SReac('R_bind_IP3_b', surfsys, \
		slhs=[RIP3], orhs=[IP3], srhs=[R])
	RIP3_bind_Ca_b = smodel.SReac('RIP3_bind_Ca_b', surfsys, \
		slhs=[Ropen], orhs=[Ca], srhs=[RIP3])
	R_bind_Ca_b = smodel.SReac('R_bind_Ca_b', surfsys, \
		slhs=[RCa], orhs=[Ca], srhs=[R])
	RCa_bind_Ca_b = smodel.SReac('RCa_bind_Ca_b', surfsys, \
		slhs=[R2Ca], orhs=[Ca], srhs=[RCa])
	R2Ca_bind_Ca_b = smodel.SReac('R2Ca_bind_Ca_b', surfsys, \
		slhs=[R3Ca], orhs=[Ca], srhs= [R2Ca])
	R3Ca_bind_Ca_b = smodel.SReac('R3Ca_bind_ca_b', surfsys, \
		slhs=[R4Ca], orhs=[Ca], srhs=[R3Ca])
	
	# Ca ions passing through open IP3R channel
	R_Ca_channel_f = smodel.SReac('R_Ca_channel_f', surfsys, \
		ilhs=[Ca], slhs=[Ropen], orhs=[Ca], srhs=[Ropen])
	R_Ca_channel_b = smodel.SReac('R_Ca_channel_b', surfsys, \
		olhs=[Ca], slhs=[Ropen], irhs=[Ca], srhs=[Ropen])
	
	# The reaction constants
	R_bind_IP3_f.setKcst(1000e6)
	R_bind_IP3_b.setKcst(25800)
	RIP3_bind_Ca_f.setKcst(8000e6)
	RIP3_bind_Ca_b.setKcst(2000)
	R_bind_Ca_f.setKcst(8.889e6)
	R_bind_Ca_b.setKcst(5)
	RCa_bind_Ca_f.setKcst(20e6)
	RCa_bind_Ca_b.setKcst(10)
	R2Ca_bind_Ca_f.setKcst(40e6)
	R2Ca_bind_Ca_b.setKcst(15)
	R3Ca_bind_Ca_f.setKcst(60e6)
	R3Ca_bind_Ca_b.setKcst(20)
	
	# Corresponds to Ca input ~ 20000/ms for open receptor
	R_Ca_channel_f.setKcst(8e6)          
	R_Ca_channel_b.setKcst(8e6)           
	
	# return model container object
	return mdl

###############################################################################

def getGeom():

    # Create geometry container object
    geom  = sgeom.Geom()
    
    # Create cytosol compartment
    cyt = sgeom.Comp('cyt', geom)
    # Assign volume to cytosol 
    cyt.vol = 0.1e-18
    
    # Create ER compartment
    ER = sgeom.Comp('ER', geom)					
    # Assign volume to ER 
    ER.setVol(0.02e-18)
    
    # Create the ER membrane
    ERmemb = sgeom.Patch('memb', geom, ER, cyt)
    
    # Add group of surface reaction rules to the ER membrane
    ERmemb.addSurfsys('ssys')
    
    # return geometry container object
    return geom

###############################################################################
