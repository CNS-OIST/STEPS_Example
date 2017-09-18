import math	

# Return all the tetrahedra with their barycenter inside a cylinder, centred on z axis.
# params are 1) steps.Tetmesh object, 2) radius of the cylinder, 
#			 3) the x coordinate of the centre of the cylinder (check if this is 0.0 in cubit)
#			 4) the y coordinate of the centre of the cylinder (again check output from cubit)
def getcyl(tetmesh, rad,  zmin, zmax, binnum=120, x = 0.0, y = 0.0):
	rad = float(rad)
	ntets = tetmesh.countTets()
	keeptets = []
	for t in range(ntets):
		baryc = tetmesh.getTetBarycenter(t)
		if (baryc[0] < x-rad or baryc[0] > x+rad): continue
		if (baryc[1] < y-rad or baryc[1] > y+rad) : continue
		if (baryc[2] < zmin or baryc[2] > zmax): continue
		if (math.sqrt(math.pow(baryc[0]-x, 2) + math.pow(baryc[1]-y,2)) > rad) : continue
		keeptets.append(t)
	binmin = tetmesh.getBoundMin()[2]
	binmax = tetmesh.getBoundMax()[2]
	binlength = (binmax-binmin)/binnum
	binned = []
	for b in range(binnum):
		binned.append([])
	assert(binned.__len__() == binnum)
	for t in keeptets:
		zcoord = binmin
		tz = tetmesh.getTetBarycenter(t)[2]
		for b in range(binnum):
			znew = zcoord+binlength
			if (tz >= zcoord and tz < znew): 
				binned[b].append(t)
				break
			zcoord = znew
	
	
	return (keeptets, binned, binmin, binmax, binlength)
	
											