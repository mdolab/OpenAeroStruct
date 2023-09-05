'''
geomEngine3 port to Python for OAS
by Safa Bakhshi
WIP
'''

#For testing
import numpy as np
import matplotlib.pyplot as plt


def userGeom4(sections,data,bPanelsL,cPanels,plotCont,symmetry):
	"""
    Quickly generates user mesh for analysis based on a few parameters

    Parameters
    ----------
    sections : int
        Integer for number of wing section specified
    data : numpy array
        sectionsx4 array consisting of taper, root chord length, aspect ratio, and leading edge sweep in column. Each row corresponds to each specified wing section.
    bPanelsL: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section
	cPanels: int
        Integer for number of chord-wise panels
	symmetry : boolean
        Flag set to True if surface is reflected about y=0 plane.

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the aerodynamic surface.
    """	

	if np.shape(data) != (sections,4):
		raise Exception("Data mismatch")

	if len(bPanelsL) != sections:
		raise Exception("Define Panels")


	tquartX = []
	panelQC = []
	panelGY = []
	panelGX = []


	Stot = 0
	span = 0
	Npanels = 2*np.sum(bPanelsL)*cPanels

	for sec in range(sections):
		# Select panel number set
		bPanels = bPanelsL[sec]

		#GeomEngine Baseline

		if sec == 0:
			taper = data[sec,0]
			rootC = data[sec,1]
			AR = data[sec,2]
			tipC = rootC*taper
			S = AR*((tipC+rootC)/2)**2
			b = np.sqrt(AR*S)
			leLambda = np.deg2rad(data[sec,3])
			Stot = Stot + S
			span = span + b
			mainrootC = rootC
			maintipC = tipC
		else:
			taper = data[sec,0]
			rootC = np.abs(panelGX[sec-1][0,0] - panelGX[sec-1][cPanels,0])
			AR = data[sec,2]
			tipC = rootc*taper;
			S = AR*((tipC+rootC)/2)**2
			b = sqrt(AR*S)
			leLambda = np.deg2rad(data[sec,3])
			Stot = Stot + S;
			span = span + b;
			if sec == sections-1:
				maintipC = tipC

		CAve = (maintipC + mainrootC)/2


		if sec == 0:
			rootEnd = 0;
			rootStart = rootC + rootEnd

			tipStart = rootStart - (b/2)*np.tan(leLambda)
			tipEnd = tipStart - tipC
		else:
			rootEnd = panelGX[sec-1][cPanels,0]
			rootStart = rootC + rootEnd

			tipStart = rootStart - (b/2)*np.tan(leLambda)
			tipEnd = tipStart - tipC

		#Generate Panels
		#Descretize wing root and wing tip into N+1 points  
		rootChordPoints = np.arange(rootStart,rootEnd,-(rootStart-rootEnd)/cPanels)
		tipChordPoints = np.arange(tipStart,tipEnd,-(tipStart-tipEnd)/cPanels)

		if not tipChordPoints.any():
			tipChordPoints = np.arange(tipStart*np.ones(1,len(rootChordPoints)))

		K = []
		K.append(b/(2*bPanels))
		if sec == 0:
			panelGeomY = np.arange(-b/2,b/2,K[sec])
		else:
			panelGeomY = np.zeros(1,2*bPanels)
			panelGeomY[0:bPanels-1] = np.arange( panelGY[sec-1][0]-(b/2), panelGY[sec-1][1]-K[sec], K[sec]  )


		if sec == 0:
			centreIndex = bPanels;
		else:
			centreIndex = bPanels - 1;

		panelGeomX = np.zeros([len(rootChordPoints),len(panelGeomY)])
		if sec == 0:
			for i in range(len(rootChordPoints)):
				#Left Wing
				panelGeomX[i,0:centreIndex] = rootChordPoints[i] + ((tipChordPoints[i]-rootChordPoints[i])/(-b/2))*(panelGeomY[0:centreIndex])
				#Right Wing
				panelGeomX[i,centreIndex+1:] = rootChordPoints[i] + ((tipChordPoints[i]-rootChordPoints[i])/(b/2))*(panelGeomY[centreIndex+1:])

		panelQuarterC = np.zeros([cPanels,len(panelGeomY)])
		tquarterPointsX = np.zeros([cPanels,len(panelGeomY)])

		for i in range(len(panelGeomY)):
			for j in range(cPanels-1):
				panelQuarterC[j,i] = panelGeomX[j,i] + (panelGeomX[j+1,i] - panelGeomX[j,i])/4
				tquarterPointsX[j,i] = panelGeomX[j,i] + 1*(panelGeomX[j+1,i]-panelGeomX[j,i])/2

		panelQC.append(panelQuarterC)
		panelGY.append(panelGeomY)
		tquartX.append(tquarterPointsX)
		panelGX.append(panelGeomX)

	#Stitch the results

	for i in range(sections):
		if i == 0:
			mainQuarterC = panelQC[i]
			mainPanelGeomY = panelGY[i]
			mainPanelGeomX = panelGX[i]
			maintquarterPointsX = tquartX[i]
		else:
			mainQuarterC = np.concatenate((panelQC[i][:,0:bPanelsL[i]-1],mainQuarterC,panelQC[i][:,bPanelsL[i]:]),axis=1)
			maintquarterPointsX = np.concatenate((tquartX[i][:,0:bPanelsL[i]-1],maintquarterPointsX,tquartX[i][:,bPanelsL[i]:]),axis=1)
			mainPanelGeomX = np.concatenate((panelGY[i][:,0:bPanelsL[i]-1],mainPanelGeomY,panelGY[i][:,bPanelsL[i]:]),axis=1)
			mainPanelGeomY = np.concatenate((panelGX[i][:,0:bPanelsL[i]-1],mainPanelGeomX,panelGX[i][:,bPanelsL[i]:]),axis=1)


	return span, rootStart, rootEnd, tipEnd, tipStart, rootStart, tipStart, tipEnd, rootEnd



#TEST
sections = 1
bPanels = np.array([8])
cPanels = 8

data = np.array([[1,1,10,0]])
[b, rootStart, rootEnd, tipEnd, tipStart, rootStart, tipStart, tipEnd, rootEnd] = userGeom4(sections,data,bPanels,cPanels,'true',True)


plt.figure(1)
plt.plot(np.array([0,0,-b/2,-b/2,0,0,b/2,b/2,0]),np.array([rootStart,rootEnd,tipEnd,tipStart,rootStart,rootStart,tipStart,tipEnd,rootEnd ]))
plt.show()
