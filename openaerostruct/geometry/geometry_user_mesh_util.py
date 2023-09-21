'''
geomEngine3 straight port to Python for OAS
by Safa Bakhshi
WIP
'''

#For testing
import numpy as np
import matplotlib.pyplot as plt


def userGeom4(sections,data,bPanelsL,cPanels,plotCont,symmetry=True):
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

    #Data Initialization
	if np.shape(data) != (sections,4):
		raise Exception("Data mismatch")

	if len(bPanelsL) != sections:
		raise Exception("Define Panels")


	tquartX = []
	panelQC = []
	panelGY = []
	panelGX = []
	K = []

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
			tipC = rootC*taper;
			S = AR*((tipC+rootC)/2)**2
			b = np.sqrt(AR*S)
			leLambda = np.deg2rad(data[sec,3])
			Stot = Stot + S;
			span = span + b;
			if sec == sections-1:
				maintipC = tipC

		CAve = (maintipC + mainrootC)/2

		#Generate geomtery
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
		rootChordPoints = np.arange(rootStart,rootEnd-(rootStart-rootEnd)/cPanels,-(rootStart-rootEnd)/cPanels)
		rootChordPoints = np.linspace(rootStart,rootEnd,cPanels+1)

		if tipStart == tipEnd:
			tipChordPoints = tipStart*np.ones(len(rootChordPoints))
		else:
			tipChordPoints = np.linspace(tipStart,tipEnd,len(rootChordPoints))

		
		K.append(b/(2*bPanels))
		if sec == 0:
			panelGeomY = np.arange(-b/2,b/2+K[sec],K[sec])
		else:
			panelGeomY = np.zeros(2*bPanels)
			panelGeomY[0:bPanels] = np.linspace(panelGY[sec-1][0]-(b/2),panelGY[sec-1][0]-K[sec],bPanels)
			panelGeomY[bPanels:2*bPanels] = np.linspace(panelGY[sec-1][-1]+K[sec],panelGY[sec-1][-1]+(b/2),bPanels)

		if sec == 0:
			centreIndex = bPanels;
		else:
			centreIndex = bPanels -1;

		panelGeomX = np.zeros([len(rootChordPoints),len(panelGeomY)])
		if sec == 0:
			for i in range(len(rootChordPoints)):
				#Left Wing
				panelGeomX[i,0:centreIndex+1] = rootChordPoints[i] + ((tipChordPoints[i]-rootChordPoints[i])/(-b/2))*(panelGeomY[0:centreIndex+1])
				#Right Wing
				panelGeomX[i,centreIndex+1:] = rootChordPoints[i] + ((tipChordPoints[i]-rootChordPoints[i])/(b/2))*(panelGeomY[centreIndex+1:])
		else:
			for i in range(len(rootChordPoints)):
				#Left Wing
				panelGeomX[i,0:centreIndex+1] = rootChordPoints[i] + ((tipChordPoints[i]-rootChordPoints[i])/(-b/2))*((panelGeomY[0:centreIndex+1])-panelGY[sec-1][0])
				#Right Wing
				panelGeomX[i,centreIndex+1:] = rootChordPoints[i] + ((tipChordPoints[i]-rootChordPoints[i])/(b/2))*(panelGeomY[centreIndex+1:] - panelGY[sec-1][-1])
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

	#Stitch the results into a mesh
	for i in range(sections):
		if i == 0:
			mainQuarterC = panelQC[i]
			mainPanelGeomY = panelGY[i]
			mainPanelGeomX = panelGX[i]
			maintquarterPointsX = tquartX[i]
		else:
			mainQuarterC = np.concatenate((panelQC[i][:,0:bPanelsL[i]-1],mainQuarterC,panelQC[i][:,bPanelsL[i]:]),axis=1)
			maintquarterPointsX = np.concatenate((tquartX[i][:,0:bPanelsL[i]-1],maintquarterPointsX,tquartX[i][:,bPanelsL[i]:]),axis=1)
			mainPanelGeomY = np.concatenate((panelGY[i][0:bPanelsL[i]-1],mainPanelGeomY,panelGY[i][bPanelsL[i]:]))
			mainPanelGeomX = np.concatenate((panelGX[i][:,0:bPanelsL[i]-1],mainPanelGeomX,panelGX[i][:,bPanelsL[i]:]),axis=1)

	#Use stiched geomtery to calculate control point locations and chord distribution

	mainControlPointsX = np.zeros([cPanels,2*np.sum(bPanelsL)])
	mainVControlPointsX = np.zeros([cPanels,2*np.sum(bPanelsL)])
	mainControlPointsY = np.zeros(2*np.sum(bPanelsL))
	chordDistGeom = np.zeros(len(mainPanelGeomY))
	chordDistCont = np.zeros(len(mainControlPointsY))

	for i in range(len(mainPanelGeomY)):
		chordDistGeom[i] = mainPanelGeomX[0,i] - mainPanelGeomX[-1,i]

	for i in range(len(mainPanelGeomY)-1):
		for j in range(cPanels):
			mainControlPointsX[j,i] = (maintquarterPointsX[j,i+1]+maintquarterPointsX[j,i])/2
			mainVControlPointsX[j,i] = (mainQuarterC[j,i+1]+mainQuarterC[j,i])/2
		mainControlPointsY[i] = mainPanelGeomY[i] + (mainPanelGeomY[i+1]-mainPanelGeomY[i])/2
		chordDistCont[i] = (chordDistGeom[i+1] + chordDistGeom[i])/2

	#Cut the planform if symmetry is set

	return span, Stot, panelGX, panelGY


#TEST

sections = 4;
#data = np.array([[0.648294, .698608, .345206, .205671],[187.23, 121.37, 84.79, 29.27],[.075305,.089833,.573734,3.46289],[84.9578,82.3695,68,40]]);
data = np.array([[0.648294,187.23,.075305,84.9578],[.698608,121.37,.089833,82.3695],[.345206,84.79,.573734,68],[.205671,29.27,3.46289,40]])
bPanels = np.array([3,2,8,15]);
cPanels = 16;

#sections = 2
#bPanels = np.array([8,8])
#cPanels = 8

#data = np.array([[0.3,20,9,33]])
#data = np.array([[0.3,20,9,33],[0.3,20,9,50]])
[b, Stot, panelGX, panelGY] = userGeom4(sections,data,bPanels,cPanels,'true',True)
print('The area of this planform is: {}'.format(Stot))

#Dictionary Idea
Np3SST = {'sections': 4, 
		  'bPanels': [3,2,8,15], 
		  'cPanels': 16,
		  'taper': [0.648294,.698608,.345206,.205671],
		  'rootC': [187.23,121.37,84.79,29.27],
		  'LEsweep': [84.9578,82.3695,68,40],
		  'symmetry': True,
		  'plotSymmetry': 'Left'}



def plotPlanform(sections,panelGX,panelGY,plotSymmetry='Left'):
	"""
    The the outline of all definied planform sections.

    Parameters
    ----------
    sections : int
        Integer for number of wing section specified
    panelGX: list[numpy array]
    	Panel X-coordinate data
    panelGY: list[numpy array]
    	Panel Y-coordinate data
	symmetry : string
        Flag set to 'Left' if only y<=0 plotted
        Flag set to 'Right' if only y>=0 plotted

    Returns
    -------
    Figure
    """	
	np.random.seed(91)
	colorSet = np.random.rand(sections,3)
	plt.figure()
	for i in range(sections):
		if i == 0:
			
			rootEnd = panelGX[i][cPanels,bPanels[i]]
			rootStart = panelGX[i][0,bPanels[i]]
			#tipStart = panelGX[i][0,0]
			#tipEnd = panelGX[i][cPanels,0]
			
			if plotSymmetry == 'Left':
				plt.plot([ 0,0,panelGY[i][0],panelGY[i][0],0 ],[rootStart,rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart],c=colorSet[i,:])
			elif plotSymmetry == 'Right':
				plt.plot([ 0,panelGY[i][-1],panelGY[i][-1],0,0 ],[rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart,rootStart],c=colorSet[i,:])
			else:
				plt.plot([ 0,panelGY[i][0],panelGY[i][0],0 ],[rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart],c=colorSet[i,:])
				plt.plot([ 0,panelGY[i][-1],panelGY[i][-1],0],[rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart],c=colorSet[i,:])

		else:
			
			rootEnd = panelGX[i-1][cPanels,0]
			rootStart = panelGX[i-1][0,0]
			#tipStart = panelGX[i][0,0]
			#tipEnd = panelGX[i][cPanels,0]
			
			if plotSymmetry == 'Left':
				plt.plot([ panelGY[i-1][0],panelGY[i-1][0],panelGY[i][0],panelGY[i][0],panelGY[i-1][0] ],[ rootStart,rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart ],c=colorSet[i,:])
			elif plotSymmetry == 'Right':
				plt.plot([ panelGY[i-1][-1],panelGY[i][-1],panelGY[i][-1],panelGY[i-1][-1],panelGY[i-1][-1] ],[ rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart, rootStart ],c=colorSet[i,:])
			else:
				plt.plot([ panelGY[i-1][0],panelGY[i-1][0],panelGY[i][0],panelGY[i][0],panelGY[i-1][0] ],[ rootStart,rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart ],c=colorSet[i,:])
				plt.plot([ panelGY[i-1][-1],panelGY[i][-1],panelGY[i][-1],panelGY[i-1][-1],panelGY[i-1][-1] ],[ rootEnd,panelGX[i][cPanels,0],panelGX[i][0,0],rootStart, rootStart ],c=colorSet[i,:])

plotPlanform(sections,panelGX,panelGY,'Left')


def plotPanels(sections,panelGX,panelGY):
	np.random.seed(151)
	colorSet = np.random.rand(sections,3)
	plt.figure()
	for i in range(sections):
		if i == 0:
			for j in range(len(panelGX[i][:,bPanels[i]])):
				plt.plot([0,panelGY[i][0]],[panelGX[i][j,bPanels[i]],panelGX[i][j,0]],c=colorSet[i,:])
				#plt.plot([0,panelGY[i][-1]],[panelGX[i][j,bPanels[i]],panelGX[i][j,0]],c=colorSet[i,:])

			for k in range(int(len(panelGY[i])/2)):
				plt.plot([panelGY[i][k],panelGY[i][k]],[panelGX[i][0,k],panelGX[i][-1,k]],c=colorSet[i,:])		
		else:
			for j in range(len(panelGX[i-1][:,0])):
				plt.plot([panelGY[i-1][0],panelGY[i][0]],[panelGX[i-1][j,0],panelGX[i][j,0]],c=colorSet[i,:])
				#plt.plot([panelGY[i-1][-1],panelGY[i][-1]],[panelGX[i-1][j,0],panelGX[i][j,0]],c=colorSet[i,:])

			for k in range(int(len(panelGY[i])/2)):
				plt.plot([panelGY[i][k],panelGY[i][k]],[panelGX[i][0,k],panelGX[i][-1,k]],c=colorSet[i,:])

plotPanels(sections,panelGX,panelGY)

plt.show()