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