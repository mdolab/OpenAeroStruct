#VORTEX for Python Modernized
import numpy as np
def VORTEX(G, R):
    TOL  = 1e-10

    if len(G) != 3:
         raise Exception('G must be a 1x3 vector. G was size {}'.format(len(G)))
    if len(R) != 3:
         raise Exception('R must be a 1x3 vector. R was size {}'.format(len(G)))

    #R = -1*R #Kroo convention: Vortex Root minus Control Point
    Gsq = np.sum(G**2)
    Rsq = np.sum(R**2) 
    R2 = R-G
    RXR2 = np.cross(R,R2)
    GXR = np.cross(G,R)
    sig1 = (np.dot(G,R))/np.linalg.norm(R)
    sig2 = (np.dot(G,R2))/np.linalg.norm(R2)
    
    E1 = np.sum(GXR**2)
    E2 = np.sum((G+R)**2)

    #Set equal to 0 if within TOL 
    if ((Rsq <= Gsq*TOL) or (E2 <= Gsq*TOL) or (E1 <= Gsq*TOL*Rsq)): #Check for colinearity (core model)
        return np.zeros(3)
    else:
        V = (sig1-sig2) / (4.0*np.pi)
        return V*np.array([(RXR2[0]/np.sum(RXR2**2)),(RXR2[1]/np.sum(RXR2**2)),(RXR2[2]/np.sum(RXR2**2))])
    


def VORTEX2(R1, R2):
    TOL  = 1e-4

    if len(R1) != 3:
         raise Exception('R1 must be a 1x3 vector. R1 was size {}'.format(len(R1)))
    if len(R2) != 3:
         raise Exception('R2 must be a 1x3 vector. R2 was size {}'.format(len(R2)))

    #R = R #Kroo convention Vortex Root minus Control Point

    G = R1 - R2
    Gsq = np.sum(G**2)
    R1sq = np.sum(R1**2) 
    GXR1 = np.cross(G,R1)

    R1XR2 = np.cross(R1,R2)
    R1XR22 = np.linalg.norm(R1XR2)**2
    kfid = np.dot((R1-R2),(R1/np.linalg.norm(R1) - R2/np.linalg.norm(R2)))
    
    E1 = np.sum(GXR1**2)
    E2 = np.sum((G+R1)**2)

    if ((R1sq <= Gsq*TOL) or (E2 <= Gsq*TOL) or (E1 <= Gsq*TOL*R1sq)): #Check for colinearity (core model)
        return np.zeros(3)
    else:
        V = (kfid) / (4.0*np.pi)
        return V*(R1XR2/R1XR22)