""" Utility for quickly generating a multi-section user specificed OAS mesh
"""

import numpy as np
import matplotlib.pyplot as plt


def generateMesh(
    surface
):
    """
    Quickly generates user multi-section mesh based on a few parameters

    Parameters
    ----------
    sections : int
        Integer for number of wing sections specified
    data : numpy array
        sectionsx4 array consisting of taper, root chord length, aspect ratio, and leading edge sweep in columns. Each row corresponds to each specified wing section.
    bPanels: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section
    cPanels: int
        Integer for number of chord-wise panels
    symmetry : boolean
        Flag set to True if surface is reflected about y=0 plane.
    controlPoints: boolean
        Flag set to True if both quarter chord and three-quarter chord control points should be computed
    generatePlots: boolean
        Flag set to True if planform and mesh plots should be created and saved.
    plotOptions: Dict
        Dictionary containing user options for plotting
            'type': 'plan', 'mesh', 'all'
                Plots the planform geometry only, mesh only, or both

            'symmetry': 'Left', 'Right', 'Full
                Plots the left, right, or full spans

            'name': string
                Name of wing

    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the aerodynamic surface.
    """

    #Parse Data
    numSections = surface["num_sections"]
    symmetry = surface["symmetry"]


    if symmetry or numSections == 1:
        rootSection = numSections-1
    else:
        if "rootSection" not in surface.keys():
            raise Exception("The root section of an asymmetrical mesh needs to be identified")
        else:
            rootSection = surface["rootSection"]
            

    sectionData = {'taper':surface["sec_taper"],'sweep':surface["sec_sweep"],'span':surface["sec_span"],'rootChord':surface["root_chord"]}

    if 'bPanels' in surface.keys():
        ny = surface['bPanels'] + 1
    else:
        ny = surface['sec_ny']

    if 'cPanels' in surface.keys():
        nx = surface['cPanels'] + 1
    else:
        nx = surface['nx']

    panelGX, panelGY  = generateSectionGeometry(numSections, symmetry, sectionData, ny, nx, rootSection)
    panelGeomX, panelGeomY = stitchSectionGeometry(numSections, panelGY, panelGX)
    mesh = outputOASMesh(panelGeomX, panelGeomY)

    #Produce meshes for each section
    sec_meshes = []
    for section in range(numSections):
        secMesh = outputOASMesh(panelGX[section],panelGY[section])
        sec_meshes.append(secMesh)


    return mesh, sec_meshes


def generateSectionGeometry(sections, symmetry, sectionData, ny, nx, rootSection):
    """
    Constructs the multi-section wing geometry specified by the user and generates a mesh for each section.

    Parameters
    ----------
    sections : int
        Integer for number of wing sections specified
    data : numpy array
        sectionsx4 array consisting of taper, root chord length, aspect ratio, and leading edge sweep in each column. Each row corresponds to each specified wing section.
    bPanels: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section
    cPanels: int
        Integer for number of chord-wise panels


    Returns
    -------
    panelGX : List
         List containing the mesh x-coordinates for each section
    panelGY : List
        List containing the mesh y-coordinates for each section
    panelQC : List
         List containing the quarter chord x-positions along each panel edge. Not used by OAS but information is too good to get rid of at this point.
    tquartX : : List
         List containing the three-quarter chord x-positions along each panel edge. Not used by OAS but information is too good to get rid of at this point.
    Npanels : int
        Total number of panels in the mesh
    """

    panelGY = [None] * sections
    panelGX = [None] * sections
    Stot = 0

    
    #Jump to root section and build left wing
    for sec in np.arange(rootSection,-1,-1):
        taper = sectionData['taper'][sec]
        b = sectionData['span'][sec]
        leLambda = sectionData['sweep'][sec]

        if sec == rootSection:
            rootC = sectionData['rootChord']
            rootTE = 0
            rootY = 0
        else:
            rootC = np.abs(panelGX[sec + 1][0, 0] - panelGX[sec + 1][nx-1, 0])
            rootTE = panelGX[sec + 1][nx-1, 0]
            rootY = panelGY[sec + 1][0]
        tipC = rootC*taper
        Stot = Stot + b * ((tipC + rootC) / 2)

        rootLE = rootC + rootTE
        tipLE = rootLE - b * np.tan(leLambda)
        tipTE = tipLE - tipC

        rootX = np.linspace(rootLE, rootTE, nx)
        
        if tipLE == tipTE:
            tipX = tipLE * np.ones(len(rootX))
        else:
            tipX = np.linspace(tipLE, tipTE, len(rootX))


        panelGeomY = np.zeros(ny[sec])
        panelGeomX = np.zeros([nx, ny[sec]])

        panelGeomY = np.linspace(rootY - b, rootY, ny[sec])

        for i in range(len(rootX)):
            panelGeomX[i, :] = rootX[i] - ((tipX[i] - rootX[i]) / b) * (panelGeomY - rootY)

        panelGY[sec] = panelGeomY
        panelGX[sec] = panelGeomX


    #Build the right wing if asymmetrical
    if not symmetry:
        for sec in np.arange(rootSection+1,sections):
            taper = sectionData['taper'][sec]
            b = sectionData['span'][sec]
            leLambda = sectionData['sweep'][sec]

            rootC = np.abs(panelGX[sec - 1][0, -1] - panelGX[sec - 1][nx-1, -1])
            tipC = rootC*taper
            rootY = panelGY[sec - 1][-1]
            Stot = Stot + b * ((tipC + rootC) / 2)

            rootTE = panelGX[sec - 1][nx-1, 0]
            rootLE = rootC + rootTE
            tipLE = rootLE + b * np.tan(leLambda)
            tipTE = tipLE - tipC

            rootX = np.linspace(rootLE, rootTE, nx)
            
            if tipLE == tipTE:
                tipX = tipLE * np.ones(len(rootX))
            else:
                tipX = np.linspace(tipLE, tipTE, len(rootX))


            panelGeomY = np.zeros(ny[sec])
            panelGeomX = np.zeros([nx, ny[sec]])

            panelGeomY[0:ny[sec]] = np.linspace(rootY, rootY + b, ny[sec])

            for i in range(len(rootX)):
                panelGeomX[i, :] = rootX[i] + ((tipX[i] - rootX[i]) / (b / 2)) * (panelGeomY - rootY)

        panelGY[sec] = panelGeomY
        panelGX[sec] = panelGeomX


    return panelGX, panelGY


def stitchSectionGeometry(sections, panelGY, panelGX):
    """
    Combines the split section array into singular unified mesh

    Parameters
    ----------
    sections : int
        Integer for number of wing sections specified
    panelGX : List
         List containing the mesh x-coordinates for each section
    panelGY : List
        List containing the mesh y-coordinates for each section
    bPanels: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section


    Returns
    -------
    panelGeomX : numpy array
         Array of the mesh x-coordinates
    panelGeomY : numpy array
        Array of the mesh y-coordinates
    """
    # Stitch the results into a singular mesh

    if sections > 1:
        panelGeomY = panelGY[0][:-1]
        panelGeomX = panelGX[0][:,:-1]
        for i in np.arange(1,sections-1):
            panelGeomY = np.concatenate((panelGeomY,panelGY[i][:-1]))
            panelGeomX = np.concatenate((panelGeomX,panelGX[i][:, :-1]), axis=1)
        panelGeomY = np.concatenate((panelGeomY,panelGY[sections-1]))
        panelGeomX = np.concatenate((panelGeomX,panelGX[sections-1]), axis=1)
    else:
        panelGeomY = panelGY[0]
        panelGeomX = panelGX[0]
    return panelGeomX, panelGeomY


def reflectSymmetric(panelGeomX, panelGeomY):
    """
    Reflect the mesh over y=0 is symmetry conditions are used

    Parameters
    ----------
    panelGeomX : numpy array
         Array of the mesh x-coordinates
    panelGeomY : numpy array
        Array of the mesh y-coordinates

    Returns
    -------
    panelGeomX : numpy array
         Array of the mesh x-coordinates
    panelGeomY : numpy array
        Array of the mesh y-coordinates
    """
    panelGeomX = np.hstack((panelGeomX,panelGeomX))
    panelGeomY = np.hstack((panelGeomY,-panelGeomY))
    return panelGeomX, panelGeomY


def outputOASMesh(panelGeomX, panelGeomY):
    """
    Outputs the mesh in OAS format

    Parameters
    ----------
    panelGeomX : numpy array
         Array of the mesh x-coordinates
    panelGeomY : numpy array
        Array of the mesh y-coordinates



    Returns
    -------
    mesh : numpy array
         3-D array with the OAS format mesh
    """
    panelGeomY = np.broadcast_to(panelGeomY, (panelGeomX.shape[0], len(panelGeomY)))
    mesh = np.zeros((panelGeomX.shape[0], panelGeomY.shape[1], 3))
    mesh[:, :, 0] = panelGeomX[::-1]
    mesh[:, :, 1] = panelGeomY
    return mesh


if __name__ == "__main__":

    surface = {
        # Wing definition

        #Basic surface parameters
        "name":"surface",
        "num_sections": 2, #The number of sections in the multi-section surface
        "sec_name": ["sec0","sec1"],  # names of the individual sections
        "symmetry": True,  # if true, model one half of wing. reflected across the midspan of the root section
        "S_ref_type": "wetted",  # how we compute the wing area,
        # can be 'wetted' or 'projected'
        "rootSection": 1,


        #Geometry Parameters
        "sec_taper": np.array([1.0,1.0]), #Wing taper for each section
        "sec_span":np.array([1.0,1.0]), #Wing span for each section
        "sec_sweep":np.array([0.0,0.0]), #Wing sweep for each section
        "sec_chord_cp": [np.array([1,1]),np.array([1.0,0.2])],
        #"sec_chord_cp": [np.ones(1),2*np.ones(1),3*np.ones(1)], #Chord B-spline control points for each section
        "root_chord" : 1.0, #Wing root chord for each section

        #Mesh Parameters
        "meshes": "gen-meshes", #Supply a mesh for each section or "gen-meshes" for automatic mesh generation
        "nx" : 2, #Number of chordwise points. Same for all sections
        "sec_ny" : np.array([2,2]), #Number of spanwise points for each section
        
        #Aerodynamic Parameters
        "sec_CL0": np.array([0.0,0.0]),  # CL of the surface at alpha=0
        "sec_CD0": np.array([0.015,0.015]),  # CD of the surface at alpha=0
        # Airfoil properties for viscous drag calculation
        "k_lam": 0.05,  # percentage of chord with laminar
        # flow, used for viscous drag
        "sec_t_over_c_cp": [np.array([0.15]),np.array([0.15])],  # thickness over chord ratio (NACA0015)
        "sec_c_max_t": [0.303,0.303],  # chordwise location of maximum (NACA0015)
        # thickness
        "with_viscous": False,  # if true, compute viscous drag
        "with_wave": False,  # if true, compute wave drag
        "groundplane":False,
    }



    meshT, sec_meshes = generateMesh(surface)


    def plot_meshes(meshes):
        """ this function plots to plot the mesh """
        plt.figure(figsize=(8, 4))
        for i,mesh in enumerate(meshes):
            mesh_x = mesh[:, :, 0]
            mesh_y = mesh[:, :, 1]
            color = 'k'
            for i in range(mesh_x.shape[0]):
                plt.plot(mesh_y[i, :], mesh_x[i, :], color, lw=1)
                #plt.plot(-mesh_y[i, :], mesh_x[i, :], color, lw=1)   # plots the other side of symmetric wing
            for j in range(mesh_x.shape[1]):
                plt.plot(mesh_y[:, j], mesh_x[:, j], color, lw=1)
                #plt.plot(-mesh_y[:, j], mesh_x[:, j], color, lw=1)   # plots the other side of symmetric wing
        plt.axis('equal')
        plt.xlabel('y (m)')
        plt.ylabel('x (m)')
        #plt.legend()
  

    plot_meshes(sec_meshes)
    plt.show()

    
