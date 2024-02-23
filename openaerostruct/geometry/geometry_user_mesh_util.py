""" Utility for quickly generating a multi-section user specificed OAS mesh
"""

import numpy as np
import matplotlib.pyplot as plt


def generateMesh(
    sections, data, bPanels, cPanels, symmetry=True, controlPoints=False, generatePlots=False, plotOptions={}
):
    """
    Quickly generates user multi-section mesh for analysis based on a few parameters

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




    Returns
    -------
    mesh[nx, ny, 3] : numpy array
        Nodal mesh defining the aerodynamic surface.
    """

    # Data Initialization
    if np.shape(data) != (sections, 4):
        raise Exception("Data not specified for every section")

    if len(bPanels) != sections:
        raise Exception("Number of spanwise panels not specified for every section")

    panelGX, panelGY, panelQC, tquartX, Npanels = generateSectionGeometry(sections, data, bPanels, cPanels)
    panelGeomX, panelGeomY = stitchSectionGeometry(sections, panelGY, panelGX, bPanels)

    if symmetry:
        panelGeomX, panelGeomY = planformSymmetric(panelGeomX, panelGeomY, bPanels)
    if controlPoints:
        quarterC, tquarterPointsX = stitchPanelChordGeometry(sections, panelQC, tquartX, bPanels)
        controlPointsX, qControlPointsX, controlPointsY, chordDistGeom, chordDistCont = calcControlPoints(
            panelGeomY, panelGeomX, tquartX, panelQC, bPanels, cPanels
        )
    if generatePlots:
        if plotOptions["type"] == "plan":
            plotPlanform(sections, panelGX, panelGY, plotOptions["symmetry"], plotOptions["name"])
        elif plotOptions["type"] == "mesh":
            plotPanels(sections, panelGX, panelGY, plotOptions["symmetry"], plotOptions["name"])
        else:
            plotPlanform(sections, panelGX, panelGY, plotOptions["symmetry"], plotOptions["name"])
            plotPanels(sections, panelGX, panelGY, plotOptions["symmetry"], plotOptions["name"])

    mesh = outputOASMesh(panelGeomX, panelGeomY)

    return mesh


def generateSectionGeometry(sections, data, bPanels, cPanels):
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
    tquartX = []
    panelQC = []
    panelGY = []
    panelGX = []
    K = []

    tipStart = []
    tipEnd = []
    rootStart = []
    rootEnd = []

    Stot = 0
    span = 0
    Npanels = 2 * np.sum(bPanels) * cPanels

    for sec in range(sections):
        taper = data[sec, 0]
        AR = data[sec, 2]
        leLambda = np.deg2rad(data[sec, 3])

        if sec == 0:
            rootC = data[sec, 1]
        else:
            rootC = np.abs(panelGX[sec - 1][0, 0] - panelGX[sec - 1][cPanels, 0])

        tipC = rootC * taper
        S = AR * ((tipC + rootC) / 2) ** 2
        b = np.sqrt(AR * S)
        Stot = Stot + S
        span = span + b

        # Generate geomtery
        if sec == 0:
            rootEnd = 0
        else:
            rootEnd = panelGX[sec - 1][cPanels, 0]
        rootStart = rootC + rootEnd
        tipStart = rootStart - (b / 2) * np.tan(leLambda)
        tipEnd = tipStart - tipC

        secGeom = {}
        secGeom["rootStart"] = rootStart
        secGeom["rootEnd"] = rootEnd
        secGeom["tipStart"] = tipStart
        secGeom["tipEnd"] = tipEnd

        panelGY, panelGX, tquartX, panelQC, K = generatePanelsSection(
            sec, secGeom, cPanels, bPanels[sec], panelGX, panelGY, tquartX, panelQC, K, b
        )

    return panelGX, panelGY, panelQC, tquartX, Npanels


def generatePanelsSection(sec, secGeom, cPanels, bPanels, panelGX, panelGY, tquartX, panelQC, K, b):
    """
    Generates the mesh coordinates for a each wing section and appends them to the list of coordinates for each section.

    Parameters
    ----------
    sec : int
        The section number panels are being generated for
    secGeom : Dict
        Dictionary containing the rootStart, rootEnd, tipStart, and tipEnd x-coordinates
    bPanels: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section
    cPanels: int
        Integer for number of chord-wise panels
    panelGX : List
         List containing the mesh x-coordinates for each section
    panelGY : List
        List containing the mesh y-coordinates for each section
    panelQC : List
         List containing the quarter chord x-positions along each panel edge. Not used by OAS but information is too good to get rid of at this point.
    tquartX : : List
         List containing the three-quarter chord x-positions along each panel edge. Not used by OAS but information is too good to get rid of at this point.
    K : numpy array
        Array containing the panel widths for each section
    b : float
        Span length of the section(full span not half span)


    Returns
    -------
    panelGX : Dict
         Dictionary containing the mesh x-coordinates for each section
    panelGY : Dict
        Dictionary containing the mesh y-coordinates for each section
    panelQC : Dict
         Dictionary containing the quarter chord x-positions along each panel edge. Not used by OAS but information is too good to get rid of at this point.
    tquartX : : Dict
         Dictionary containing the three-quarter chord x-positions along each panel edge. Not used by OAS but information is too good to get rid of at this point.
    K : numpy array
        Array containing the panel widths for each section
    """
    # Generate Panels
    # Descretize wing root and wing tip into N+1 points

    rootStart = secGeom["rootStart"]
    rootEnd = secGeom["rootEnd"]
    tipStart = secGeom["tipStart"]
    tipEnd = secGeom["tipEnd"]

    # rootChordPoints = np.arange(rootStart,rootEnd-(rootStart-rootEnd)/cPanels,-(rootStart-rootEnd)/cPanels)
    rootChordPoints = np.linspace(rootStart, rootEnd, cPanels + 1)

    if tipStart == tipEnd:
        tipChordPoints = tipStart * np.ones(len(rootChordPoints))
    else:
        tipChordPoints = np.linspace(tipStart, tipEnd, len(rootChordPoints))

    K.append(b / (2 * bPanels))
    if sec == 0:
        panelGeomY = np.arange(-b / 2, b / 2 + K[sec], K[sec])
    else:
        panelGeomY = np.zeros(2 * bPanels)
        panelGeomY[0:bPanels] = np.linspace(panelGY[sec - 1][0] - (b / 2), panelGY[sec - 1][0] - K[sec], bPanels)
        panelGeomY[bPanels : 2 * bPanels] = np.linspace(
            panelGY[sec - 1][-1] + K[sec], panelGY[sec - 1][-1] + (b / 2), bPanels
        )

    if sec == 0:
        centreIndex = bPanels
    else:
        centreIndex = bPanels - 1

    panelGeomX = np.zeros([len(rootChordPoints), len(panelGeomY)])
    if sec == 0:
        for i in range(len(rootChordPoints)):
            # Left Wing
            panelGeomX[i, 0 : centreIndex + 1] = rootChordPoints[i] + (
                (tipChordPoints[i] - rootChordPoints[i]) / (-b / 2)
            ) * (panelGeomY[0 : centreIndex + 1])
            # Right Wing
            panelGeomX[i, centreIndex + 1 :] = rootChordPoints[i] + (
                (tipChordPoints[i] - rootChordPoints[i]) / (b / 2)
            ) * (panelGeomY[centreIndex + 1 :])
    else:
        for i in range(len(rootChordPoints)):
            # Left Wing
            panelGeomX[i, 0 : centreIndex + 1] = rootChordPoints[i] + (
                (tipChordPoints[i] - rootChordPoints[i]) / (-b / 2)
            ) * ((panelGeomY[0 : centreIndex + 1]) - panelGY[sec - 1][0])
            # Right Wing
            panelGeomX[i, centreIndex + 1 :] = rootChordPoints[i] + (
                (tipChordPoints[i] - rootChordPoints[i]) / (b / 2)
            ) * (panelGeomY[centreIndex + 1 :] - panelGY[sec - 1][-1])

    panelQuarterC = np.zeros([cPanels, len(panelGeomY)])
    tquarterPointsX = np.zeros([cPanels, len(panelGeomY)])

    for i in range(len(panelGeomY)):
        for j in range(cPanels - 1):
            panelQuarterC[j, i] = panelGeomX[j, i] + (panelGeomX[j + 1, i] - panelGeomX[j, i]) / 4
            tquarterPointsX[j, i] = panelGeomX[j, i] + 1 * (panelGeomX[j + 1, i] - panelGeomX[j, i]) / 2

    panelQC.append(panelQuarterC)
    panelGY.append(panelGeomY)
    tquartX.append(tquarterPointsX)
    panelGX.append(panelGeomX)

    return panelGY, panelGX, tquartX, panelQC, K


def stitchSectionGeometry(sections, panelGY, panelGX, bPanels):
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
    for i in range(sections):
        if i == 0:
            panelGeomY = panelGY[i]
            panelGeomX = panelGX[i]
        else:
            panelGeomY = np.concatenate((panelGY[i][0 : bPanels[i]], panelGeomY, panelGY[i][bPanels[i] :]))
            panelGeomX = np.concatenate(
                (panelGX[i][:, 0 : bPanels[i]], panelGeomX, panelGX[i][:, bPanels[i] :]), axis=1
            )
    return panelGeomX, panelGeomY


def stitchPanelChordGeometry(sections, panelQC, tquartX, bPanels):
    """
    Combines the split section arrays for quarter chord and three-quarters chord x-coordintes into singular unified mesh

    Parameters
    ----------
    sections : int
        Integer for number of wing sections specified
    panelQC : List
         List containing the mesh quarter chord x-coordinates for each section
    tquartX : List
        List containing the mesh three-quarter chord x-coordinates for each section
    bPanels: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section


    Returns
    -------
    quarterC : numpy array
         Array of the quarter chord x-coordinates
    tquarterPointsX : numpy array
        Array of the three-quarter chord x-coordinates
    """
    for i in range(sections):
        if i == 0:
            quarterC = panelQC[i]
            tquarterPointsX = tquartX[i]
        else:
            quarterC = np.concatenate((panelQC[i][:, 0 : bPanels[i]], quarterC, panelQC[i][:, bPanels[i] :]), axis=1)
            tquarterPointsX = np.concatenate(
                (tquartX[i][:, 0 : bPanels[i]], tquarterPointsX, tquartX[i][:, bPanels[i] :]), axis=1
            )
    return quarterC, tquarterPointsX


def calcControlPoints(panelGeomY, panelGeomX, tquartX, panelQC, bPanels, cPanels):
    """
    Calculates the control point locations on each panel

    Parameters
    ----------
    panelGeomX : numpy array
         Array of the mesh x-coordinates
    panelGeomY : numpy array
        Array of the mesh y-coordinates
    bPanels: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section
    cPanels: int
        Integer for number of chord-wise panels
    panelQC : numpy array
         Array of the quarter chord x-coordinates
    tquartX : : numpy array
         Array of the three-quarter chord x-coordinates


    Returns
    -------
    controlPointsX : numpy array
         Array of the three-quarter chord control point x-coordinates
    qControlPointsX : numpy array
        Array of the quarter chord control point x-coordinates
    controlPointsY : numpy array
        Array of the control point y-coordinates
    chordDistGeom : numpy array
        The geometric chord distribution
    chordDistCont : numpy array
        Array of chord distribution at the control points
    """
    controlPointsX = np.zeros([cPanels, 2 * np.sum(bPanels)])
    qControlPointsX = np.zeros_like(controlPointsX)
    controlPointsY = np.zeros(2 * np.sum(bPanels))
    chordDistGeom = np.zeros_like(panelGeomY)
    chordDistCont = np.zeros_like(controlPointsY)

    for i in range(len(panelGeomY)):
        chordDistGeom[i] = panelGeomX[0, i] - panelGeomX[-1, i]

    for i in range(len(panelGeomY) - 1):
        for j in range(cPanels):
            controlPointsX[j, i] = (tquartX[j, i + 1] + tquartX[j, i]) / 2
            qControlPointsX[j, i] = (panelQC[j, i + 1] + panelQC[j, i]) / 2
        controlPointsY[i] = panelGeomY[i] + (panelGeomY[i + 1] - panelGeomY[i]) / 2
        chordDistCont[i] = (chordDistGeom[i + 1] + chordDistGeom[i]) / 2

    return controlPointsX, qControlPointsX, controlPointsY, chordDistGeom, chordDistCont


def planformSymmetric(panelGeomX, panelGeomY, bPanels):
    """
    Cuts the mesh in half is symmetry conditions are used

    Parameters
    ----------
    panelGeomX : numpy array
         Array of the mesh x-coordinates
    panelGeomY : numpy array
        Array of the mesh y-coordinates
    bPanels: numpy array
        1 x sections 1-D array consisting of integers that specify the number of spanwise panels corresponding to each wing section


    Returns
    -------
    panelGeomX : numpy array
         Array of the mesh x-coordinates
    panelGeomY : numpy array
        Array of the mesh y-coordinates
    """
    panelGeomX = panelGeomX[:, 0 : np.sum(bPanels)]
    panelGeomY = panelGeomY[0 : np.sum(bPanels)]
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
    mesh[:, :, 0] = panelGeomX
    mesh[:, :, 1] = panelGeomY
    return mesh


def plotPlanform(sections, panelGX, panelGY, plotSymmetry="Left", wingName="CustomUserWing"):
    """
    Plots the outline of all definied planform sections.

    Parameters
    ----------
    sections : int
        Integer for number of wing section specified
    panelGX: list[numpy array]
        Panel X-coordinate data
    panelGY: list[numpy array]
        Panel Y-coordinate data
    plotSymmetry : string
        Flag set to 'Left' if only y<=0 plotted
        Flag set to 'Right' if only y>=0 plotted
        Flag set to 'None' if full planform plotted
    wingName: string
        Name for plot title

    Returns
    -------
    Figure
    """
    np.random.seed(91)
    colorSet = np.random.rand(sections, 3)
    plt.figure()
    for i in range(sections):
        if i == 0:

            rootEnd = panelGX[i][cPanels, bPanels[i]]
            rootStart = panelGX[i][0, bPanels[i]]
            # tipStart = panelGX[i][0,0]
            # tipEnd = panelGX[i][cPanels,0]

            if plotSymmetry == "Left":
                plt.plot(
                    [0, 0, panelGY[i][0], panelGY[i][0], 0],
                    [rootStart, rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart],
                    c=colorSet[i, :],
                )
            elif plotSymmetry == "Right":
                plt.plot(
                    [0, panelGY[i][-1], panelGY[i][-1], 0, 0],
                    [rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart, rootStart],
                    c=colorSet[i, :],
                )
            else:
                plt.plot(
                    [0, panelGY[i][0], panelGY[i][0], 0],
                    [rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart],
                    c=colorSet[i, :],
                )
                plt.plot(
                    [0, panelGY[i][-1], panelGY[i][-1], 0],
                    [rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart],
                    c=colorSet[i, :],
                )

        else:

            rootEnd = panelGX[i - 1][cPanels, 0]
            rootStart = panelGX[i - 1][0, 0]
            # tipStart = panelGX[i][0,0]
            # tipEnd = panelGX[i][cPanels,0]

            if plotSymmetry == "Left":
                plt.plot(
                    [panelGY[i - 1][0], panelGY[i - 1][0], panelGY[i][0], panelGY[i][0], panelGY[i - 1][0]],
                    [rootStart, rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart],
                    c=colorSet[i, :],
                )
            elif plotSymmetry == "Right":
                plt.plot(
                    [panelGY[i - 1][-1], panelGY[i][-1], panelGY[i][-1], panelGY[i - 1][-1], panelGY[i - 1][-1]],
                    [rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart, rootStart],
                    c=colorSet[i, :],
                )
            else:
                plt.plot(
                    [panelGY[i - 1][0], panelGY[i - 1][0], panelGY[i][0], panelGY[i][0], panelGY[i - 1][0]],
                    [rootStart, rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart],
                    c=colorSet[i, :],
                )
                plt.plot(
                    [panelGY[i - 1][-1], panelGY[i][-1], panelGY[i][-1], panelGY[i - 1][-1], panelGY[i - 1][-1]],
                    [rootEnd, panelGX[i][cPanels, 0], panelGX[i][0, 0], rootStart, rootStart],
                    c=colorSet[i, :],
                )
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("{} Planform Geometry".format(wingName))
    plt.savefig("{}PlanformGeom.pdf".format(wingName))


def plotPanels(sections, panelGX, panelGY, plotSymmetry="Left", wingName="CustomUserWing"):
    """
    Plots the mesh(panels) of each planform section

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
        Flag set to 'None' if full planform plotted
    wingName: string
        Name for plot title

    Returns
    -------
    Figure
    """
    np.random.seed(151)
    colorSet = np.random.rand(sections, 3)
    plt.figure()
    for i in range(sections):
        if i == 0:
            for j in range(len(panelGX[i][:, bPanels[i]])):
                if plotSymmetry == "Left":
                    plt.plot([0, panelGY[i][0]], [panelGX[i][j, bPanels[i]], panelGX[i][j, 0]], c=colorSet[i, :])
                    for k in range(int(np.ceil(len(panelGY[i]) / 2))):
                        plt.plot(
                            [panelGY[i][k], panelGY[i][k]], [panelGX[i][0, k], panelGX[i][-1, k]], c=colorSet[i, :]
                        )
                elif plotSymmetry == "Right":
                    plt.plot([0, panelGY[i][-1]], [panelGX[i][j, bPanels[i]], panelGX[i][j, 0]], c=colorSet[i, :])
                    for k in range(int(np.floor(len(panelGY[i]) / 2)), len(panelGY[i])):
                        plt.plot(
                            [panelGY[i][k], panelGY[i][k]], [panelGX[i][0, k], panelGX[i][-1, k]], c=colorSet[i, :]
                        )
                else:
                    plt.plot([0, panelGY[i][0]], [panelGX[i][j, bPanels[i]], panelGX[i][j, 0]], c=colorSet[i, :])
                    plt.plot([0, panelGY[i][-1]], [panelGX[i][j, bPanels[i]], panelGX[i][j, 0]], c=colorSet[i, :])
                    for k in range(len(panelGY[i])):
                        plt.plot(
                            [panelGY[i][k], panelGY[i][k]], [panelGX[i][0, k], panelGX[i][-1, k]], c=colorSet[i, :]
                        )

        else:
            for j in range(len(panelGX[i - 1][:, 0])):
                if plotSymmetry == "Left":
                    plt.plot(
                        [panelGY[i - 1][0], panelGY[i][0]], [panelGX[i - 1][j, 0], panelGX[i][j, 0]], c=colorSet[i, :]
                    )
                    plotRange = range(int(np.ceil(len(panelGY[i]) / 2)))
                    for k in plotRange:
                        plt.plot(
                            [panelGY[i][k], panelGY[i][k]], [panelGX[i][0, k], panelGX[i][-1, k]], c=colorSet[i, :]
                        )
                elif plotSymmetry == "Right":
                    plt.plot(
                        [panelGY[i - 1][-1], panelGY[i][-1]], [panelGX[i - 1][j, 0], panelGX[i][j, 0]], c=colorSet[i, :]
                    )
                    for k in range(int(np.floor(len(panelGY[i]) / 2)), len(panelGY[i])):
                        plt.plot(
                            [panelGY[i][k], panelGY[i][k]], [panelGX[i][0, k], panelGX[i][-1, k]], c=colorSet[i, :]
                        )
                else:
                    plt.plot(
                        [panelGY[i - 1][0], panelGY[i][0]], [panelGX[i - 1][j, 0], panelGX[i][j, 0]], c=colorSet[i, :]
                    )
                    plt.plot(
                        [panelGY[i - 1][-1], panelGY[i][-1]], [panelGX[i - 1][j, 0], panelGX[i][j, 0]], c=colorSet[i, :]
                    )
                    for k in range(len(panelGY[i])):
                        plt.plot(
                            [panelGY[i][k], panelGY[i][k]], [panelGX[i][0, k], panelGX[i][-1, k]], c=colorSet[i, :]
                        )
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("{} Mesh".format(wingName))
    plt.savefig("{}Mesh.pdf".format(wingName))


if __name__ == "__main__":
    """Generate example meshes"""

    plotOptions = {}
    plotOptions["type"] = "all"
    plotOptions["symmetry"] = "Left"
    plotOptions["name"] = "Rectangular Wing"

    # Simple rectangular wing
    # 32 x 16 mesh
    sections = 1
    data = np.array([[1, 1, 10, 0]])
    bPanels = np.array([16])
    cPanels = 16
    mesh = generateMesh(sections, data, bPanels, cPanels, True, False, True, plotOptions)

    # Double Delta Wing 75/65
    plotOptions["name"] = "Double Delta Wing 75\\65"

    sections = 2
    bPanels = np.array([6, 7])
    cPanels = 16
    data = np.array([[0.39319, 14.4533, 0.466822, 75], [0, 5.6829, 1.86524, 65]])
    mesh = generateMesh(sections, data, bPanels, cPanels, True, False, True, plotOptions)

    # NASA N+3 SST
    plotOptions["name"] = "NASA N+3 SST"
    sections = 4
    data = np.array(
        [
            [0.648294, 187.23, 0.075305, 84.9578],
            [0.698608, 121.37, 0.089833, 82.3695],
            [0.345206, 84.79, 0.573734, 68],
            [0.205671, 29.27, 3.46289, 40],
        ]
    )
    bPanels = np.array([3, 2, 8, 15])
    cPanels = 16

    mesh = generateMesh(sections, data, bPanels, cPanels, True, False, True, plotOptions)
    plt.show()
