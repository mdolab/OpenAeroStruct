# Main file for coupled system components

from __future__ import print_function, division
import warnings
import numpy
import sys
import os
from openmdao.api import IndepVarComp, Problem, Group
from geometry import GeometryMesh, Bspline, gen_crm_mesh, gen_mesh, get_inds
from transfer import TransferDisplacements, TransferLoads
from materials import MaterialsTube
from b_spline import get_bspline_mtx
from spatialbeam import radii, SpatialBeamStates
from vlm import VLMStates

warnings.filterwarnings("ignore") # to disable OpenMDAO warnings which will create an error in Matlab

# In Matlab code, add this before calling Python functions:
#   if count(py.sys.path,'') == 0
#       insert(py.sys.path,int32(0),'');
#   end

def setup(num_inboard=2, num_outboard=3, check=False, out_stream=sys.stdout):
    ''' Setup the aerostruct mesh using OpenMDAO'''

    # Define the aircraft properties
    from CRM import span, v, alpha, rho

    # Define spatialbeam properties
    from aluminum import E, G, stress, mrho

    # Create the mesh with 2 inboard points and 3 outboard points.
    # This will be mirrored to produce a mesh with 7 spanwise points,
    # or 6 spanwise panels
    # print(type(num_inboard))
    mesh = gen_crm_mesh(int(num_inboard), int(num_outboard), num_x=2)
    num_x, num_y = mesh.shape[:2]
    num_twist = numpy.max([int((num_y - 1) / 5), 5])
    r = radii(mesh)

    # Set the number of thickness control points and the initial thicknesses
    num_thickness = num_twist
    t = r / 10

    mesh = mesh.reshape(-1, mesh.shape[-1])
    aero_ind = numpy.atleast_2d(numpy.array([num_x, num_y]))
    fem_ind = [num_y]
    aero_ind, fem_ind = get_inds(aero_ind, fem_ind)

    # Set additional mesh parameters
    dihedral = 0.  # dihedral angle in degrees
    sweep = 0.  # shearing sweep angle in degrees
    taper = 1.  # taper ratio

    # Initial displacements of zero
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    disp = numpy.zeros((tot_n_fem, 6))

    # Define Jacobians for b-spline controls
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    num_surf = fem_ind.shape[0]
    jac_twist = get_bspline_mtx(num_twist, num_y)
    jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem-num_surf)

    # Define ...
    twist_cp = numpy.zeros(num_twist)
    thickness_cp = numpy.ones(num_thickness)*numpy.max(t)

    # Define the design variables
    des_vars = [
        ('twist_cp', twist_cp),
        ('dihedral', dihedral),
        ('sweep', sweep),
        ('span', span),
        ('taper', taper),
        ('v', v),
        ('alpha', alpha),
        ('rho', rho),
        ('disp', disp),
        ('aero_ind', aero_ind),
        ('fem_ind', fem_ind)
    ]

    root = Group()

    root.add('des_vars',
         IndepVarComp(des_vars),
         promotes=['twist_cp','span','v','alpha','rho','disp','dihedral'])
    root.add('mesh',  # This component is needed, otherwise resulting loads matrix is NaN
         GeometryMesh(mesh, aero_ind), # changes mesh given span, sweep, twist, and des_vars
         promotes=['span','sweep','dihedral','twist','taper','mesh'])
    root.add('def_mesh',
         TransferDisplacements(aero_ind, fem_ind),
         promotes=['mesh','disp','def_mesh'])

    prob = Problem()
    prob.root = root
    prob.setup(check=check, out_stream=out_stream)
    prob.run()

    # Output the def_mesh for the aero modules
    def_mesh = prob['def_mesh']

    # Other variables needed for aero and struct modules
    params = {
        'mesh': mesh,
        'num_x': num_x,
        'num_y': num_y,
        'span': span,
        'twist_cp': twist_cp,
        'thickness_cp': thickness_cp,
        'v': v,
        'alpha': alpha,
        'rho': rho,
        'r': r,
        't': t,
        'aero_ind': aero_ind,
        'fem_ind': fem_ind,
        'num_thickness': num_thickness,
        'num_twist': num_twist,
        'sweep': sweep,
        'taper': taper,
        'dihedral': dihedral,
        'E': E,
        'G': G,
        'stress': stress,
        'mrho': mrho,
        'tot_n_fem': tot_n_fem,
        'num_surf': num_surf,
        'jac_twist': jac_twist,
        'jac_thickness': jac_thickness,
        'out_stream': out_stream,
        'check': check
    }

    return (def_mesh, params)

def aero(def_mesh=None, params=None):

    if (params is None) or (def_mesh is None):
        tmpmesh, params = setup()
        if not def_mesh:
            def_mesh = tmpmesh

    # Unpack variables
    mesh = params.get('mesh')
    num_x = params.get('num_x')
    num_y = params.get('num_y')
    span = params.get('span')
    twist_cp = params.get('twist_cp')
    thickness_cp = params.get('thickness_cp')
    v = params.get('v')
    alpha = params.get('alpha')
    rho = params.get('rho')
    r = params.get('r')
    t = params.get('t')
    aero_ind = params.get('aero_ind')
    fem_ind = params.get('fem_ind')
    num_thickness = params.get('num_thickness')
    num_twist = params.get('num_twist')
    sweep = params.get('sweep')
    taper = params.get('taper')
    disp = params.get('disp')
    dihedral = params.get('dihedral')
    check = params.get('check')
    out_stream = params.get('out_stream')

    # Define Jacobians for b-spline controls
    tot_n_fem = numpy.sum(fem_ind[:, 0])
    num_surf = fem_ind.shape[0]
    jac_twist = get_bspline_mtx(num_twist, num_y)
    jac_thickness = get_bspline_mtx(num_thickness, tot_n_fem-num_surf)

    disp = numpy.zeros((num_y, 6))  # for display?

    # Define the design variables
    des_vars = [
        ('twist_cp', numpy.zeros(num_twist)),
        ('dihedral', dihedral),
        ('sweep', sweep),
        ('span', span),
        ('taper', taper),
        ('v', v),
        ('alpha', alpha),
        ('rho', rho),
        ('disp', numpy.zeros((tot_n_fem, 6))),
        ('def_mesh', def_mesh)
    ]

    root = Group()

    root.add('des_vars',
             IndepVarComp(des_vars),
             promotes=['twist_cp','span','v','alpha','rho','disp','dihedral','def_mesh'])
    # root.add('twist_bsp',  # What is this doing?
    #          Bspline('twist_cp', 'twist', jac_twist),
    #          promotes=['*'])
    # root.add('thickness_bsp',    # What is this doing?
    #          Bspline('thickness_cp', 'thickness', jac_thickness),
    #          promotes=['*'])
    # root.add('tube',
    #          MaterialsTube(fem_ind),
    #          promotes=['*'])

    root.add('VLMstates',
             VLMStates(aero_ind),
             promotes=[
                'def_mesh','b_pts','mid_b','c_pts','widths','normals','S_ref', # VLMGeometry
                'alpha','circulations','v',  # VLMCirculations
                'rho','sec_forces'           # VLMForces
             ])
    root.add('loads',
             TransferLoads(aero_ind, fem_ind),
             promotes=['def_mesh','sec_forces','loads'])

    prob = Problem()
    prob.root = root
    prob.setup(check=check, out_stream=out_stream)
    prob.run()

    return prob['loads']  # Output the Loads matrix


def struct(loads, params):

    # Unpack variables
    mesh = params.get('mesh')
    num_x = params.get('num_x')
    num_y = params.get('num_y')
    span = params.get('span')
    twist_cp = params.get('twist_cp')
    thickness_cp = params.get('thickness_cp')
    v = params.get('v')
    alpha = params.get('alpha')
    rho = params.get('rho')
    r = params.get('r')
    t = params.get('t')
    aero_ind = params.get('aero_ind')
    fem_ind = params.get('fem_ind')
    num_thickness = params.get('num_thickness')
    num_twist = params.get('num_twist')
    sweep = params.get('sweep')
    taper = params.get('taper')
    disp = params.get('disp')
    dihedral = params.get('dihedral')
    E = params.get('E')
    G = params.get('G')
    stress = params.get('stress')
    mrho = params.get('mrho')
    tot_n_fem = params.get('tot_n_fem')
    num_surf = params.get('num_surf')
    jac_twist = params.get('jac_twist')
    jac_thickness = params.get('jac_thickness')
    check = params.get('check')
    out_stream = params.get('out_stream')

    # Define the design variables
    des_vars = [
        ('twist_cp', twist_cp),
        ('thickness_cp', thickness_cp),
        ('r', r),
        # ('dihedral', dihedral),
        # ('sweep', sweep),
        # ('span', span),
        # ('taper', taper),
        # ('v', v),
        # ('alpha', alpha),
        # ('rho', rho),
        # ('disp', disp),
        ('loads', loads)
    ]

    root = Group()

    root.add('des_vars',
             IndepVarComp(des_vars),
            #  promotes=['thickness_cp','r','loads'])
            promotes=['*'])
    # root.add('twist_bsp',  # What is this doing?
    #          Bspline('twist_cp', 'twist', jac_twist),
    #          promotes=['*'])
    root.add('twist_bsp',
             Bspline('twist_cp', 'twist', jac_twist),
             promotes=['*'])
    root.add('thickness_bsp',    # What is this doing?
             Bspline('thickness_cp', 'thickness', jac_thickness),
            #  promotes=['thickness'])
            promotes=['*'])
    root.add('mesh',
             GeometryMesh(mesh, aero_ind),
             promotes=['*'])
    root.add('tube',
             MaterialsTube(fem_ind),
            #  promotes=['r','thickness','A','Iy','Iz','J'])
            promotes=['*'])
    root.add('spatialbeamstates',
             SpatialBeamStates(aero_ind, fem_ind, E, G),
            #  promotes=[
            #     'mesh', # ComputeNodes
            #     'A','Iy','Iz','J','loads', # SpatialBeamFEM
            #     'disp' # SpatialBeamDisp
            #  ])
            promotes=['*'])
    root.add('transferdisp',
             TransferDisplacements(aero_ind, fem_ind),
            #  promotes=['mesh','disp','def_mesh'])
            promotes=['*'])

    prob = Problem()
    prob.root = root
    prob.setup(check=check, out_stream=out_stream)
    prob.run()

    return prob['def_mesh']  # Output the def_mesh matrix
