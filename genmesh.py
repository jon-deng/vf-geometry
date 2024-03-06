"""
Generate some standard vocal fold (VF) meshes using GMSH
"""

# from typing import

from argparse import ArgumentParser

import gmsh

import numpy as np

def get_option_string(option: type):
    clscale = option.get_number('Mesh.MeshSizeFactor')
    return f'CL{clscale:.2e}'

def get_extrude_string(z_extrude, n_extrude):
    if z_extrude == 0:
        return f'DZ{z_extrude:.2f}'
    else:
        return f'DZ{z_extrude:.2e}--NZ{n_extrude:d}'
    

def gen_M5(medial_angle: float=0.0, z_extrude: float=0.0, n_extrude: int=1):
    """
    Generate a mesh for the M5_CB_GA*.STEP geometry
    """
    if z_extrude < 0:
        raise ValueError("`z_extrude` must be > 0")
    
    gmsh.clear()
    gmsh.model.add('main')

    gmsh.option.set_string('Geometry.OCCTargetUnit', 'CM')
    gmsh.merge(f'stp/M5_CB--GA{medial_angle:.0f}.STEP')

    medial_curve = [10]
    inf_curve = [12, 11]
    sup_curve = [9, 8]
    if z_extrude == 0:
        gmsh.model.add_physical_group(2, [2], name='body')
        gmsh.model.add_physical_group(2, [1], name='cover')

        gmsh.model.add_physical_group(
            1, medial_curve+inf_curve+sup_curve, 
            name='pressure'
        )
        gmsh.model.add_physical_group(1, [13, 7, 1], name='fixed')

        gmsh.model.add_physical_group(0, [10], name='separation-inf')
        gmsh.model.add_physical_group(0, [9], name='separation-sup')

    else:
        offset = len(gmsh.model.get_entities(dim=2))
        _medial_surf = [n+offset for n in medial_curve]
        _inf_surf = [n+offset for n in inf_curve]
        _sup_surf = [n+offset for n in sup_curve]

        _out_dim_tags = gmsh.model.occ.extrude(
            [(2, 1), (2, 2)],
            0.0, 0.0, z_extrude,
            [n_extrude]
        )
        gmsh.model.occ.synchronize()

        gmsh.model.add_physical_group(3, [2], name='body')
        gmsh.model.add_physical_group(3, [1], name='cover')

        gmsh.model.add_physical_group(
            2, [14, 13, 12, 11, 10],
            name='pressure'
        )

        surfs_anterior = [15, 17]
        surfs_posterior = [2, 1]
        surfs_lateral = [3, 16, 9]
        gmsh.model.add_physical_group(
            2, surfs_anterior+surfs_posterior+surfs_lateral,
            name='fixed'
        )

        gmsh.model.add_physical_group(1, [31], name='separation-inf')
        gmsh.model.add_physical_group(1, [29], name='separation-sup')
    
    # Add refinement at the medial surface
    gmsh.model.mesh.field.add('Distance', 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", medial_curve)
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    lc = 0.025
    gmsh.model.mesh.field.add('Threshold', 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc )
    gmsh.model.mesh.field.setNumber(2, "SizeMax", 5*lc)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 5*lc)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    if z_extrude == 0:
        gmsh.model.mesh.generate(2)
    else:
        gmsh.model.mesh.generate(3)

    option_string = get_option_string(gmsh.option)
    extrude_string = get_extrude_string(z_extrude, n_extrude)
    gmsh.write(f'M5_BC--GA{medial_angle:.2f}--{extrude_string}--{option_string}.msh')

def gen_LiEtal2020(medial_angle: float=0.0, z_extrude: float=0.0, n_extrude: int=1):
    """
    Generate a mesh for the 'LiEtal2020*.STEP' geometry
    """
    gmsh.clear()
    gmsh.model.add('main')

    gmsh.option.set_string('Geometry.OCCTargetUnit', 'CM')
    gmsh.merge(f'stp/LiEtal2020--GA{medial_angle:d}.STEP')

    if z_extrude == 0:
        gmsh.model.add_physical_group(2, [1], name='body')

        gmsh.model.add_physical_group(1, [1, 6, 5, 4, 3], name='pressure')
        gmsh.model.add_physical_group(1, [2], name='fixed')

        gmsh.model.add_physical_group(0, [5], name='separation-inf')
        gmsh.model.add_physical_group(0, [4], name='separation-sup')

        gmsh.model.mesh.generate(2)
    elif z_extrude > 0:
        out_dim_tags = gmsh.model.occ.extrude(
            [(2, 1)],
            0.0, 0.0, z_extrude,
            [n_extrude]
        )
        gmsh.model.occ.synchronize()

        gmsh.model.add_physical_group(3, [1], name='body')

        gmsh.model.add_physical_group(2, [2, 7, 6, 5, 4], name='pressure')
        gmsh.model.add_physical_group(2, [3], name='fixed')

        gmsh.model.add_physical_group(1, [14], name='separation-inf')
        gmsh.model.add_physical_group(1, [12], name='separation-sup')

        gmsh.model.mesh.generate()
    else:
        raise ValueError("`z_extrude` must be > 0")

    option_string = get_option_string(gmsh.option)
    gmsh.write(f'LiEtal2020--GA{medial_angle:.2f}--DZ{z_extrude:.2f}--{option_string}.msh')

def gen_M5_split(medial_angle: float=0.0, z_extrude: float=0, n_extrude: int=1):
    """
    Generate a mesh for the M5_CB_GA*_split.STEP geometry
    """
    gmsh.clear()
    gmsh.model.add('main')

    gmsh.option.set_string('Geometry.OCCTargetUnit', 'CM')
    gmsh.merge(f'stp/M5_CB_GA{medial_angle:d}_split.STEP')

    gmsh.model.add_physical_group(2, [3], name='body')
    gmsh.model.add_physical_group(2, [1, 2], name='cover')

    gmsh.model.add_physical_group(1, [6, 5, 4, 3, 15, 14], name='pressure')
    gmsh.model.add_physical_group(1, [13, 7, 16], name='fixed')

    gmsh.model.add_physical_group(0, [4], name='separation-inf')
    gmsh.model.add_physical_group(0, [3], name='separation-mid')
    gmsh.model.add_physical_group(0, [14], name='separation-sup')

    gmsh.model.mesh.generate(2)
    option_string = get_option_string(gmsh.option)
    gmsh.write(f'M5_CB_GA{medial_angle:.2f}_split--{option_string}.msh')

def gen_trapezoid(
        medial_angle: float=0.0,
        medial_surface_length: float=0.5,
        z_extrude: float=0,
        n_extrude: int=1
    ):
    """
    Generate a mesh for a trapezoidal geometry
    """
    gmsh.clear()
    gmsh.model.add('main')

    gmsh.option.set_string('Geometry.OCCTargetUnit', 'CM')

    # Origin
    gmsh.model.occ.add_point(0.0, 0.0, 0.0, tag=1)
    # Superior corner at base
    gmsh.model.occ.add_point(1.0, 0.0, 0.0, tag=2)
    # Superior point of medial surface (apex)
    coord_med_sup = np.array([1.0, 0.5, 0.0])
    gmsh.model.occ.add_point(*coord_med_sup, tag=3)
    # Inferior point of medial surface
    DEG = np.pi/180.0
    _dir = np.array([-1.0, -np.tan(medial_angle*DEG), 0.0])
    unit_dir = _dir/np.linalg.norm(_dir)
    coord_med_inf = coord_med_sup+medial_surface_length*unit_dir
    gmsh.model.occ.add_point(*coord_med_inf, tag=4)

    gmsh.model.occ.add_line(1, 2, tag=1)
    gmsh.model.occ.add_line(2, 3, tag=2)
    gmsh.model.occ.add_line(3, 4, tag=3)
    gmsh.model.occ.add_line(4, 1, tag=4)

    gmsh.model.occ.add_curve_loop([1, 2, 3, 4], tag=1)

    gmsh.model.occ.add_plane_surface([1], tag=1)

    gmsh.model.occ.synchronize()

    gmsh.model.add_physical_group(2, tags=[1], tag=1, name='VocalFold')

    gmsh.write("Trapezoid.geo_unrolled")
    gmsh.model.mesh.generate(2)

    option_string = get_option_string(gmsh.option)
    gmsh.write(f'Trapezoid{medial_angle:.2f}--{option_string}.msh')

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--geometry-name', type=str, default='M5')
    parser.add_argument('--medial-angle', type=float, default=0.0)
    parser.add_argument('--z-extrude', type=float, default=0.0)
    parser.add_argument('--n-extrude', type=int, default=1)
    parser.add_argument('--gmsh-args', type=str, default='')
    clargs = parser.parse_args()

    gmsh_args = ['gmsh'] + clargs.gmsh_args.split(' ')
    gmsh.initialize(gmsh_args)

    ## Parse the geometry name and any geometry/meshing parameters
    if clargs.geometry_name == 'M5':
        gen_mesh = gen_M5
    elif clargs.geometry_name == 'LiEtal2020':
        gen_mesh = gen_LiEtal2020
    elif clargs.geometry_name == 'M5Split':
        gen_mesh = gen_M5_split
    elif clargs.geometry_name == 'Trapezoid':
        gen_mesh = gen_trapezoid
    else:
        raise ValueError(f"Unknown 'geometry-name', {clargs.geometry_name}")

    mesh_params = {
        'medial_angle': clargs.medial_angle,
        'z_extrude': clargs.z_extrude,
        'n_extrude': clargs.n_extrude
    }

    gen_mesh(**mesh_params)
