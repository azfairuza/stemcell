"""Module that contain the collection of forces for each scenario of
integrin.
"""
# built-in import
from __future__ import annotations

# third party imports
import numpy as np

# local imports
import physica as psc


def total_force_1(
    position_object,
    velocity_object,
    surface_integrins,
    neighbor_integrins,
    normal_length,
    spring_constant,
    damping_coefficient,
    viscocity,
    epsilon=1,
    sigma=1,
    dim=3
):
    """The total force for integrin in condition as follows:
    1. Surface integrin.
    2. The number of cell is more than 1.

    The forces
    ----------
    1. Lennard-Jones force from another surface integrin in another cell
    2. Spring force from neighboring integrins in the same cell
    """
    forces_a = force_objects_lennardjones(
        position_object, surface_integrins, epsilon=epsilon, sigma=sigma, dim=dim
    )
    forces_b = force_objects_spring(
        position_object,
        velocity_object,
        neighbor_integrins,
        normal_length,
        spring_constant,
        damping_coefficient,
        dim=dim
    )
    forces_c = psc.force.drag(velocity_object, sigma, viscocity)
    total_force = forces_a + forces_b + forces_c
    return total_force


def total_force_2(
    position_object,
    velocity_object,
    nearest_ligands,
    neighbor_integrins,
    normal_length,
    spring_constant,
    damping_coefficient,
    viscocity,
    epsilon=1,
    sigma=1,
    dim=3
):
    """The total force for integrin in condition as follows:
    1. Non surface integrin, or
    2. The number of cell is only one.

    The forces
    ----------
    1. Lennard-Jones force from nearest unoccupied ligands
    2. Spring force from neighboring integrins in the same cell
    """
    forces_a = force_objects_lennardjones(
        position_object, nearest_ligands, epsilon=epsilon, sigma=sigma, dim=dim
    )
    forces_b = force_objects_spring(
        position_object,
        velocity_object,
        neighbor_integrins,
        normal_length,
        spring_constant,
        damping_coefficient,
        dim=dim
    )
    forces_c = psc.force.drag(velocity_object, sigma, viscocity)
    total_force = forces_a + forces_b + forces_c
    return total_force


def force_objects_lennardjones(
    position_object, collection_target: list[psc.ObjBase], epsilon=1, sigma=1, dim=3
):
    """Function to calculate the total force from collection of objects
    using Lennard-Jones 6-12 potential.

    Parameter
    --------
    position_object: numpy.ndarray
        The position vector of the object.
    collection_target: list
        The list of object's targets.
    epsilon: float
        The depth of LJ potential.
    sigma: float
        The position where the potential is equal to Zero.

    Return
    ------
    totalForce:
        The total force from the objects.

    """
    total_force = np.zeros(dim)
    for target in collection_target:
        total_force = total_force + psc.force.lj_6_12(
            target.position, position_object, epsilon, sigma
        )
    return total_force


def force_objects_spring(
    position_object,
    velocity_object,
    collection_target: list[psc.ObjBase],
    normal_length,
    spring_constant,
    damping_coefficient,
    dim=3
):
    """Function to calculate the total force from collection of objects
    using spring potential.

    Parameter
    --------
    position_object: numpy.ndarray
        The position vector of the object.
    velocity_object: numpy.ndarray
        The velocity vector of the object.
    collection_target: :obj: list[Integrin] | list[Ligand]
        The list of positions of object's target.
    normal_length: float
        The normal_length of the spring.
    spring_constant: float
        The spring constant value
    damping_coefficient: float
        The damping coefficient value

    Return
    ------
    totalForce:
        The total force from the objects

    """
    total_force = np.zeros(dim)
    for target in collection_target:
        total_force = total_force + psc.force.spring(
            target.position,
            position_object,
            target.velocity,
            velocity_object,
            spring_constant,
            normal_length,
            damping_coefficient,
        )
    return total_force
