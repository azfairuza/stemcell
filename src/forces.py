import numpy as np
import src.physica as psc


def total_force_1(position_object, 
                  velocity_object, 
                  surface_integrins,
                  neighbor_integrins,
                  normal_length,
                  spring_constant,
                  damping_coefficient,
                  epsilon=1, 
                  sigma=1
                  ):
    """The total force for integrin in condition as follows:
    1. Surface integrin.
    2. The number of cell is more than 1.

    The forces
    ----------
    1. Lennard-Jones force from another surface integrin in another cell
    2. Spring force from neighboring integrins in the same cell 
    """
    forces_A = forceObjects_LJ(position_object, 
                                surface_integrins,
                                epsilon=epsilon,
                                sigma=sigma
                                )
    forces_B = forceObjects_Spring(position_object,
                                   velocity_object,
                                   neighbor_integrins,
                                   normal_length,
                                   spring_constant,
                                   damping_coefficient
                                   )
    total_force = forces_A + forces_B
    return total_force
            
def total_force_2(position_object, 
                  velocity_object, 
                  nearest_ligands,
                  neighbor_integrins,
                  normal_length,
                  spring_constant,
                  damping_coefficient,
                  epsilon=1, 
                  sigma=1
                  ):
    """The total force for integrin in condition as follows:
    1. Non surface integrin, or
    2. The number of cell is only one.

    The forces
    ----------
    1. Lennard-Jones force from nearest unoccupied ligands
    2. Spring force from neighboring integrins in the same cell 
    """
    forces_A = forceObjects_LJ(position_object, 
                                nearest_ligands,
                                epsilon=epsilon,
                                sigma=sigma
                                )
    forces_B = forceObjects_Spring(position_object,
                                   velocity_object,
                                   neighbor_integrins,
                                   normal_length,
                                   spring_constant,
                                   damping_coefficient
                                   )
    total_force = forces_A + forces_B
    return total_force

def forceObjects_LJ(position_object, 
                     collection_target, 
                     epsilon=1, 
                     sigma=1):
    """Function to calculate the total force from collection of objects
    using Lennard-Jones 6-12 potential.

    Parameter
    --------
    position_object: numpy.ndarray
        The position vector of the object.
    collection_target: physica.MultiObjBase
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
    total_Force = np.array([0, 0, 0])
    for target in collection_target:
        total_Force = total_Force + psc.force_LJ_6_12(target.position, 
                                                    position_object, 
                                                    epsilon, 
                                                    sigma)
    return total_Force

def forceObjects_Spring(position_object, 
                          velocity_object, 
                          collection_target, 
                          normal_length, 
                          spring_constant, 
                          damping_coefficient):
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
    total_Force = np.array([0, 0, 0])
    for target in collection_target:
        total_Force = total_Force + psc.force_spring(target.position,
                                                   position_object,
                                                   target.speed,
                                                   velocity_object,
                                                   spring_constant,
                                                   normal_length,
                                                   damping_coefficient)
    return total_Force
