"""module for contain other uncategorized methods"""

# Third party module
import numpy as np

# local module
import ligand
import integrin


def filter_by_dist(list_of_objects, max_dist, origin):
    """function to filter object base on maximum distance

    Parameter
    --------
    list_of_object:
        compilation list of the Ligands or Integrins
    max_dist:
        The maximum distance accepted
    origin:
        The position of the center/origin

    Return
    ------
    filtered list of objects

    """

    origin = np.array(origin)
    if list_of_objects:
        if isinstance(list_of_objects[0], ligand.Ligand):
            filtered1: list[ligand.Ligand] = []
            for member in list_of_objects:
                dist = np.linalg.norm(member.position - origin)
                if dist < max_dist:
                    filtered1.append(member)
            return filtered1
        if isinstance(list_of_objects[0], integrin.Integrin):
            filtered2: list[integrin.Integrin] = []
            for member in list_of_objects:
                dist = np.linalg.norm(member.position - origin)
                if dist < max_dist:
                    filtered2.append(member)
            return filtered2
    return []
