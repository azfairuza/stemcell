import src.physica as psc


def readFile(filename: str):
    """ Procedure to read input file. 
    
    There are three input file that can be read.
    PATCON: contains the configuration for building nanopattern.
    CELCON: contains the configuration for building cells.
    SIMCON: contains the configuration for running the simulation.

    parameter
    ---------
    filename: str
        the name of the file, may also the full path of the file.

    return
    ------
    stripped_srting: list of str
    """
    if filename == "PATCON":
        try:
            data_file = open("input\PATCON.txt", "r", encoding="utf-8")
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = psc.filterItem(
                lst_strng
            ) 
        except:
            stripped_strng = (
                "Error in opening PATCON file"
            )
        finally:
            return stripped_strng

    elif filename == "CELCON":
        try:
            data_file = open("input\CELCON.txt", "r", encoding="UTF-8")
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = psc.filterItem(
                lst_strng
            )
        except:
            stripped_strng = (
                "Error in opening CELCON file" 
            )
        finally:
            return stripped_strng

    elif filename == "SIMCON":
        try:
            data_file = open("input\SIMCON.txt", "r", encoding="utf-8")
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = psc.filterItem(
                lst_strng
            )
        except:
            stripped_strng = (
                "Error in opening SIMCON file"
            )
        finally:
            return stripped_strng

    else:
        print("ERROR: File name is not correct")
        return 0

def getValue(lst_strng: list, property_name: str):
    """Procedure to get the value of a property in XXXCON file.
    There are special property name to be used in this function.
    PATCON: property name = 'LIGAND'.
    CELCON: property name = 'CELL'.
    SIMCON: property name = 'METADATA'.

    parameters
    ----------
    lst_strng: list
        the list of string which contain the property and value
    property_name: str
        the name of property or attribute

    return
    ------
    value: float or else
        depend on the value of the property
    """
    # 'Ligand'
    # if property_name == 'LIGAND':
    #     start_index = (lst_strng.index('#LIGAND') + 2)
    #     end_index = lst_strng.index('#END')
    #     ligand_position = []
    #     for i in range(start_index, end_index):
    #         dot_position = [float(j) for j in lst_strng[i].split()]
    #         ligand_position.append(dot_position)
    #     return ligand_position
    if property_name == "xdist" or property_name == "ydist":
        start_index = lst_strng.index("#CONFIG") + 1
        end_index = lst_strng.index("#END")
        value = None
        for i in range(start_index, end_index):
            data = lst_strng[i].split()
            raw_value = data[1:]
            value = [float(i) for i in raw_value]
        return value
    elif property_name == "CELL":
        start_index = lst_strng.index("#CELL") + 2
        end_index = lst_strng.index("#END")
        cell_properties = []
        for i in range(start_index, end_index):
            cell_property = [float(j) for j in lst_strng[i].split()]
            cell_properties.append(cell_property)
        return cell_properties
    elif property_name == "METADATA":
        start_index = lst_strng.index("#METADATA") + 1
        end_index = lst_strng.index("#CONFIG")
        cell_properties = {"username": None, "title": None}
        for i in range(start_index, end_index):
            stringdata = lst_strng[i].split()
            if stringdata[0] == "username":
                cell_properties["username"] = " ".join(stringdata[1:])
            elif stringdata[0] == "title":
                cell_properties["title"] = " ".join(stringdata[1:])
        return cell_properties
    else:
        start_index = lst_strng.index("#CONFIG") + 1
        end_index = lst_strng.index("#END")
        value = None
        for i in range(start_index, end_index):
            data = lst_strng[i].split()
            if data[0] == property_name:
                if len(data) > 2:
                    raw_value = data[1:]
                    value = [float(i) for i in raw_value]
                elif len(data) == 2:
                    value = float(data[1])
                break
        return value
