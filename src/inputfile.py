"""Module related to read the input files"""

class Read:
    def __init__(self, filename: str, stripped_strng = None) -> None:
        """Procedure to read input file.

        There are three input file that can be read.
        PATCON: contains the configuration for building nanopattern.
        CELCON: contains the configuration for building cells.
        SIMCON: contains the configuration for running the simulation.

        parameter
        ---------
        filename: str
            the name of the file, may also the full path of the file.
        """
        if filename == "PATCON":
            try:
                with open("./input/PATCON.txt", "r", encoding="utf-8") as data_file:
                    lst_strng = data_file.readlines()
                stripped_strng = filter_item(lst_strng)
            except:
                print("Error in opening PATCON file")
                stripped_strng = "Error in opening PATCON file"

        elif filename == "CELCON":
            try:
                with open("./input/CELCON.txt", "r", encoding="utf-8") as data_file:
                    lst_strng = data_file.readlines()
                stripped_strng = filter_item(lst_strng)
            except:
                print("Error in opening CELCON file")
                stripped_strng = "Error in opening CELCON file"

        elif filename == "SIMCON":
            try:
                with open("./input/SIMCON.txt", "r", encoding="utf-8") as data_file:
                    lst_strng = data_file.readlines()
                stripped_strng = filter_item(lst_strng)
            except:
                print("Error in opening SIMCON file")
                stripped_strng = "Error in opening SIMCON file"
        elif filename == "OUTPUT":
            pass
        
        else:
            print("ERROR: File name is not correct")
            stripped_strng = "ERROR: File name is not correct"
        self.contents = stripped_strng
    
    def get(self, property_name: str):
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

        # xdist and ydits
        if property_name in ("xdist", "ydist"):
            start_index = self.contents.index("#CONFIG") + 1
            end_index = self.contents.index("#END")
            value = None
            for i in range(start_index, end_index):
                data = self.contents[i].split()
                raw_value = data[1:]
                value = [float(i) for i in raw_value]
            return value

        # Cell
        if property_name == "CELL":
            start_index = self.contents.index("#CELL") + 2
            end_index = self.contents.index("#END")
            cell_properties = []
            for i in range(start_index, end_index):
                cell_property = [float(j) for j in self.contents[i].split()]
                cell_properties.append(cell_property)
            return cell_properties

        # Metadata
        if property_name == "METADATA":
            start_index = self.contents.index("#METADATA") + 1
            end_index = self.contents.index("#CONFIG")
            cell_properties = {"username": None, "title": None}
            for i in range(start_index, end_index):
                stringdata = self.contents[i].split()
                if stringdata[0] == "username":
                    cell_properties["username"] = " ".join(stringdata[1:])
                elif stringdata[0] == "title":
                    cell_properties["title"] = " ".join(stringdata[1:])
            return cell_properties

        # default mode
        start_index = self.contents.index("#CONFIG") + 1
        end_index = self.contents.index("#END")
        value = None
        for i in range(start_index, end_index):
            data = self.contents[i].split()
            if data[0] == property_name:
                if len(data) > 2:
                    raw_value = data[1:]
                    value = [float(i) for i in raw_value]
                elif len(data) == 2:
                    value = float(data[1])
                break
        return value
 

def filter_item(input_lst: list):
    """procedure to remove empty item in a list and remove the '\n' character."""
    new_lst = []
    for item in input_lst:
        new_item = item.rstrip()
        if new_item:
            new_lst.append(new_item)
    return new_lst

