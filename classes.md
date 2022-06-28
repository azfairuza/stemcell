# Classes

This are the classes used in the program

Some notations in this file are:
- `-`: private instance attribute
- `~`: private class attribute


```mermaid
classDiagram
    class Ligand{
        -float x_position
        -float y_position
        -bool bound
        -bool targeted
        -int id
        -list~int~ target_id
        ~int ligand_number
        resetNumber(cls)
    }
    class Nanopattern{
        -float height
        -float width
        -float grid_height
        -float grid_width
        -List~Position~ position_seed
        -float dot_size
        -int row_number
        -int col_number
        -string type_class
        -List~Ligand~ ligands
        getNotTargetedLigand(self) List~Ligand~
        getXPositionLigand(self) List~float~
        getYPositionLigand(self) List~float~
        getLigandById(self, id) List~float~
        show(self)
    }
    class Integrin{
        -int cell_id
        -int id
        -string type_class
        -float x_position
        -float y_position
        -bool surface
        -bool bound
        -bool targeting
        -int target_type
        -int cell_target_id
        -int integrin_target_id
        -int ligand_target_id
        -float x_target
        -float y_target
        -float mass
        ~int integrin_number
        resetNumber(cls)
        getInformation(self) string
        getTargetDistance(self) float
        getLigandDistance(self, ligand) float
        getIntegrinDistance(self, integrin) float
        targetingProcedure1(self, cells, substrate) target
        updatingProcedure1(self, target)
        movingProcedure1(self, cells, substrate, dstlimit, movespeed)
        move(self, movespeed)
        updateTarget(self, cells)
    }
    class Cell{
        -int id
        -string type_class
        -float mass
        -float x_center
        -float y_center
        -float integrin_size
        -float radius
        -List~Integrin~ integrins
        ~int cell_number
        getIntegrinInfo(self) string
        getXPositionIntegrin(self) List~float~
        getYPositionIntegrin(self) List~float~
        getIntegrinById(self, id) Integrin
        getAvailableIntegrin(self) List~Integrin~
        getCenterofMass(self) Position
        totalIntegrinBound(seld) int
        totalMass float
        show(self, Nanopattern)
        resetNumber(cls)
    }
    class Main{
        readFile(filename)
        getValue(lst_strng, property_name)
        filterElement(input_lst)
        getCellbyId(cells, id) Cell
        excludeCellById(cells, id) List~Cell~
        showAll(cells, Nanopattern, show_substrate, save, number, folder, line)
        saveCenterOfMass(cells, num_iteration) CELLCM.txt
        simulate1(cells, Nanopattern)
        circles(x, y, s, c='b', vmin, vmax, **kwargs)
    }
    
    direction TD
    Nanopattern "1" --* "1..*" Ligand : Contains
    Cell "1..*" --* "1..*" Integrin : Contains
    Main --* Cell : Procedure
    Main --* Nanopattern : Procedure
