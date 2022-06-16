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
        -bool bound_status
        -bool targeted_status
        -int ligand_id
        -int integrin_id
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
        -List~Ligand~ ligand
        getNotTargetedLigand(self) List~Ligand~
        getXPositionLigand(self) List~float~
        getYPositionLigand(self) List~float~
        show(self)
    }
    class Integrin{
        -float x_position
        -float y_position
        -bool surface
        -bool bound_status
        -int object_type
        -int cell_target_id
        -int integrin_target_id
        -int ligand_target_id
        -float x_target
        -float y_target
        -float mass
        -int cell_id
        -int integrin_id
        ~int integrin_number
        resetNumber(cls)
        getInformation(self) string
    }
    class Cell{
        -float mass
        -float x_center_of_mass
        -float y_center_of mass
        -float integrin_size
        -float radius
        -List~Integrin~ integrin
        ~int cell_number
        getIntegrinList(self) string
        getXPositionIntegrin(self) List~float~
        getYPositionIntegrin(self) List~float~
        show(self, Nanopattern)
        resetNumber(cls)
    }
