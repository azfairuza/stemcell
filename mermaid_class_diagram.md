# Classes

This are the classes used in the program

Some notation in this file are:
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
        getNotTargetedLigand(self): List~Ligand~
        getXPositionLigand(self): List~float~
        getYPositionLigand(self): List~float~
        show(self)

    }
