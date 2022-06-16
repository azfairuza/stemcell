# Simulation Flowchart

This page contains the flowchart of the simulation. Because the main program has not had any simulation procedure, 
this flowchart will be the planned procedure for the future.


This flowchart occurs after all the elements has created. The element are nanopattern (contains ligands) and 
cells (contain integrins)

```mermaid
  flowchart TD
  %% components of the flowchart
  Start([Start]) 
  End([End])
  getCell[Pick a cell]
  getIntegrin1[Pick an integrin from the cell]
  checkSurface{is the integrin in the cell surface?}
  getIntegrin2[find nearest untargeted integrin from other cells]
  getLigand1[find nearest untargeted ligand]
  getLigand2[find nearest untargeted ligand]
  compare1[find the nearest between the integrin and the ligand]
  update1[update the attribute of main integrin and targeted object]
  update2[update the attribute of main integrin and targeted object]
  checkIntegrin1{does each integrin in the cell have a target?}
  checkCell1{does all cell have been picked?}
  
  
  %% relations of the flowchart
  Start --> getCell
  getCell --> getIntegrin1
  getIntegrin1 --> checkSurface
  checkSurface --> |Yes| getIntegrin2
  getIntegrin2 --> getLigand1
  getLigand1 --> compare1
  compare1 --> update1
  update1 --> checkIntegrin1
  checkIntegrin1 --> |No| getIntegrin1
  checkSurface --> |No| getLigand2
  getLigand2 ----> update2
  update2 --> checkIntegrin1
  checkIntegrin1 --> |Yes| checkCell1
  checkCell1 --> |No| getCell
  
  
  
  
  
```
