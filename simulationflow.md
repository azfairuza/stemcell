# Simulation Flowchart

This page contains the flowchart of the simulation.

## Procedure Mark-I
this flowchart is named as simulation number 1. This simulation account the surface integrin but have strict rules on targeting procedure. The ligand only can be targeted by single object, same as the integrin-integrin interaction. 

This flowchart occurs after all the elements has been created. The element are nanopattern (contains ligands) and 
cells (contain integrins)

```mermaid
  flowchart TD
  %% components of the flowchart
  checkBounding1{does the distance </br> between integrin and target </br> is close enough?}
  checkCell1{does all cells </br> have been picked?}
  checkCell2{does all cells </br> have been updated?}
  checkDistance1{does the distance </br> more than the limit?}
  checkIntegrin1{does each integrin </br> in the cell  </br> have a target?}
  checkIntegrin2{does each integrin </br> in the cell  </br> have moved?}
  compare1[find the nearest </br> between the integrin </br> and the ligand]
  End([End])
  getCell1[Pick a cell]
  getCell2[Pick a cell]
  getIntegrin1[Pick an integrin </br> from the cell]
  getIntegrin2[find nearest untargeted </br> integrin from other cells]
  getIntegrin3[Pick an integrin </br> from the cell]
  checkSurface{is the integrin </br> in the cell surface?}
  getLigand1[find nearest </br> untargeted ligand]
  getLigand2[find nearest </br> untargeted ligand]
  iterate1{iterate again?}
  Move1[Move the integrin] 
  Start([Start]) 
  update1[update the attribute </br> of main integrin </br> and targeted object]
  update2[update the attribute </br> of main integrin </br> and targeted object]
  update3[update the attribute </br> of main integrin </br> and targeted object]
  
  
  
  %% relations of the flowchart
  
  subgraph targetting
  Start --> getCell1
  getCell1 --> getIntegrin1
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
  checkCell1 --> |No| getCell1
  end
  
  subgraph moving
  checkCell1 ---> |Yes| getCell2
  getCell2 --> getIntegrin3
  getIntegrin3 --> checkDistance1
  checkDistance1 --> |Yes| checkIntegrin2
  checkDistance1 --> |No| checkBounding1
  checkBounding1 --> |Yes| update3
  update3 --> checkIntegrin2
  checkBounding1 --> |No| Move1
  Move1 --> checkIntegrin2
  checkIntegrin2 --> |No| getIntegrin3
  checkIntegrin2 --> |Yes| checkCell2
  checkCell2 --> |No| getCell2
  checkCell2 --> |Yes| iterate1
  iterate1 --> |Yes| getCell2
  iterate1 --> |No| End
  end  
```

These flowchart ia a simple mechanisme and has not include any law of physics

## Procedure Mark-II

This procedure is named by simulation2. The differences is this procedure allow multi-targeting to occur, in which a ligand can be targeted by many-objects. Another upgrade that is introduced in this procedure are as follows:
  1. The distance of targeting will be made short.
  2. Surface integrin prioritize integrin-integrin interaction, as long as below the limit distance, surface integrin doesn't search any nearest ligand afterwards.
  3. The integrin does validating routine which checks the target position (if the target is integrin) and also bound status of the target. If the target has bound to another object, the integrin will do the searching routine to find another nearest target. 
  4. The limit of surface can be set through celcon.txt
  5. Integrin who have not found any nearest object, will do random move. 

