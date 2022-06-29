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

```mermaid
  flowchart TD
  %% components of the flowchart
  bound1[/the integrin </br> bound to target/]
  checkBound1{is the integrin </br> bound to another </br> object?}
  checkBound2{is the integrin </br> bound to another </br> object?}
  checkNumberOfCell1{is the number if cell </br> is more than 1?}
  checkSurface1{is the integrin </br> at surface?}
  checkTarget1{does the integrin </br> have target?}
  End([end])
  findNearestIntegrin1[Find nearest </br> unbound-integrin </br> from other cell]
  findNearestLigand1[Find nearest </br> unbound ligand]
  foundTarget1[/Found a target/]
  foundTarget2[/Found a target/]
  getCell1[Pick an unchoosen cell]
  getCell2[Pick an unchoosen cell]
  getIntegrin1[Pick unchoosen integrin </br> from the cell]
  getIntegrin2[Pick unchoosen integrin </br> from the cell]
  iterCell1{Are all </br> of the cell </br> picked?}
  iterCell2{Are all </br> of the cell </br> picked?}
  iterIntegrin1{Are all </br> of the integrin </br> picked?}
  iterIntegrin2{Are all </br> of the integrin </br> picked?}
  iterSimulation1{iterate again?}
  moveIntegrin1[move 1 step </br> towards the target]
  removeTarget1[remove </br> current target]
  targeting1[find a target </br> through targeting </br> procedure]
  updateAttribute1[Update the integrin </br> and the target attribut]
  updateAttribute2[Update the integrin </br> and the target attribut]
  validateTarget1{does the target </br> bound to </br> another object?}
  validateBound1{is the distance </br> close enough?}
  
  start([Start])
  
  %% diagram scheme
  subgraph A [Targeting]
    start --> getCell1
    getCell1 --> getIntegrin1
    getIntegrin1 --> checkBound1
    checkBound1 --> |YES| iterIntegrin1
    checkBound1 --> |NO| checkNumberOfCell1
    checkNumberOfCell1 --> |YES| checkSurface1
    checkNumberOfCell1 --> |NO| findNearestLigand1
    checkSurface1 --> |YES| findNearestIntegrin1
    checkSurface1 --> |NO| findNearestLigand1
    findNearestIntegrin1 --> foundTarget1
    findNearestLigand1 --> foundTarget1
    foundTarget1 --> updateAttribute1
    updateAttribute1 --> iterIntegrin1
    iterIntegrin1 --> |YES| iterCell1
    iterIntegrin1 --> |NO| getIntegrin1
    iterCell1 --> |NO| getCell1
  end
  subgraph B [Move Procedure]
    iterCell1 --> |YES| getCell2
    getCell2 --> getIntegrin2
    getIntegrin2 --> checkBound2
    checkBound2 --> |YES| iterIntegrin2
    checkBound2 --> |NO| checkTarget1
    checkTarget1--> |YES| validateTarget1
    checkTarget1 --> |NO| targeting1
    validateTarget1 --> |YES| removeTarget1
    validateTarget1 --> |NO| moveIntegrin1
    moveIntegrin1 --> validateBound1
    validateBound1 --> |YES| bound1
    validateBound1 --> |NO| iterIntegrin2
    bound1 --> updateAttribute2
    removeTarget1 --> targeting1
    targeting1 --> foundTarget2
    foundTarget2 --> updateAttribute2
    updateAttribute2 --> iterIntegrin2
    iterIntegrin2 --> |YES| iterCell2
    iterIntegrin2 --> |NO| getIntegrin2
    iterCell2 --> |YES| iterSimulation1
    iterCell2 --> |NO| getCell2
    iterSimulation1 --> |YES| getCell2
    iterSimulation1 --> |NO| End
   end
    
```

