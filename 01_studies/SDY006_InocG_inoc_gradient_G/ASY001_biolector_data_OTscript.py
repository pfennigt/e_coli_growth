from opentrons import robot, labware, instruments, protocol_api, util
import numpy as np

# metadata
metadata = {
    'protocolName': 'OT_GlcGradientONECELL_BigFalconRack',
    'authors': 'Tobias Pfennig (topfe101@hhu.de)',
    'description': '\
       Protocol used for pipetting Biolecor flowerplates with each six replicate rows\
       of a Glucose gradient. Wells in the same column have equal carbon molarity.\
       All columns: glucose (0 mM - 19.98 mM)\
       Last column with highest glucose concenration\
       Last three wells in last column: Blank without E. coli and highest glucose concentration',
}


# ============================= custom functions =============================

def fillhgt_falcon50000(volume, diameter=27, cone_height=14.5, cone_tipdiameter=5):
    """
    Describes 50000 ul falcons geometrically and 
    calculates the respective fill-height to a fill-volume.
    """
    # measurements: 30 x 150 mm

    cone_slope = (diameter / 2 - cone_tipdiameter / 2) / cone_height
    cone_maxvol = 0.33 * np.pi * cone_height * ((diameter / 2) ** 2 + (diameter / 2) * (cone_tipdiameter / 2) + (cone_tipdiameter / 2) ** 2)
    cone_heighttip = cone_tipdiameter / (2 * cone_slope)
    cone_voltip = 0.33 * np.pi * cone_heighttip * (cone_tipdiameter / 2) ** 2

    ##catch errors##
    if volume > 50000:
        raise ValueError("volume exceeds 50 ml falcon volume")
    if volume < 0:
        raise ValueError("negative volumes not possible")

    ##distribute volume to modeled parts##
    #volume in cone
    if volume <= cone_maxvol:
        cone_vol = volume
    else:
        cone_vol = cone_maxvol
    
    #volumein cylinder
    cyl_vol = volume - cone_vol

    ##calculate height in cone##
    # volume: 1/3 * pi * h * r^2, r = m * h
    if cone_vol == cone_maxvol:
        fillhgt_cone = cone_height # if the cone is full, then the fill-height in the cone is the cone height
    else:
        cone_fullvol = volume + cone_voltip 
        fillhgt_conefull = ((3 * cone_fullvol) / (np.pi * cone_slope ** 2)) ** 0.33
        fillhgt_cone = fillhgt_conefull - cone_heighttip

    ##top cylinder##
    if cyl_vol == 0:
        fillhgt_cyl = 0
    else:
        fillhgt_cyl = cyl_vol / (np.pi * (diameter / 2) ** 2)

    ##add heights##
    fillhgt = fillhgt_cone + fillhgt_cyl
    
    return fillhgt


def fillhgt_falcon15000(volume, diameter=14, cone_height=21, cone_tipdiameter=4):
    """
    Describes 15000 ul falcons geometrically and 
    calculates the respective fill-height to a fill-volume.
    """
    # measurements: 20 x 120 mm

    cone_slope = (diameter / 2 - cone_tipdiameter / 2) / cone_height
    cone_maxvol = 0.33 * np.pi * cone_height * ((diameter / 2) ** 2 + (diameter / 2) * (cone_tipdiameter / 2) + (cone_tipdiameter / 2) ** 2)
    cone_heighttip = cone_tipdiameter / (2 * cone_slope)
    cone_voltip = 0.33 * np.pi * cone_heighttip * (cone_tipdiameter / 2) ** 2

    ##catch errors##
    if volume > 15000:
        raise ValueError("volume exceeds 15 ml falcon volume")
    if volume < 0:
        raise ValueError("negative volumes not possible")

    ##distribute volume to modeled parts##
    #volume in cone
    if volume <= cone_maxvol:
        cone_vol = volume
    else:
        cone_vol = cone_maxvol
    
    #volumein cylinder
    cyl_vol = volume - cone_vol

    ##calculate height in cone##
    # volume: 1/3 * pi * h * r^2, r = m * h
    if cone_vol == cone_maxvol:
        fillhgt_cone = cone_height # if the cone is full, then the fill-height in the cone is the cone height
    else:
        cone_fullvol = volume + cone_voltip 
        fillhgt_conefull = ((3 * cone_fullvol) / (np.pi * cone_slope ** 2)) ** 0.33
        fillhgt_cone = fillhgt_conefull - cone_heighttip

    ##top cylinder##
    if cyl_vol == 0:
        fillhgt_cyl = 0
    else:
        fillhgt_cyl = cyl_vol / (np.pi * (diameter / 2) ** 2)

    ##add heights##
    fillhgt = fillhgt_cone + fillhgt_cyl
    
    return fillhgt


def transfer_fillhgt(volume, source, dest, pipette, fillhgtfunc, fillvol, pipParams, **kwargs):
    """
    Changed version of default "transfer" function.
    Transfers static volumes from one source to one or more destinations.
    During transfer the fillheight in the source is calculated after every pipetting.
    This fillheight is used to ensure minimal tip immersion of around tipDepth.
    """
    
    global _new_fillvol
    
    liqHeightAdjust, tipDepth, lowestHeight = pipParams #unpack parameters
    
    source = source.top() # force inclusion of coordinates in transfer plan
    
    
    #adjust destination if vector is given
    destVector = None

    if type(dest) is tuple:
        if type(dest[1]) is util.vector.Vector:
            destVector = dest[1]
            
            dest=dest[0]
    
    
    #cited from "pipette.transfer", substituting 'self' with 'pipette'
    kwargs['mode'] = kwargs.get('mode', 'transfer')

    touch_tip = kwargs.get('touch_tip', False)
    if touch_tip is True:
        touch_tip = -1
    kwargs['touch_tip'] = touch_tip

    tip_options = {
        'once': 1,
        'never': 0,
        'always': float('inf')
    }
    tip_option = kwargs.get('new_tip', 'once')
    tips = tip_options.get(tip_option)
    if tips is None:
        raise ValueError('Unknown "new_tip" option: {}'.format(tip_option))

    plan = pipette._create_transfer_plan(volume, source, dest, **kwargs)
    #end citation
    
    
    #adjust pipetting plan
    for i in range(len(plan)):
        curr_vol = plan[i]["aspirate"]["volume"]
        
        fillvol -= curr_vol  # update fillvolume to after the following pipetting
        
        currdist = fillhgtfunc(fillvol) + liqHeightAdjust - tipDepth
        if currdist < lowestHeight: 
            currdist = lowestHeight 
        
        newLoc = (plan[i]["aspirate"]["location"][0],
                  util.vector.Vector(
                          plan[i]["aspirate"]["location"][1][0],
                          plan[i]["aspirate"]["location"][1][1],
                          currdist))
        
        plan[i]["aspirate"]["location"] = newLoc
        
        if destVector is not None:
            plan[i]["dispense"]["location"] = (plan[i]["dispense"]["location"], destVector)
        
        _new_fillvol = fillvol #save the updated fillvol as _new_fillvol to use further
        
    #cited from "pipette.transfer"
    pipette._run_transfer_plan(tips, plan, **kwargs)
    
    return(pipette)


# =============================== Load labware ===============================

tip_name = 'Biozym-tiprack-200ul'       # Add Biozym 200 ul tiprack to labware library
if tip_name not in labware.list():
    custom_tip = labware.create(
        tip_name,                       # name of you labware
        grid=(12, 8),                   # specify amount of (columns, rows)
        spacing=(9, 9),                 # distances (mm) between each (column, row)
        diameter=5.23,                  # diameter (mm) of each well on the plate
        depth=50,                       # depth (mm) of each well on the plate
        volume=200)

plate_name = 'biolector-flower-plate-48-well'
if plate_name not in labware.list():
    custom_plate = labware.create(
        plate_name,                    # name of labware
        grid=(8, 6),                    # amount of (columns, rows)
        spacing=(13, 13),               # distances (mm) between each (column, row)
        diameter=10.85,                     # diameter (mm) of each well on the plate
        depth=33,                       # depth (mm) of each well on the plate
        volume=3200)

tiprack = labware.load(tip_name, '9')   # Load tiprack

trough = labware.load('opentrons_6_tuberack_falcon_50ml_conical', '6')  # Load reservoir

troughCell = labware.load('opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', '5')  

pipette = instruments.P300_Single(mount='right', tip_racks=[tiprack], max_volume=200)  # Load pipette

plate = labware.load(plate_name, '3')  # Load wellplate


# =========================== Set source locations ============================

#Only use Ace cell culture -> Glc and Ace cultures are the same!!!
CellSource =  troughCell.wells("A1")          #50 ml falcon with dye added as baseline and  gradient
GlcSource =   trough.wells("B2")
#AceSource =   trough.wells("B3")
M9Source = troughCell.wells("A2")

pipette.start_at_tip(tiprack['A1'])


# ============================ Define Parameters =============================

liqHeightAdjust = 1.2                   # added to calculated fillheight to adjust for position of the falcon in mm (from fillhgt50000_calibrate)
tipDepth = 3                            # desired tipdepth in liquid AFTER aspiring in mm
lowestHeight = 5                        # lowest accepted tip height above the falcon bottom in mm

pipettingParameters = liqHeightAdjust, tipDepth, lowestHeight

sourceVolumes = {"CellSource"   : 10000,
                 "GlcSource"    : 50000,
                 #"AceSource"   : 30000,
                 "M9Source"     : 10000} #initial volume of the source falcons in ul


#%%
whichSteps = "all" #"all": complete experiment, "medium": all but cells, "cells": only cells


# ================================= Pipetting =================================
cellGradient = np.exp(np.linspace(np.log(0.02), np.log(0.5), 15))*1000 #cell volume for gradient, spaced evenly on log scle to investigate low and high 
cellGradient = np.append(cellGradient,0)

if whichSteps == "medium" or whichSteps == "all":
    
    # First distribute M9, then glucose/ actetate and cells last, best cooled
    pipette.pick_up_tip()
    
    # Distribute M9 in all since concentration isn't varied
    transfer_fillhgt(100,
                     M9Source,
                     plate.wells().top(),
                     pipette,
                     fillhgt_falcon15000,
                     sourceVolumes["M9Source"],
                     pipettingParameters,
                     new_tip = "once"
                     )
    sourceVolumes["M9Source"] = _new_fillvol
    
    
    # Distribute glucose everywhere, blanks treated through appended value
    for curRow in [(0,3),(1,4),(2,5)]:
        curWells = plate.rows(curRow[0]) + plate.rows(curRow[1])
        
        transfer_fillhgt((900-cellGradient),
                         GlcSource,
                         curWells.top(),
                         pipette,
                         fillhgt_falcon50000,
                         sourceVolumes["GlcSource"],
                         pipettingParameters,
                         new_tip = "once"
                         )
        sourceVolumes["GlcSource"] = _new_fillvol

if whichSteps == "all": robot.pause() #if the whole experiment is done, the calls can be cooled untl this point and then added into the robot

if whichSteps == "cells" or whichSteps == "all":

    # Distribute cells over all but last row
    for curRow in [(0,3),(1,4),(2,5)]:
        curWells = plate.rows(curRow[0]) + plate.rows(curRow[1])
        
        transfer_fillhgt(cellGradient,
                         CellSource,
                         curWells.top(),
                         pipette,
                         fillhgt_falcon15000,
                         sourceVolumes["CellSource"],
                         pipettingParameters,
                         new_tip = "once"
                         )
        sourceVolumes["CellSource"] = _new_fillvol