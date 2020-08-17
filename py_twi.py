'''
This script calculates the Topographic Wetness Index (TWI)
of a given Area of Interest.
'''

import os
import arcpy
from arcpy.sa import*
import math

def getDEM (roi):
    '''
    This function returns the DEM file that will be used to calculate the TWI.
    It will check if LIDAR coverage is available for a user provided ROI otherwise uses
    the TRIM DEM.
    '''
    # location of LIDAR and TRIM DEMs
    lidar_dir = r'\...\LiDAR\LiDAR_Products.gdb'
    lidar_outline = os.path.join (lidar_dir, 'lidar_utm_outline')
    dem_lidar = os.path.join (lidar_dir, 'dem_mosaic')
    dem_trim = r'\....\trim\...\bc_elevation_25m_bcalb.tif'

    # Check if the ROI is fully, partially or not covered by LIDAR
    area_roi = sum([row[0] for row in arcpy.da.SearchCursor(roi, ("SHAPE@AREA"))])

    out_intersect = r'in_memory/intersect'
    arcpy.Intersect_analysis([roi,lidar_outline], out_intersect)
    area_intersect = sum([row[0] for row in arcpy.da.SearchCursor(out_intersect, ("SHAPE@AREA"))])

    coverage_percent = int((area_intersect/area_roi)*100)

    # Decide which DEM to use and clip to ROI extent
    print ("Clipping DEM to ROI extent")
    if coverage_percent == 0:
        print ('Your ROI is not covered by LIDAR.')
        print ('TRIM DEM is used instead.')
        dem_roi = ExtractByMask(dem_trim, roi)

    elif coverage_percent == 100:
        print ('Your ROI is fully covered by LIDAR.')
        print ('LIDAR DEM is used.')
        dem_roi = ExtractByMask(dem_lidar, out_intersect)

    else:
        print ('Your ROI is partially ({}%) covered by LIDAR' .format(coverage_percent))
        print ('LIDAR DEM is used anyway.')
        dem_roi = ExtractByMask(dem_lidar, out_intersect)

    return dem_roi

def calculateTWI (dem):
    '''
    This function will calculate the TWI based on DEM dataset
    '''
    # Calculate the Flow Direction and Accumulation
    print ("Calculating Flow Direction and Accumulation...")
    flow_dir = FlowDirection(dem)
    flow_acc = FlowAccumulation(flow_dir)

    # Calculate the Slope and convert from Degree to Radians
    print ("Calculating Slope in radians...")
    slope_deg = os.path.join (arcpy.env.scratchWorkspace, 'slope_deg') # could'nt write to in_memory (too big?!)
    arcpy.Slope_3d (dem, slope_deg, "DEGREE")

    slope_rad = Raster(slope_deg) * (math.pi/180)

    # Get the DEM pixel area
    PxsizeX = float(arcpy.GetRasterProperties_management(dem, "CELLSIZEX").getOutput(0))
    PxsizeY = float(arcpy.GetRasterProperties_management(dem, "CELLSIZEY").getOutput(0))
    pxArea = PxsizeX*PxsizeY

    # Apply the TWI equation
    print ("Calculating TWI...")
    c = 0.0001
    exp_1 = Ln(((flow_acc + 1) * pxArea) / Tan(slope_rad))
    exp_2 = Ln(((flow_acc + 1) * pxArea) / c + Tan(slope_rad))

    TWI = Con(slope_rad > 0, exp_1, exp_2)

    # Apply a low-pass filter to reduce noise. Only if the DEM Cell size is small (very high resolution)
    if (PxsizeX <= 2 or PxsizeY <= 2 ):
        print ("The DEM cell size is {}. A low-pass filter is applied." .format(int(PxsizeX)))
        TWI_filtered = FocalStatistics(TWI, NbrRectangle(3, 3, "CELL"), "MEAN")
        TWI_filtered.save(os.path.join (Workspace, 'p_TWI_filtered'))

    else:
        TWI.save(os.path.join (Workspace, 'p_TWI'))

#----------#
# Run!
#----------#

Workspace = r'\...\TWI.gdb'
ROI = os.path.join (Workspace, 'ROI_demo')

def main():
    arcpy.env.overwriteOutput=True
    arcpy.env.scratchWorkspace = arcpy.env.scratchGDB
    arcpy.CheckOutExtension("Spatial")
    arcpy.CheckOutExtension("3D")

    dem_roi = getDEM (ROI)
    calculateTWI (dem_roi)

    print ("Deleting temporary files...")
    arcpy.Delete_management("in_memory")
    #arcpy.Delete_management(arcpy.env.scratchWorkspace)

    print ('Processing Completed!')

if __name__ == "__main__":
    main()
