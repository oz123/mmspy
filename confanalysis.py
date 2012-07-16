#!/usr/bin/env python
# -*- coding: utf-8 -*-  
#       This file is the main code of the conflictAnalysis tool developed in DSITE group
#       Center for Applied Geosciences, University of Tuebingen, Germany
#       The algorithm was created by Max Morio, and written in python by Oz Nahum
#       and Max Morio.
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.


"""
This is the main script that runs the conflict analysis.
most of the functios are in the module mmsca.py
"""
import stat,sys,os,shutil 
import mmsca
import numpy as np
from osgeo import ogr

def parseArgs():
    """
    Do some basic checks about the command line arguments.
    - Verify the Projekt.ini exists.
    - Verify that conflict type is specified.  
    """
    if len(sys.argv) < 2:
        print "No arguments given ... quiting ..."
        sys.exit(1)
    if not os.path.isdir(sys.argv[1]):
        print "ERROR: The first argument must be a directory which " \
        + "holds the Projekt.ini"
        sys.exit(1)
    try:
        open(os.path.abspath(sys.argv[1])+'/DATA/Projekt.ini')
    except IOError:
        print "ERROR: Could not find the Projekt.ini under " \
        + os.path.abspath(sys.argv[1])+'/DATA/'
    try: 
        conflicttype = sys.argv[2]
        try: 
            conflicttype=int(sys.argv[2])
            print "conflicttype ", conflicttype
        except ValueError:
            print "Error: Conflict Type must be an integer"
    except IndexError:      
        print "ERROR: You didn't specify conflict type."
        sys.exit(1)
    return conflicttype



def populateShpfileDbase(shpfile,zwerts):
    """
    Populate the database of the shapefile with allowed contaminants 
    thresholds for each land use.
    """
    for polygon in range(shpfile.NPolygons):
        pol = shpfile.layer.GetNextFeature()
        luc = shpfile.get_value(pol, "Kategorie")
        # land uses code range is from 100 to 900 which translates
        # to 1st or 9th column in allowed thresholds in zwert
        # keep in mind, column indecies start from 0 ... 
        luc = (int(luc) / 100) - 1
        for field, contname in zip(shpfile.fields[3:],
            zwerts.contnames):
            values=zwerts.targets_LUT[contname]
            value=values[luc]
            shpfile.set_value(pol, field, value)
    shpfile.layer.ResetReading()

def create_target(shpfile, proj, contname,xres, yres):
    """
    create a raster of thresholds for each contaminant
    the target rasters are created in the size of the bounding box 
    of area of interest
    """
    traster =  proj.aktscenario+'/'+proj.aktlayout+'/'
    target = shpfile.rasterize_field(xres,yres,fieldname=contname, 
    rasterfilepath = traster+contname+'_target.asc')
    print "Target Raster created in: ", os.path.abspath(traster+contname+'_target.asc')
    return target
    

def cut_to_cutline(proj,contraster):
    """
    cut raster to cutline
    """
    craster = mmsca.MaskRaster()
    craster.reader("DATA/"+contraster.replace('aux', 'asc'))
    xres, yres = craster.extent[1], craster.extent[1]
    craster.fillrasterpoints(xres, yres)
    craster.getareaofinterest("DATA/area_of_interest.shp")
    area_of_interest_polygon=craster.boundingvertices
    craster.clip_to_cutline(xres,yres)
    minX, maxY = craster.xllcorner, craster.yurcorner
    minX, maxY = round(minX,-1), round(maxY,-1)
    clipped = os.path.abspath( proj.aktscenario+'/'+proj.aktlayout+'/'+contraster.replace('.aux','_cutline.asc'))
    craster.data = np.ma.MaskedArray(craster.data, mask=craster.mask)
    craster.data = np.ma.filled(craster.data, fill_value=-9999)
    craster.writer(clipped,craster.data, (minX-xres, maxY+yres), xres, yres, Flip = False)
    print "Cutline  Raster created in: ", clipped    
    
    return craster 

def calculate_exceedance(proj, cont, cont_name, target):
    """
    devide the data from each cont raster by the 
    data from each target raster 
    """ 
    #for cont, targets in conts.keys():
    eraster = mmsca.MaskRaster()
    eraster.data = np.flipud(cont.data)/target.data
    eraster.data = np.flipud(eraster.data)
    # replace all negative values with zeros
    eraster.data = np.where(eraster.data > 0, eraster.data, 0)
    eraster.mask = cont.mask
    xres, yres =  cont.extent[1],  cont.extent[1]
    minX, maxY = cont.xllcorner, cont.yurcorner
    minX, maxY = round(minX,-1), round(maxY,-1)
    # put mask, e.g. remove all values outside of Area Of Interest
    eraster.data = np.ma.MaskedArray(eraster.data, mask=eraster.mask)
    # exceedance absolute is the ratio of contaminant level to allowed threshold
    exceedance = os.path.abspath( proj.aktscenario+'/'+proj.aktlayout+'/'+cont_name+'_exceedance_abs.asc')
    print "Absolute exceedance  raster created in: ", exceedance
    eraster.writer(exceedance, np.ma.filled(eraster.data, fill_value=-9999), (minX-xres, maxY+yres), xres, yres, Flip = False)
    # calucalte binary exceedance raster
    # if the ratio is higher than 1 than replace with True (1), if ratio is smaller
    # than replace with False (0)
    eraster.data = np.ma.where(eraster.data > 1, 1, 0)
    exceedance = os.path.abspath( proj.aktscenario+'/'+proj.aktlayout+'/'+cont_name+'_exceedance_bool.asc')
    print "Boolean exceedance  raster created in: ", exceedance
    eraster.writer(exceedance, np.ma.filled(eraster.data, fill_value=-9999), (minX-xres, maxY+yres), xres, yres, Flip = False)
    return eraster

def polygonize_exceedance(cont):
    pass
    
def main(conflicttype):
    """
    The main algorithm that drives the conflict analysis.
    """
    
    # initialize all project properties
    workdir=sys.argv[1]
    proj =  mmsca.Project()
    proj.getconfig(os.path.abspath(workdir)+"/DATA/Projekt.ini")
    proj.cleanup()
    print proj.aktscenario, proj.aktlayout
    layout = proj.aktscenario+'/'+proj.aktlayout+'.shp'
    layout_tgl = proj.aktscenario+'/'+proj.aktlayout+'/' \
        +proj.aktlayout+'_tgl'+'.shp'
    zwert = mmsca.ZielWerte(proj)
    
    # create shape file object , since we specify a copy
    # all further actions are preformed on the copy with _tgl suffix
    scenario = mmsca.LandUseShp(layout,layout_tgl)
    xres, yres = 10, 10 
    luraster =  proj.aktscenario+'/'+proj.aktlayout+'/'+proj.aktlayout+'.asc'
    # rasterize the layers
    layout_raster=scenario.rasterize_field(xres, yres, rasterfilepath=os.path.abspath(luraster))
    print "Land Uses Raster created in: ", os.path.abspath(luraster)
    # add column for each contaminant 
    for contaminant, component in zip(zwert.contnames, zwert.compartments):
        #print contaminant, component
        if component == "Boden":  scenario.addfield(contaminant+"_B")
        elif component == "GW": scenario.addfield(contaminant+"_in_GW")
    # for each polygon fill in the allowed threshold for each contaminant
    # based on the land use codeG
    scenario.layer.ResetReading()
    populateShpfileDbase(scenario,zwert) 
    # import pdb; pdb.set_trace()
    # create targets, cut raster to cutline, calculate the exceedances   
    # write ascii raster for each exceedance and convert that raster to
    # shape file
    for cont in zwert.contrasternames:
        cont_target = create_target(scenario, proj, cont.replace(".aux",""), xres, yres)
        cont_raster = cut_to_cutline(proj, cont)
        # calculate exceedance now works directly on the arrays in memory
        calculate_exceedance(proj, cont_raster,  cont.replace(".aux","") , cont_target)
        wdir= proj.aktscenario+'/'+proj.aktlayout#+'/'+proj.aktlayout
        layer = cont.replace(".aux","")+"_exceedance"
        # TOOD: add srs information !!!
        exceedance_shp = mmsca.ShapeFile(wdir,layer,fields={"ID":ogr.OFTInteger, "Exceedance":ogr.OFTInteger})
        exceedance_shp.dst_layer.SyncToDisk()
        exceedance_ras = mmsca.ASCIIRaster()
        exceedance_ras.polygonize("SzenarioA/ScALayout2/" \
                +cont.replace(".aux","_exceedance_bool.asc"), exceedance_shp.dst_layer, 1)
        print "**************\n"
        


if __name__ == '__main__':
    conflicttype=parseArgs()
    main(conflicttype)
    




#usage:
#:~/Desktop/PROJECT_MMSpy/Projekt$ ./confanalysis.py . 0
## output:

#conflicttype  0
#SzenarioA ScALayout2
#Anzahl Kontaminanten (B und GW):  4
#Shape file copied to : SzenarioA/ScALayout2/ScALayout2_tgl.shp
#Land Uses Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/ScALayout2.asc
#Target Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK30_B_target.asc
#Extent des Bewertungsgebiets (area_of_interest):
#(3366307.654296875, 3367789.362121582, 5813454.151306152, 5814844.734924316)
#UR: 3367789.36212 5814844.73492
#LL: 3366307.6543 5813454.15131
#Cutline  Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK30_B_cutline.asc
#Absolute exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK30_B_exceedance_abs.asc
#Boolean exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK30_B_exceedance_bool.asc
#created SzenarioA/ScALayout2/PAK30_B_exceedance
#**************

#Target Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK60_B_target.asc
#Extent des Bewertungsgebiets (area_of_interest):
#(3366307.654296875, 3367789.362121582, 5813454.151306152, 5814844.734924316)
#UR: 3367789.36212 5814844.73492
#LL: 3366307.6543 5813454.15131
#Cutline  Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK60_B_cutline.asc
#Absolute exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK60_B_exceedance_abs.asc
#Boolean exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PAK60_B_exceedance_bool.asc
#created SzenarioA/ScALayout2/PAK60_B_exceedance
#**************

#Target Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PCE_in_gw_target.asc
#Extent des Bewertungsgebiets (area_of_interest):
#(3366307.654296875, 3367789.362121582, 5813454.151306152, 5814844.734924316)
#UR: 3367789.36212 5814844.73492
#LL: 3366307.6543 5813454.15131
#Cutline  Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PCE_in_gw_cutline.asc
#Absolute exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PCE_in_gw_exceedance_abs.asc
#Boolean exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/PCE_in_gw_exceedance_bool.asc
#created SzenarioA/ScALayout2/PCE_in_gw_exceedance
#**************

#Target Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/TCE_in_gw_target.asc
#Extent des Bewertungsgebiets (area_of_interest):
#(3366307.654296875, 3367789.362121582, 5813454.151306152, 5814844.734924316)
#UR: 3367789.36212 5814844.73492
#LL: 3366307.6543 5813454.15131
#Cutline  Raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/TCE_in_gw_cutline.asc
#Absolute exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/TCE_in_gw_exceedance_abs.asc
#Boolean exceedance  raster created in:  /home/ozdeb/Desktop/PROJECT_MMSpy/Projekt/SzenarioA/ScALayout2/TCE_in_gw_exceedance_bool.asc
#created SzenarioA/ScALayout2/TCE_in_gw_exceedance
#**************
