#!/usr/bin/env python
# -*- coding: utf-8 -*-  
#        This file is the main code of the conflictAnalysis tool developed in DSITE group
#        Center for Applied Geosciences, University of Tuebingen, Germany
#        The algorithm was created by Max Morio, and written in python by Oz Nahum
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

def clipPollutionRasters(proj,zwert):
    """
    Iterate over the list of pollution rasters and clip them using 
    original demensions.
    """
    conts = {}
    for i, (contraster,contname) in enumerate(zip(zwert.contrasternames, zwert.contnames)):
    # clip each pollution raster so it is in the size of our masks
        if i == 0:
            craster = mmsca.MaskRaster()
            craster.reader("DATA/"+contraster.replace('aux','asc'))
            xres, yres = craster.extent[1], craster.extent[1]
            craster.fillrasterpoints(xres, yres)
            craster.getareaofinterest("DATA/area_of_interest.shp")
            craster.clip2(new_extent_polygon=craster.boundingvertices)
            
            area_of_interest_polygon=craster.boundingvertices
            craster.data = np.ma.MaskedArray(craster.data, mask=craster.mask)
            craster.data = np.ma.filled(craster.data, fill_value=-9999)
            clipped = os.path.abspath( proj.aktscenario+'/'+proj.aktlayout+'/'+contname+'_clipped.asc')
            craster.writer(clipped, craster.data, (craster.extent[0], craster.extent[3]+yres*0.5),10,10)
            print "Clipped Raster created in: ", clipped    
            craster.writer(clipped, craster.data, (craster.extent[0], craster.extent[3]+yres*0.5),10,10,Flip=False)
            conts[contname] = craster 
        else:
            craster.reader("DATA/"+contraster.replace('aux','asc'))
            craster.clip2(new_extent_polygon=area_of_interest_polygon)
            craster.data = np.ma.MaskedArray(craster.data, mask=craster.mask)
            craster.data = np.ma.filled(craster.data, fill_value=-9999)
            clipped = os.path.abspath( proj.aktscenario+'/'+proj.aktlayout+'/'+contname+'_clipped.asc')
            craster.writer(clipped, craster.mask, (craster.extent[0], craster.extent[3]+yres*0.5),10,10)
            print "Clipped Raster created in: ", clipped    
            craster.writer(clipped, craster.data, (craster.extent[0], craster.extent[3]+yres*0.5),10,10,Flip=False)
            conts[contname] = craster 
    return conts
    
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
    # based on the land use code
    scenario.layer.ResetReading()
    populateShpfileDbase(scenario,zwert)
    
    # create a raster of thresholds for each contaminant
    traster =  proj.aktscenario+'/'+proj.aktlayout+'/'+'target_'
    targets = {}
  
    for field, contname  in zip(scenario.fields[3:],  zwert.contnames) :
        targets[contname] = scenario.rasterize_field(xres,yres,fieldname=field, 
        rasterfilepath = traster+contname+'.asc')
        print "Target Raster created in: ", os.path.abspath(traster+contname+'.asc')
    
    
    cont_rasters=clipPollutionRasters(proj,zwert)
    craster=cont_rasters["PCE"]
    craster.clip_to_cutline(xres,yres)
    minX,maxY = craster.xllcorner, craster.yurcorner
    minX , maxY = round(minX,-1), round(maxY,-1)
    craster.writer("ccdata2m_cutline.asc",craster.data, (minX-xres, maxY+yres), 10,10,Flip=False)
    
if __name__ == '__main__':
    conflicttype=parseArgs()
    main(conflicttype)
    




#usage:
#:~/Desktop/PROJECT_MMSpy/Projekt$ ./confanalysis.py . 0
# output:
#conflicttype  0
#SzenarioA ScALayout1
#Shape file copied to : SzenarioA/ScALayout1/ScALayout1_tgl.shp
#Anzahl Kontaminanten (B und GW):  4
