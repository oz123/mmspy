#!/usr/bin/env python
# -*- coding: utf-8 -*-  
'''
Created on 12.10.2010

@author: Max Morio
Creates a shapefile from a rasterfile using gdal_polygonize.py, rsp. 
the according function CreateTempExceedPolygon from tools.py
'''
import sys, os
try:
	from osgeo import gdal, ogr, osr
except ImportError:
    import gdal, ogr, osr

# get the path to the project-folder from commandline argument
# if no command-line-argument is passed, the script is expected to be located
# in the projects directory
try:
	if (len(sys.argv) > 1):
		projectinifile = sys.argv[1]
		projectinifile.replace(':', '//')
		projectinifile.replace('\\', '/')
except:
	print 'Kein Pfadargument bei Aufruf uebergeben. \n \
	 Nehme daher an, wir befinden uns im Verzeichnis des aktuellen Projekts! \
	 \n----------------------------------------------'
	# Define the path to the ini-file Projekt.ini
	projectinifile = os.getcwd() + '/DATA/Projekt.ini'

# this tools import * reads out the actual project's *.ini file and returns all relevant parameters
from tools import *

src_fn = [LAYOUTDIR + '/opttemp/optweighted/optmap_1_final.asc', 
		LAYOUTDIR + '/opttemp/optcost/optmap_1_final.asc', 
		LAYOUTDIR + '/opttemp/optMV/optmap_1_final.asc',
		LAYOUTDIR + '/opttemp/optSSI/optmap_1_final.asc',
		LAYOUTDIR + '/opttemp/optcompact/optmap_1_final.asc']

for i in range(0,len(src_fn)):
	# 1. remove optmap_1_final.shp, 2. remove sf.shp, 3. create new sf.shp, 4. cp sf.shp to optmap_1_final.shp
	try:
		# Remove the shapefiles named after the final map in the optimization out directory
		if os.path.isfile(src_fn[i].replace('.asc', '.shp')):
			os.remove(src_fn[i].replace('.asc', '.shp'))
		if os.path.isfile(src_fn[i].replace('.asc', '.shx')):
			os.remove(src_fn[i].replace('.asc', '.shx'))
		if os.path.isfile(src_fn[i].replace('.asc', '.dbf')):
			os.remove(src_fn[i].replace('.asc', '.dbf'))
		# Get the files from the outdir/test/ into the list
		for root,dir, files in os.walk(src_fn[i].replace('optmap_1_final.asc', '/test')):
			try:
				# print "Deleting /test/sf.shp from previous optimization runs ..."
				for name in files:
					suffix = os.path.splitext(name)
					os.remove(os.path.join(root,name)) #, src_fn[i].replace('.asc', suffix[1]))
				print "Convert opt-run's raster to shapefiles..."
			except: print 'Error in deleting old or converting new shapefiles....'
		CreateTempExceedPolygon(src_fn[i], src_fn[i].replace('optmap_1_final.asc', '/test'), 'sf', 'Kategorie')
		# Get the updated list again
		for root,dir, files in os.walk(src_fn[i].replace('optmap_1_final.asc', '/test')):
			try:
				for name in files:
					suffix = os.path.splitext(name)
					# Move the files from /test to the optimization out dir and name them after the final asc map
					os.rename(os.path.join(root,name), src_fn[i].replace('.asc',  suffix[1]))
			except:
				print 'Erzeugen/Ueberschreiben der Karten fehlgeschlagen.'
	except:
		print 'Fehler bei Erzeugen oder Kopieren der Vektorkarten (Shapefiles)!'

