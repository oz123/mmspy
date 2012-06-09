#!/usr/bin/env python
# -*- coding: utf-8 -*-  
'''
Created on 12.01.2010

@author: Max Morio
Creates a jpeg out of a layout shapefile using shp2img.exe the mapfile parts header.map and
tail.map 
'''
import sys, os
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

#contname, compartment, targ_comm, targ_resid, targ_recreat, targ_agric, contrastername, \
#    no_contaminants = readZielwert(aktlayout)
#Replaced 2011.2.09: contname,compartment,contrastername,no_contaminants , targets_LUT = readZielwert(aktlayout)
#contname,compartment,contrastername,no_contaminants , targets_LUT, boolConflictType, iContRealisations \
#		, no_contaminants, selcont = readZielwert(aktlayout)
#aktscenario, aktlayout, layoutdir , scen_landuseratio, aoiraster, LUTcolours, LUTName, n_Landuses, \
#			no_contaminants, selcont = readZielwert( aktlayout )
contname,compartment,contrastername,no_contaminants,targets_LUT,boolConflictType,iContRealisations = readZielwert(aktlayout)
#Call the map file creator:
mapFileName = CreateMAP(aktscenario, aktlayout, '')
os.system('shp2img -m ' + mapFileName + ' -o ../' + aktlayout + '.jpg')
print '\n----------------------------------------------\n ---> JPEG aus Layoutshapefile wurde erzeugt.  \
	\n \nFenster bitte schliessen, um zum MMT zurueckzukehren !!!'

