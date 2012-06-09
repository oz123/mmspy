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

import stat,sys,os,shutil #commands #stat,os,sys, 
#import shutil #pdb

# get the path to the project-folder from commandline argument
# if no command-line-argument is passed, the script is expected to be located
# in the projects directory and targets table (Zielwertset.rtv) is used

# script runs with 
# ./MMSConflictAnalysis.py . 1

print len(sys.argv)
try:
	print "inside"
        if (len(sys.argv) > 1):
		pathtoproject = sys.argv[1]
		conflicttype = sys.argv[ 2 ]
                print "conf", conflicttype 
		# get integer 0 for target based conflict analysis 1 for exceedance raster based
		pathtoproject.replace(':', '//')
		pathtoproject.replace('\\', '/')
                print pathtoproject
        os.chdir(pathtoproject)
except:
    print 'FEHLER: Datei zum Projektpfad nicht uebergeben oder nicht vorhanden! \n ####################################'
    sys.exit()

print 'Pfad: ', pathtoproject
print conflicttype

# put the project ini reading stuff here:
# Define the path to the ini-file 'Projekt.ini
projectinifile = os.getcwd() + '/DATA/Projekt.ini'
riskfilespath = os.getcwd() + '/DATA/RISK/'

# Read ./DATA/Project.ini, returning:
# aktscenario, aktlayout, LAYOUTDIR, aoiraster, LUTcolours, LUTName, n_landuses, no_contaminants, selcont
from tools import * #this already does all the reading of the ini file
#print selcont
dataSource = aktscenario + '/' + aktlayout + '.shp' #'SzenarioA/ScALayout1.shp'

#Max: added these two functions in order to erase results
#    from previous runs. Makes to code more stable on Windows
#    and also avoids to keep relictic results from previous runs
#    of the conflict analysis with different settings.
#    Also backs up the *.cost / *cosk , *snh and *WE files
mmsfile2keep = [aktscenario + '/' + aktlayout + '/' + aktlayout + '.cost'\
              , aktscenario + '/' + aktlayout + '/' + aktlayout + '.cosk'\
              , aktscenario + '/' + aktlayout + '/' + aktlayout + '.snh'\
              , aktscenario + '/' + aktlayout + '/' + aktlayout + '.WE'\
              , aktscenario + '/' + aktlayout + '/' + 'opttemp']
for x in mmsfile2keep:
    if os.path.isfile(x):
        shutil.copy2(x, x.replace(aktlayout + '/' + aktlayout, aktlayout)) #, ignore_errors=False, onerror=handleRemoveReadonly)
    elif os.path.isdir( x ): 
        shutil.copytree( x, x.replace( aktlayout + '/opttemp', '/opttemp' ))
#delete the dir including files and subdirs
shutil.rmtree(aktscenario + '/' + aktlayout, ignore_errors=True)#, onerror=handleRemoveReadonly)
#create the layout dir, again
if not os.path.isdir(aktscenario + '/' + aktlayout):
    os.mkdir(aktscenario + '/' + aktlayout)
#move back the *.cost, *cosk , *snh and *WE files to the aktlayout dir
for x in mmsfile2keep:
    if os.path.isfile(x.replace(aktlayout + '/' + aktlayout, aktlayout)):
        #shutil.move(x.replace(aktlayout + '/' + aktlayout, aktlayout), x)
        shutil.copy2(x.replace(aktlayout + '/' + aktlayout, aktlayout), x)
    elif os.path.isdir( x.replace(aktlayout + '/opttemp', '/opttemp' )):
		shutil.copytree(x.replace(aktlayout + '/opttemp', '/opttemp' ), x)

#pdb.set_trace()
driver = ogr.GetDriverByName('ESRI Shapefile')
dataSource = driver.Open(dataSource, 0)
layer = dataSource.GetLayer()
LandUses, Codes, Boundaries = getLayerProperties(layer)
#in order to read the rasters of contaminations we need to read the ZielWert file
#Which stores this info ! This is done by the call to function readZielwert(aktlayout)
#contname, compartment, targ_comm, targ_resid, targ_recreat, targ_agric, contrastername, \
# Max, 8.Dec2010: Added two new parameters for probabil. handling...
contname,compartment,contrastername,no_contaminants,targets_LUT,boolConflictType,iContRealisations = readZielwert(aktlayout)
#, no_contaminants,selcont
print 'Beruecksichtigte Schadstoffe: ', contrastername
#from getMask3 import getMask
if conflicttype == 0:
	print '\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
	print 'Konfliktanalyse wird anhand nutzungsspezifischer Zielwerttabelle ausgefuehrt!\n'
	print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
else:
	print '\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
	print 'Konfliktanalyse wird anhand risikobasierter Expositionsszenarien ausgefuehrt!\n'
	print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'

contSource = os.getcwd() + '/DATA/' + contrastername[0].replace('aux', 'asc') 
#here naming of rasters is case sensitive. If name in Zielwert.rtv is PAK_B.asc, a file pak_b.asc will result in error


#Max. 2009-12-15: Commented out the (re)recreation of the area_of_interest.asc,
# since we assume, it is alredy there!
# aoiraster is in global scope (returned from tools import)
#===============================================================================
# try:
#    aoits = os.stat('DATA/area_of_interest.shp')
#    aoiTime = aoits.st_mtime
#    aoirs = os.stat('DATA/area_of_interest.asc')
#    aoirsTime = aoirs.st_mtime
#    if aoiTime > aoirsTime:
#        print "Area of interest changed, recreating mask raster..."
#        nameOfMask = getMask(aoiraster.replace('asc', 'shp'), contSource)
#        RasterSource = nameOfMask
#    else: 
#        RasterSource = 'DATA/' + aoiraster
#        nameOfMask = RasterSource
# except:
#    print  "Could not find mask raster, creating now..."
#    nameOfMask = getMask(aoiraster.replace('asc', 'shp'), contSource)
#    RasterSource = nameOfMask
#===============================================================================
#RasterSource = 'DATA/mask_' + aoiraster
nameOfMask = 'DATA/mask_' + aoiraster
#nameOfMask = RasterSource

#Make sure that landUses raster is older than Projekt.ini. If older than project.ini
#re-create the landuse raster !
from time import *
a = time()
print '\n-------------------------------------------'

try:
    luts = os.stat(aktscenario + '/' + aktlayout + '.asc')
    lutsTime = luts.st_mtime
    pits = os.stat("DATA/Projekt.ini")
    pitsTime = pits.st_mtime #unix time of the Project.ini
    #if pitsTime>lutsTime:
    if pitsTime - 90 > lutsTime: #since we are writing the Project.ini section about 90 sec later
        print "Projekt/Layoutdaten geaendert!? Erzeuge neues Layout Raster..."
        #Oz, 23.03.2010: Removed obsolete call to RasterSource. 
        #createLandUsesRaster(Boundaries, Codes, LandUses, RasterSource, aktlayout, nameOfMask)
        createLandUsesRaster(Boundaries, Codes, LandUses, aktlayout, nameOfMask)
    else:
        print "Einstellungen bzgl. Konfliktanalyse unveraendert. \n Skip Erzeugung neues Layout Raster"
except:
    #print "Konnte kein Layout Raster finden... \nErzeuge ein Raster aus dem Layout Shapefile!"
    if not os.path.exists(aktscenario + '/' + aktlayout + '/'):
        os.mkdir(aktscenario + '/' + aktlayout)
    #createLandUsesRaster(Boundaries, Codes, LandUses, RasterSource, aktscenario, aktlayout, nameOfMask)
    createLandUsesRaster(Boundaries, Codes, LandUses, aktscenario, aktlayout, nameOfMask)

#nameOfMask='DATA/area_of_interest.asc'
#nameOfLandUsesRaster=aktscenario+'/'+aktlayout+'/'+aktlayout+'.asc'
nameOfLandUsesRaster = aktscenario + '/' + aktlayout + '.asc'
###add Project.ini section
projectFile = open(projectinifile, 'r')
config.readfp(projectFile)
#s=config.sections()
#print s
if config.has_section(aktscenario + '_' + aktlayout):
    config.remove_section(aktscenario + '_' + aktlayout)
    config.add_section(aktscenario + '_' + aktlayout)
else:
    config.add_section(aktscenario + '_' + aktlayout)

#add layout raster name to Project.ini
config.set(aktscenario + '_' + aktlayout, 'layout', nameOfLandUsesRaster)
#layoutshape=aktscenario+'/'+aktlayout + '.shp'
# Create land-use shapefile with targetvalues ( tablebased, in if clause is TRUE )
if conflicttype == '0':
	read_landuseshape(aktscenario,aktlayout, no_contaminants, contname, \
		compartment, targets_LUT, contrastername)
	###Create Target value rasters
	for contNumber in range(no_contaminants):
		CreateTargets(nameOfLandUsesRaster, contname, compartment, targets_LUT, contrastername, contNumber, \
				aktscenario, aktlayout)
	print '|\n --->\tZielwertraster fertiggestellt! \n'
else:
	# use targetvalues of '1' for all land uses, since risk-based conflict-analysis uses exceedance rasters
	# as input for conflict analysis
	read_landuseshape(aktscenario,aktlayout, no_contaminants, contname, \
		compartment, targets_LUT, contrastername)
	###Create Target value rasters
	targets_LUTones =  targets_LUT[ : ]
	#for i,value in enumerate( targets_LUTones ):
	#	 value.replace('\d', '1')
	for index,item in enumerate(targets_LUTones):
		for index2, item2 in enumerate( item ):
			item2 = '1'
			targets_LUTones[ index][ index2 ] = '1'
	print 'contrastername: ', contrastername
	for contNumber in range(no_contaminants):
		CreateTargets(nameOfLandUsesRaster, contname, compartment, targets_LUTones, contrastername, contNumber, \
				aktscenario, aktlayout)
	print "I was Here"
	# Create layout-specific targetrasters based on uniform use-type exceedance rasters...
	CreateTargetsRisk( nameOfLandUsesRaster, contname, compartment, targets_LUTones, contrastername, contNumber, \
			aktscenario, aktlayout)
	print '|\n --->\tLayoutspezifisches Zielwertraster fertiggestellt! \n'	

exidanceLayers = []
iRunsPerContaminant = []

for contNumber in range(no_contaminants):
	# New, Max, 10.12.2010: Define cases for deterministic or probabilistic conflict analysis
	nameOfPollutionRaster = []
	iRunsPerContaminant.append( 1 )
	#print contname
	#print boolConflictType # this parameter checks for probabilistic type of conflict analysis
	if boolConflictType[ contNumber ] =="1":
		# Perform probabilistic conflict analysis:
		iRunsPerContaminant[ contNumber ] = int( iContRealisations[ contNumber] )
		if conflicttype == '0':
			for n in range( 0, int( iContRealisations[ contNumber ]) ):
				#n denotes the realization suffix of the contrastername
				nameOfPollutionRaster.append('DATA/'+contrastername[contNumber].replace('.aux',str(n)+'.asc' ))
		else:
			for n in range( 0, int( iContRealisations[ contNumber ]) ):	
				#for m in range( 1, 9 ):
				#	nameOfPollutionRaster.append( 'DATA/RISK/lay_'+contrastername[contNumber].replace('.aux',str(n)+'.asc'))
				nameOfPollutionRaster.append( 'DATA/RISK/lay_'+contrastername[contNumber].replace('.aux',str(n)+'.asc'))

	elif boolConflictType[ contNumber ] =="0":
		if not iContRealisations[ contNumber ].isdigit():
			if conflicttype == '0':
				# Perform deterministic analysis on one single realization which has no integer suffix 
				nameOfPollutionRaster.append('DATA/' + contrastername[contNumber].replace('aux', 'asc'))
			else:
				nameOfPollutionRaster.append( 'DATA/RISK/lay_'+contrastername[contNumber].replace('.aux','.asc'))


		elif iContRealisations[ contNumber ].isdigit():
			if conflicttype == '0':
				iRunsPerContaminant[ contNumber ] = int(  iContRealisations[ contNumber ] )
				# Perform detterministic analysis on one of the realizations named by the iContRealisations suffix
				nameOfPollutionRaster.append('DATA/'+contrastername[contNumber].replace('.aux',str(iContRealisations[ contNumber ])+'.asc' ))
			else:
				nameOfPollutionRaster.append('DATA/RISK/lay_'+contrastername[contNumber].replace('.aux',str(iContRealisations[ contNumber ])+'.asc' ))
			#
			# TODO : check else part above...	
			#
	print 'nameOfPollutionRaster: ', nameOfPollutionRaster
	

	#print "RunsPerContaminant: " , iRunsPerContaminant[ contNumber ]
	#print "nameOfPollutionRaster: " , nameOfPollutionRaster
	for i in range( 0, iRunsPerContaminant[ contNumber ] ):
		# Perform the analysis part below here with nameOfPoll... as a list...

		#commented out, 10.12.2010: nameOfPollutionRaster = 'DATA/' + contrastername[contNumber].replace('aux', 'asc')
		#print nameOfMask, nameOfPollutionRaster
		if conflicttype == '0':
			nameOfClippedRaster = ClipPollution2(nameOfMask, nameOfPollutionRaster[i] )
			threshold = aktscenario + '/' + aktlayout + '/' + contname[contNumber] + '_target.asc'
		else:
			nameOfClippedRaster = ClipPollution2(nameOfMask, nameOfPollutionRaster[i] )
			threshold = aktscenario + '/' + aktlayout + '/' + contname[contNumber] + '_target.asc'
		if iRunsPerContaminant[ contNumber ] > 1:
			if i==0:
				nameOfContaminant = contname[contNumber] + str( i )
				contrastername[ contNumber ]= contrastername[ contNumber ].replace( '.aux', str( i )+'.aux' )
			else:
				nameOfContaminant = contname[ contNumber ] + str( i )
				contrastername[ contNumber ]= contrastername[ contNumber ]\
						.replace( str(i-1)+'.aux', str( i )+'.aux' ) # + str( i )
			#print '\n\n',contrastername[ contNumber ], '<- contraster | contname ->', contname[ contNumber ]
			
			#contrastername[ contNumber ]= contrastername[ contNumber ].replace( '.aux', str( i )+'.aux' ) # + str( i )
		else:
			nameOfContaminant = contname[ contNumber ] 
		#print  "|\n ### contrastername", contrastername[ contNumber ]
		#print " ### nameOfContaminant", nameOfContaminant
		#print " ### nameOfPollutionRaster[i]", nameOfPollutionRaster[ i ]
		
		exidance_raster, exidance_raster01s = CalculateExceedances(threshold, nameOfMask, nameOfClippedRaster, \
		nameOfContaminant, compartment[contNumber], aktscenario, contrastername[contNumber] )
		#add name of exidance raster to Project.ini
		opt = 'Zielwertueberschreitung(' + str(contNumber) + ')'
		config.set(aktscenario + '_' + aktlayout, opt, exidance_raster)
		if i == 0: os.chdir(aktscenario + '/' + aktlayout + '/' + 'TEMP/')
		else: os.chdir( aktscenario + '/' + aktlayout + '/' + 'TEMP/' )
		print 'CurDir: ', os.getcwd()
		if boolConflictType[ contNumber ]=="1":
			exidance_shp = contname[contNumber]+str(i)
		else: exidance_shp = contname[ contNumber ]
		exidance_raster01s = exidance_raster01s.replace(aktscenario + '/' + aktlayout + '/' + 'TEMP/', '')
		print 'Shapefiles mit Zielwertueberschreitungen: ', exidance_shp
		print "Test: ", exidance_raster01s, " ", exidance_shp, " ", exidance_shp
		CreateTempExceedPolygon(exidance_raster01s, exidance_shp, exidance_shp)
		#pdb.set_trace()
		#tempShp=aktscenario+'/'+aktlayout+'/TEMP/'+exidance_shp + '/' + exidance_shp +'.shp'
		tempShp = exidance_shp + '/' + exidance_shp + '.shp'
		print '|\n --->\t Erzeuge Polygone mit Zielwertueberschreitung fuer: ', nameOfContaminant
		dst_filename = CreatContExceedance(tempShp, aktscenario, aktlayout, nameOfContaminant, compartment[contNumber])
		exidanceLayers.append(dst_filename)

#b = time() - a
#print b, ' seconds to run so far ...'

###write new Project.ini
projectFile = open(projectinifile, 'w') #w mode overwrites original file
config.write(projectFile)

print 'Polygonlayers fuer Visualisierung: ', exidanceLayers

# TODO: Dec. 2010: Create JPEG maps for all realizations....
# Create the map file for output visualization
mapFileName = CreateMAP(aktscenario, aktlayout, exidanceLayers)
n =0 # counter to properly cycle through exceedance layers for JPEG creation below:
for i in range(no_contaminants):
	#print 'i: ',i
	#print "iRunsPerContaminant[i]: ", iRunsPerContaminant[ i ]
	for j in range( iRunsPerContaminant[ i ] ):
		#print 'j: ', j
		#print 'contrastername[i] in j-loop: ', contrastername[ i ]
		contname = contrastername[ i ].replace((str(iRunsPerContaminant[i]-1)+'.aux' ), str( j ))
		contname = contrastername[ i ].replace(str(j)+'.aux', str( j ))
		print contname, " ", str(iRunsPerContaminant[i]-1)+'.aux'
		if j > 0:
			#print "if j> 0"
			# This is the probabilistic case ( many realizations present )
			contname = contrastername[i].replace((str(iRunsPerContaminant[i]-1)+'.aux' ),  str( j ))
			contname = contrastername[ i ].replace(str(j)+'.aux', str( j ))
		else: 
			contname = contrastername[ i ].replace((str(iRunsPerContaminant[i]-1)+'.aux' ), str( j ))
			contname = contrastername[ i ].replace(str(j)+'.aux', str( j ))
		contname = contrastername[ i ].replace('.aux', '')
		contname= contname.replace(  ( str(iRunsPerContaminant[i]-1) ), str( j ) )
		print 'contname: ', contname
#TODO: if iRunsPerContaminant > 1: loop for contrastername
		if compartment[i] == 'Boden':
			contname = 'B_' + contname
		else:
			contname = 'W_' + contname
		print 'exidanceLayers[n] :', exidanceLayers[ n ]
		os.system('shp2img -m ' + mapFileName + ' -l ' + exidanceLayers[n] + ' -o '	+ contname + '.jpg')
		n = n+1
	#print 'Ueberschreitungskarte ',i,':', contname[2:], ': ', aktscenario + '/'	+ aktlayout + '/' + contname + '.jpg'

# converting all output to DOS-type ascii files (replace "LF" by "CRLF" and line ends)
files2convert = os.listdir(os.getcwd()+'/'+aktscenario+'/'+aktlayout)
os.chdir(aktscenario+'/'+aktlayout)
convertUnix2DosText(files2convert)
os.chdir('..')
files2convert = os.listdir(os.getcwd())
convertUnix2DosText(files2convert)
b = time() - a
print '\n\t', b,' ---> seconds in total for this conflict analysis !!!'
print '\n ------------------------------------------------------------------'
print '  KONFLIKTANALYSE BEENDET! Schliessen Sie dieses Fenster, um '
print '  zur MegasiteManagementToolsuite zurueckzukehren'
print ' ------------------------------------------------------------------'