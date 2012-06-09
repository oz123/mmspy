# -*- coding: utf-8 -*-  

#        This file is part of the conflictAnalysis tool developed in DSITE group
#        Center for Applied Geosciences, University of Tuebingen, Germany
#        The algorithm was created by Max Morio <max.morio__AT__uni-tuebingen.de , 
#        and written in python by Max Morio and Oz Nahum <nahumoz__AT_gmail.com>
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
#  
#       Credits for the function cn_PnPoly():
#       Copyright 2001, softSurfer (www.softsurfer.com)
#       translated to Python by Maciej Kalisiak <mac@dgp.toronto.edu>
#       http://www.dgp.toronto.edu/~mac/e-stuff/point_in_polygon.py


#       this file just holds all the functions used by MMSConflictanalysis.py

from osgeo import ogr
from osgeo import osr
from osgeo import gdal
from osgeo import gdal_array
from numpy import ma, ones,ones_like,  zeros_like, ndenumerate, empty, dtype
import csv
import os, ConfigParser, sys, re
import pdb
config = ConfigParser.RawConfigParser()

def readprojectfile(projectinifile): #(aktscenario, aktlayout, layoutdir, projectinifile):
	'''Read the Project.ini file which determines  
	all the settings and file names required for the MMS program run.
	inputs:
	projectinifile - string - file name path such as on Linux "/home/user/Safira/DATA/Project.ini" or
	on Windows: "C:\Documents and Settings\User\Safira\DATA\Project.ini" 
	outputs:
	aktscenario - string - name of the scenario layout to process, i.e. 'SzenarioA'
	aktlayout - string - name of the layout to process, i.e. 'ScALayout1'
	layoutdir - string - layout direcory, i.e '/home/user/safira/SzenarioA/ScALayout1'
	scen_landuseratio - list - A description of landusage distribution in precentage, 
	an example would be: ['40', '38', '14', '8', '0', '0', '0', '0', '0']
	aoiraster - string - name of Area of Interest raster, i.e 'area_of_interest.asc'
	'''
	config.readfp(open(projectinifile))
	scen_landuseratio = ''    
	aktscenario = config.get('DSS_project', 'AktSzenario')
	aktlayout = config.get('DSS_project', 'AktLayout')
	aoiraster = config.get('DSS_project', 'pathstandort')
	# MAX, 29.04.2010
	# start to add new variable in order to make the code generic 
	# regarding land uses amount and names:
	n_Landuses= int(config.get('DSS_project', 'n_Landuses'))
	layoutdir = os.getcwd() + '/' + aktscenario + '/' + aktlayout
	scen_landuseratio = []
	LUTName = []
	LUTcolours=[]
	for i in range(1, n_Landuses+1):
		scen_landuseratio.append(config.get(aktscenario, 'scen_landuseratio(' + str(i) + ')')) 
		LUTName.append(config.get('DSS_project','LUTName(' + str(i) + ')'))
		LUTcolours.append([config.get('DSS_project','LUTcoloursR(' + str(i) + ')'), \
				config.get('DSS_project','LUTcoloursG(' + str(i) + ')'), \
				config.get('DSS_project','LUTcoloursB(' + str(i) + ')')])
	# get the selected contaminants amongst the available ones ( for which a conflict analysis is to be done )
	no_contaminants = int( config.get( aktscenario, 'anzahlschadstoffe' ))
	selcont=[]
	for i in range( 1, no_contaminants ):
		selcont.append(config.get(aktscenario,'selcont(' + str(i) + ')') )
	return aktscenario, aktlayout, layoutdir , scen_landuseratio, aoiraster, LUTcolours, LUTName, n_Landuses, no_contaminants, selcont
#these variables are not passed to anywhere except aktscenario ?
	#They are passed into one variable which is then accessed as tupple

###########################################################################################
#  Following is executed when "from tools.py import *" is included in python files elsewhere
###########################################################################################

# Define the path to the ini-file 'Projekt.ini
projectinifile = os.getcwd() + '/DATA/Projekt.ini'

# This function parses the Projekt.ini and returns... an inituple
inituple = readprojectfile(projectinifile) 
#pdb.set_trace()
# split the inituple with the returned values
aktscenario = inituple[0]
aktlayout = inituple[1]
#path to the current layout of the active scenario
LAYOUTDIR = inituple[2]
aoiraster = inituple[4]
LUTcolours=inituple[5]
LUTName=inituple[6]
n_Landuses=inituple[7]
no_contaminants = inituple[ 8 ]
selcont = inituple[ 9 ]
###########################################################################################

def WriteRaster(dst_filename, raster, extent): #source shape, destination file, src_filename
	'''This is a generic GDAL function to write ASCII Rasters.
	Here it is an aid function called by ClipPollution and others.
	inputs:
	dst_filename - string - relative path to write the ASCII raster
	raster - ndarray, numpy array - a data array to be written in to the file. 
	extent - list - a list of numbers containing the following data:
	 extent[0] /* top left x */
	 extent[1] /* w-e pixel resolution */
	 extent[2] /* rotation, 0 if image is "north up" */
	 extent[3] /* top left y */
	 extent[4] /* rotation, 0 if image is "north up" */
	 extent[5] /* n-s pixel resolution */
	 for more info on extent see http://www.gdal.org/gdal_tutorial.html
	
	output:
	file - This function writes an ASCII raster file to the same path as 
	dst_filename.
	'''
	format = "MEM"
	driver = gdal.GetDriverByName(format)
	dst_ds = driver.Create(dst_filename, len(raster[0]), len(raster), \
			1, gdal.GDT_Float32)
	#dst_ds.SetGeoTransform( [3366306.98802, 1.0, 0.0, 5814845.04404, 0.0, -1.0] )
	dst_ds.SetGeoTransform(extent)
	#    adfGeoTransform[0] /* top left x */
	#adfGeoTransform[1] /* w-e pixel resolution */
	#adfGeoTransform[2] /* rotation, 0 if image is "north up" */
	#adfGeoTransform[3] /* top left y */
	#adfGeoTransform[4] /* rotation, 0 if image is "north up" */
	#adfGeoTransform[5] /* n-s pixel resolution */
	# Tiff origin is at the left upper corner
	# Ascii raster origin is at the left bottom corner
	#3366306.98802 - x,5814845.04404 - y
	#SetGeoTransform () In a north up image, 
	#[1] is the pixel width, and [5] is the pixel height. 
	#The upper left corner of the upper left pixel is at position [0],[3]).
	#dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
	#pdb.set_trace()
	srs = osr.SpatialReference()
	#srs.ImportFromEPSG(4326) #WGS84 lat long.
	srs.ImportFromEPSG(31468) #Gauss Kruger DHDN zone 4
	#see http://spatialreference.org/ref/epsg/4326/
	#    http://spatialreference.org/ref/epsg/31468/
	dst_ds.SetProjection(srs.ExportToWkt())
	dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
	dst_ds.GetRasterBand(1).WriteArray(raster)
	format = 'AAIGrid'
	driver = gdal.GetDriverByName(format)
	dst_ds_new = driver.CreateCopy(dst_filename, dst_ds)
	dst_ds = None
	
def cn_PnPoly(P, V):
	#modified to works with numpy arrays arrays too
	'''This Function Calculates if a point is in Polygon.
	
	inputs:
	P - list or numpy array - the indecies of a point in the form [x1,y1]. 
	V - list or numpy array - A list containing the all vertices of a polygon in the 
	form [[x1,y1],[x2,y2],[x3,y3],...,[xn,yn]]
	
	output:
	cn - boolean integeer - The function returns 0 if P is inside the
	polygon, and 1 if it outside the polygon.
	
	Copyright 2001, softSurfer (www.softsurfer.com)
	This code may be freely used and modified for any purpose
	providing that this copyright notice is included with it.
	SoftSurfer makes no warranty for this code, and cannot be held
	liable for any real or imagined damage resulting from its use.
	 Users of this code must verify correctness for their application.
	translated to Python by Maciej Kalisiak <mac@dgp.toronto.edu>'''

	cn = 0    # the crossing number counter
	# repeat the first vertex at end
	V = tuple(V[:]) + (V[0],)
	# loop through all edges of the polygon
	for i in range(len(V) - 1):   # edge from V[i] to V[i+1]
		if ((V[i][1] <= P[1] and V[i + 1][1] > P[1])   # an upward crossing
				or (V[i][1] > P[1] and V[i + 1][1] <= P[1])):  # a downward crossing
			# compute the actual edge-ray intersect x-coordinate
			vt = (P[1] - V[i][1]) / float(V[i + 1][1] - V[i][1])
			if P[0] < V[i][0] + vt * (V[i + 1][0] - V[i][0]): # P[0] < intersect
				cn += 1  # a valid crossing of y=P[1] right of P[0]

	#return cn % 2   # 0 if even (out), and 1 if odd (in)
	cn = (cn % 2 - 1) ** 2 # This changes the return to: 0 if in, 1 if out
	return cn 

def readZielwert(aktlayout):
	'''This function reads the Zielwert file, it returns lists which describe 
	allowed threshold for the contaminations, as well as the names of the contminations rasters
	inputs:
	Aktlayout - path to Scenario/Layout which is to be calculated

	outputs:
	contname - list - This list is the 1st column in the Zielwertset.rtv file, which describes the 
	  contaminant names. 
	compartment - list - This list is the 2nd column, which describes where 
	  the contaminant is found (Soil (Boden) or Groudwater (GW)). 
	targets_LUT - list - This list is the 2nd to 2+ no_of_landuses_column in the Zielwertset.rtv file, which holds the 
	  accepted target value for the particular land uses. The first element in the list is the allowed value
	  for the first contaminant (the first row, first element in contname).
	contrastername - list - This list is the 3 previous last column in the Zielwertset.rtv file, which holds the 
	  file names of the contaminent rasters found in the DATA directory
	boolConflictType - list - THis list is the 2nd previous last column, which holds the type of conflict analysis
      ( deterministic 0 or probabilistic 1 )
    iContRealisations - list - This list is the last column, which hold the number of realisations for a 
	   contaminant distribution, or the realisation nr. for deterministic analysis, or nothing.
	no_contaminants - integer - the number of contaminents to be evaluated.
	'''
	global LAYOUTDIR
	infile_target = LAYOUTDIR.replace(aktlayout, '') + 'Zielwertset.rtv'
	
	targetreader = csv.reader(open(infile_target, "rb"))
	try:
		i = 0 
		contname = []
		compartment = []
		# Max, 30.04.2010: replace use-specific targets by generic targets_LUT(n_Landuses)
		targets_LUT=[[]]
		contrastername = []
		# Max. 09.Dec.2010: Added two new names for probabilistic conflict analysis:
		boolConflictType = []
		iContRealisations = [  ]

		for row in targetreader:
			#print 'Row: ',row, 'No.elem.p/row:', len(row), 'last element:',row[-1]
			if i == 0: 
				no_contaminants = int(row[0])
				print 'Anzahl Kontaminanten (B und GW): ', no_contaminants
			else:
				contname.append(row[0])
				compartment.append(row[1])
				# Max, 30.04.2010: read targets dependent on given n_Landuses
				targets_LUT.append([])
				for j in range(2, n_Landuses+2):
					#print 'Getting ', contname[i-1], ' targets for use type index: ', \
					#		(j-1)*100, ' : ',row[j]
					targets_LUT[i].append(row[j])
				contrastername.append(row[n_Landuses+2])
				boolConflictType.append( row[ n_Landuses +3 ] )
				iContRealisations.append( row[ n_Landuses +4 ] )
			i = i + 1
		i = 0
		#while i < len(contname):
		#	print i, contname[i], compartment[i], targets_LUT, contrastername[i] 
		#	#targ_comm[i], targ_resid[i], targ_recreat[i] , targ_agric[i], contrastername[i]
		#	i = i + 1
	except csv.Error, e:
		sys.exit('file %s, line %d: %s' % (infile_target, targetreader.line_num, e))
	#print contname, compartment, contrastername ,no_contaminants, targets_LUT, boolConflictType, iContRealisations 
	return contname, compartment, contrastername ,no_contaminants , targets_LUT, boolConflictType, iContRealisations 

def read_landuseshape(aktscenario,aktlayout, no_contaminants, contname, \
		compartment,targets_LUT, contrastername):
	'''Create shape file output of the layout shapefile which includes the 
	contaminant specific target values.
	inputs:
	aktscenario - string - name of the scenario layout to process, i.e. 'SzenarioA'
	aktlayout - string - name of the layout to process, i.e. 'ScALayout1'
	no_contaminants - integer - the number of contaminents to be evaluated
	contname - list - This list is the 2th column in the Zielwertset.rtv file, which describes the 
	  contaminant names. 
	compartment - list - This list is the 2th column in the Zielwertset.rtv file, which describes where 
	  the contaminant is found (Soil (Boden) or Groudwater (GW)). 
	targets_LUT - list - This list is the 3th column in the Zielwertset.rtv file, which holds the 
	  accepted target value for the commercial landuse. The first element in the list is the allowed value
	  for the first contaminant (the first row in the Zielwertset.rtv file, first element in contname).
	contrastername - list - This list is the last column in the Zielwertset.rtv file, which holds the 
	  file names of the contaminant rasters found in the DATA directory
	
	output:
	<layout>_tgl.shp - ESRI Shape file - A polygon shape file with identical polygon 
	of landuses as in aktlayout.shp file. This file has modified data base
	with additional entries of allowed contaminant threshold for each land use
	and each contaminant. E.g. for 3 contaminants, 3 additional columns will be added
	and for each land use row, 3 entries will be added. 
	''' 
	inFn = aktscenario+'/'+aktlayout+'.shp'
	#append *_tgl.shp for output shapefile name
	outFn = aktscenario+'/'+aktlayout+'/'+aktlayout+'_tgl.shp'
	driver = ogr.GetDriverByName('ESRI Shapefile')
	#create datasource for in-/output, get feature definition
	# open the input shapefile and get the layer and featureDefn
	inDS = driver.Open(inFn)
	if inDS is None:
		print 'Could not open ' + inFn
		return False
	inLayer = inDS.GetLayer()
	featureDefn = inLayer.GetLayerDefn()
	numFields = featureDefn.GetFieldCount()
	# delete the output shapefile if it already exists
	if os.path.exists(outFn):
		driver.DeleteDataSource(outFn)

	# create the output shapefile and get the layer
	outDS = driver.CreateDataSource(outFn)
	if outDS is None:
		print 'Could not create ' + outFn
		return False
	outLayer = outDS.CreateLayer('', geom_type=featureDefn.GetGeomType())
	#add the fields to the layer
	for i in range(numFields):
		outLayer.CreateField(featureDefn.GetFieldDefn(i))
	
	for i in range(no_contaminants):
		# Define attribute table
		field_TEXT  = ogr.FieldDefn(contrastername[i][:-4], ogr.OFTString) #[11]
		field_TEXT.SetWidth(20)
		outLayer.CreateField(field_TEXT)
		print contrastername[i][:-4]
	inFeature = inLayer.GetNextFeature()
	#outFeature = ogr.Feature(outLayer.GetLayerDefn())
	while inFeature: #outFeature: #
		# get the input geometry and reproject it
		geom = inFeature.GetGeometryRef()
		#geom.Transform(coordTrans)
		# create the new feature and set the geometry
		#create outFeature and set the geometry identical to inFeature
		outFeature = ogr.Feature(featureDefn)
		outFeature.SetGeometry(geom)
		# set attributes
		for i in range(numFields):
			outFeature.SetField(i, inFeature.GetField(i))
			outLayer.SetFeature(outFeature)
			#### THIS PUT HERE ONLY COPIES THE FILE
			#### WHEN COMMENTED OUT IT ENSURES FEATURES get extra properties
		#temp = outFeature.getAllAttributeNames()
		#print temp
		outLayer.CreateFeature(outFeature)
		outFeature.Destroy()
		inFeature.Destroy()
		# get the next feature
		inFeature = inLayer.GetNextFeature()
	# close the shapefiles
	inDS.Destroy()
	outDS.Destroy()
	#append Fields for the target values
	ds_tgl = driver.Open(outFn,1)
	tgl_layer = ds_tgl.GetLayer()
	tgl_feature = tgl_layer.GetNextFeature()
	lu_index = tgl_feature.GetFieldIndex("Kategorie") #Get the row number of the land use type index field
	#print "Fieldnames: ", tgl_feature.getName()
	# get the field columns for all contaminants and put targets in 
	while tgl_feature:
		# Max, 30.04.2010: changed in order to generalize according to n_Landuses
		lu_indexcheck = tgl_feature.GetField(lu_index)
		for n in range (100, n_Landuses*100, 100):
			if lu_indexcheck==n:
				for i in range(no_contaminants):
					searchforfieldname=contrastername[i][:-4]
					cont_index=tgl_feature.GetFieldIndex(searchforfieldname)
					#print 'indexes i, n: ', i+1, ' , ', n, ' fieldname: ', searchforfieldname, ' targets_LUT: ', targets_LUT[i+1][(n/100)-1]
					tgl_feature.SetField(cont_index,  targets_LUT[i+1][n/100-1])
					tgl_layer.SetFeature(tgl_feature)
			if (lu_indexcheck<100) or (lu_indexcheck>(n_Landuses*100)):
				searchforfieldname=contrastername[i][:-4]
				cont_index=tgl_feature.GetFieldIndex(searchforfieldname)
				tgl_feature.SetField(cont_index, 99999999.99)
				tgl_layer.SetFeature(tgl_feature)
		tgl_feature.Destroy()
		tgl_feature=tgl_layer.GetNextFeature()
	ds_tgl.Destroy()
	#else: print 'Fehler!'     
	print "Layout shapefile mit Zielwerten: ", aktscenario+'/'+aktlayout+'/'\
			+aktlayout+'_tgl.shp'

def CreateTargetsRisk(nameOfLandUsesRaster, contname, compartment, targets_LUT,\
	contrastername, contNumber,aktscenario, aktlayout):
	'''THis function reads land-use specific uniform exceedance rasters and creates
	a layout-specific exceedance file under /DATA/RISK/. named similar to the original
	contaminant raster file in /DATA/.

	Inputs:
	nameOfLandUsesRaster - string - a path to ESRI ASCII raster, i.e 'SzenarioA/ScALayout1.asc'
	contname - list - This list is the 2th column in the Zielwertset.rtv file, which describes the 
	  contaminant names.
	compartment - list - This list is the 2th column in the Zielwertset.rtv file, which describes where 
	  the contaminant is found (Soil (Boden) or Groudwater (GW)). 
	targets_LUT [[]] - list - This list is the 3rd to 3 previous the last column in the Zielwertset.rtv 
	  file, which holds the 
	  accepted target values for the landuse types. The first element in the list is the allowed value
	  for the first contaminant (the first row in the Zielwertset.rtv file, first element in contname).
	contrastername - list - This list is the last column in the Zielwertset.rtv file, which holds the 
	  file names of the contaminent rasters found in the DATA directory
	contNumber - integer - indexing value used for accessing values in the lists above.
	aktscenario - string - name of the scenario layout to process, i.e. 'SzenarioA'.
	aktlayout - string - name of the layout to process, i.e. 'ScALayout1'.
	
	Output:
	ESRI ASCII Rasters - One ascii raster with target values for the actual layout
	'''
	
	nameOfUnifExceedanceRaster=[]
	iRunsPerContaminant=[]

	# TODO open land-use raster 
	dataset = gdal.Open(nameOfLandUsesRaster)
	extent = dataset.GetGeoTransform()
	landuses = gdal_array.DatasetReadAsArray(dataset)
	lay_exceedcontraster = landuses.astype( 'float') #array( zeros_like( landuses ), dtype=float)
	for item in contrastername: # loops over contaminant
		for n in range ( 1, n_Landuses+1 ):
			# open the land use specific uniform exceedance raster for contrastername item
			exceedanceRiskRaster = 'DATA/RISK/'+ str(n) + item.replace( '.aux', '.asc' )
			exceedanceRiskRaster = os.getcwd() +'/'+ exceedanceRiskRaster
			print exceedanceRiskRaster
			exceedanceRiskRaster = exceedanceRiskRaster[:-1]+exceedanceRiskRaster[-1].lower()
			dataset2 = gdal.Open(exceedanceRiskRaster)
			
			exceedancecontraster = gdal_array.DatasetReadAsArray( dataset2 )
			# loop over i,j of the raster land uses and get the i,j from contraster, if == n*100
			# close the land use specific uniform exceedance raster for contrastername item
			for index, x in ndenumerate(landuses):
				#print index, 'index', x, 'x', 'n', n
				if x == int( n*100):
					lay_exceedcontraster[ index[0]][index[1] ] =  exceedancecontraster[ index[0]][index[1] ]
					#print "\n CHECK!!!",  exceedancecontraster[ index[0]][index[1] ]
				elif x == int(-9999):
					lay_exceedcontraster[ index[0]][index[1] ] = -9999
				#else:
				#	lay_exceedcontraster[ index[0]][index[1] ] = -1

		print 'writing raster ', 'DATA/RISK/lay_'+ item.replace( '.aux', '.asc' )
		WriteRaster('DATA/RISK/lay_'+ item.replace( '.aux', '.asc' ), lay_exceedcontraster, extent)
	return


def CreateTargets(nameOfLandUsesRaster, contname, compartment, targets_LUT,\
		contrastername, contNumber,aktscenario, aktlayout):
	'''This function reads inland use raster and creates target value rasters
	for contaminants. Each cell in the raster is check in which polygon it is, and is assigned
	the value of target value according to the land use found in the field "Kategorie" and the 
	matching target value. 
	
	Inputs:
	nameOfLandUsesRaster - string - a path to ESRI ASCII raster, i.e 'SzenarioA/ScALayout1.asc'
	contname - list - This list is the 2th column in the Zielwertset.rtv file, which describes the 
	  contaminant names.
	compartment - list - This list is the 2th column in the Zielwertset.rtv file, which describes where 
	  the contaminant is found (Soil (Boden) or Groudwater (GW)). 
	targets_LUT [[]] - list - This list is the 3rd to 3 previous the last column in the Zielwertset.rtv 
	  file, which holds the 
	  accepted target values for the landuse types. The first element in the list is the allowed value
	  for the first contaminant (the first row in the Zielwertset.rtv file, first element in contname).
	contrastername - list - This list is the last column in the Zielwertset.rtv file, which holds the 
	  file names of the contaminent rasters found in the DATA directory
	contNumber - integer - indexing value used for accessing values in the lists above.
	aktscenario - string - name of the scenario layout to process, i.e. 'SzenarioA'.
	aktlayout - string - name of the layout to process, i.e. 'ScALayout1'.
	
	Output:
	ESRI ASCII Raster - An ascii raster with target values for each land use will be
	created in the "scenario"/"layout" directory. 
	'''
	#print 'Inside Create Targets'
	#print 'current dir', os.getcwd()
	#pdb.set_trace()
	#print aktscenario, 'aktscenario'
	dataset = gdal.Open(nameOfLandUsesRaster)
	extent = dataset.GetGeoTransform()
	landuses = gdal_array.DatasetReadAsArray(dataset)
	#Creat an empty array to hold the values of thresholds
	#this array has the same shape as landuses, and is initialized with ones.
	#thresholds = numpy.ones_like(landuses) * -9999
	thresholds = ones_like(landuses) * -9999
	# codes for land uses are: e.g. 100, 200, ... 
	#print 'Creating target values raster for ', contname[contNumber]
	stor = [] #store unique values
	#for index, x in numpy.ndenumerate(landuses):
	for index, x in ndenumerate(landuses):
		# land usues codes are NOT constant any more !!!
		# Max, 30.04.2010 changed this part i.o. to generalize acc. to n_Landuses^
		for j in range(100,100+ n_Landuses*100, 100):
			#if x==900:
			if x==j:
				#print 'Check: x ', x, ' targLUT contN, j/100-1: ', contNumber+1,' ', (j/100)-1
				thresholds[index[0]][index[1]]=targets_LUT[contNumber+1][(j/100)-1]
				stor.append(x)
		if x != -9999 and x!=j and x not in stor:
			print 'FEHLER: Unbekannte LandnutzungsID', x, ' wird nicht beruecksichtigt'
			stor.append(x)
			#    thresholds[index[0]][index[1]] = x #this line is responsible of transforming all out of areas to -9999
	#pdb.set_trace()
	#WriteRaster (aktscenario+'/'+aktlayout+'/'+ contrastername[contNumber].replace('.aux','')+'_target.asc', thresholds, extent)
	WriteRaster(aktscenario + '/' + aktlayout + '/' + contname[contNumber] + '_target.asc', thresholds, extent)
	#print 'Created Raster of Thresholds in: ', aktscenario+'/'+aktlayout +'/' + contrastername[contNumber].replace('.aux','')+ '_target.asc'
	print 'Zielwertraster abgelegt in: ', aktscenario + '/' + aktlayout + '/' + contname[contNumber] + '_target.asc'

def ClipPollution2(nameOfMask, nameOfPollutionRaster):
	'''This function takes in two rasters, and returns a clipped raster, 
	inputs:
	nameOfMask - string - name pointing to  the mask raster. This is 
	a ESRI ASCII raster,  containing 0 for points INSIDE area of interest.
	For all points out side of the domain of interest 1  values. 
	nameOfPolutionRaster - string - name pointing to raster of pollution. 
	
	outputs:
	ESRI ASCII Raster - An ascii raster with contamination values clipped
	to the extent of interest. 
	'''
	nameOfPollutionRaster = os.getcwd() +'/'+ nameOfPollutionRaster
	dataset = gdal.Open(nameOfPollutionRaster)
	extent = dataset.GetGeoTransform()
	dataSetMask = gdal.Open(nameOfMask)
	areaOfPollution = gdal_array.DatasetReadAsArray(dataset)
	maska = gdal_array.DatasetReadAsArray(dataSetMask)
	#clip = numpy.ma.MaskedArray(areaOfPollution, mask=maska)
	clip = ma.MaskedArray(areaOfPollution, mask=maska)
	cclip=ma.filled(clip,fill_value=-9999)
	#cclip = numpy.ma.filled(clip, fill_value= -9999)
	#WriteRaster(nameOfPollutionRaster+'clipped', cclip, extent)
	nameOfPollutionRaster = nameOfPollutionRaster.replace('.asc', '_clipped.asc')
	WriteRaster(nameOfPollutionRaster, cclip, extent)
	#print 'Successfuly Cliped pollution Raster...', ' clipped ' + nameOfPollutionRaster
	return nameOfPollutionRaster
#scheme to implement:
#Read the masked_area_of interest
#Create an empty array 1 row array in the same size as the size of new array of valid rows and columns.
#for each point, in the valid rows and valid columns, make a test PIP to determine in which polygon it is.
#Then replace the 0 value with the Kategorie value of the polygon. 
#Reshape the vector, and write down the raster. 
#Useful functions of area :  DumpReadable(*args), GetFieldAsInteger(*args)

#extent = layer.GetExtent()
#print 'Extent:', extent
#print 'UL:', extent[0], extent[3]
#print 'LR:', extent[1], extent[2]

#print "Number of Features in layer:", layer.GetFeatureCount()
#From here it's easy to continue... we know the number of polygons in the layer. 
#So we iterate over each polygon and get it's extent. The we also use it's function GetFieldAsInteger(*args)
#Now we iterate over each node and determine in which polygon it is. When we find the poligon that it's in, we 
#we write in the node location the number from GetFieldAsInteger, thus rasterizing the Land Uses 
#from the scenarion is complete. 
#layer.GetFeature()

### READ THE SCENARIO LAND USES, MAINLY, PRODUCE AN ARRAY WITH THE BOUNDARIES OF LAND USES POLYGONS
def getLayerProperties(layer):
	'''This function gets a layer object from an ESRI Shape file and 
	reads its following properties:
	LandUses, Codes, Polygon Boundaries
	input:
	layer - osgeo.ogr.Layer - an ESRI Shape file data object, which can be read
	with the following set of commands:
	driver = ogr.GetDriverByName('ESRI Shapefile')
	dataSource = driver.Open('Scenario.shp', 0)
	layer = dataSource.GetLayer()
	
	output:
	LandUses - list of strings - A list containing the land names code of each polygon 
	in the layer.
	Codes - list of integers - A list containing the land use code of each polygon 
	in the layer.
	Boundaries - list of lists - A list containing a all the polygons in
	the layer. Each polygon is represented by a list of polygon vertices
	in the form [[x1,y1],[x2,y2],[x3,y3],...,[xn,yn]].
	'''
	LandUses = []
	Codes = [] #Store Categories
	Boundaries = [] # store boundaryies

	for i in range(layer.GetFeatureCount()):
		area = layer.GetFeature(i)
		Name = area.GetField(0)
		Kategorie = area.GetField(1)
		LandUses.append(Name)
		Codes.append(Kategorie)
		geometry = area.GetGeometryRef()
		boundary_raw = str(geometry.GetBoundary())
		#remove tail and head
		boundary = boundary_raw[12:-1]
		boundary = list(boundary.split(','))
	#print boundary
		lenboundary = len(boundary)
		for i in range(lenboundary):
			#convert the string to a list of lists (couples of x y coordinates)
			#boundary[i]=boundary[i].replace(' ',',')
			boundary[i] = list(boundary[i].split(' '))
	#convert each coordinate from string to float
		for i in range(lenboundary):
			boundary[i][0] = float(boundary[i][0])
			boundary[i][1] = float(boundary[i][1])
		Boundaries.append(boundary)
	return LandUses, Codes, Boundaries

#def readShpAsData(dataSource):
	#print (dataSource)
	#driver = ogr.GetDriverByName('ESRI Shapefile')
	##print type(driver)
	#dataSource = driver.Open(dataSource)
	####print 'asf'
	#layer = dataSource.GetLayer()
	#print type(layer)
	#return layer
#put these lines in main script
#LandUses, Codes, Boundaries = getLayerProperties(layer)
##print Boundaries
#print Codes
#print LandUses

#before this there is a need to have already the ascii raster of the pollution or aoi. 
#dataSource = 'DATA/area_of_interest_raster_numpy.asc'
##driver = ogr.GetDriverByName('ESRI Shapefile')
#dataSource = gdal.Open(dataSource, 0)
#dataArray = gdal_array.DatasetReadAsArray(dataSource)
#contArray=dataArray

def createLandUsesRaster(layerBoundaries, LayerCodes, LandUses, aktscenario, aktlayout, nameOfMask):
	'''This function creates a raster with the landuses codes, with the
	extent of the polution rasters.
	This function uses the functions: cn_PnPoly and WriteRaster
	inputs:
	layerBoundaries - list of lists - A list containing a all the polygons in
	  the layer. Each polygon is represented by a list of polygon vertices
	  in the form [[x1,y1],[x2,y2],[x3,y3],...,[xn,yn]]
	LayerCodes - list of integers - A list containing the land use codes of each polygon 
	  in the layer.
	LandUses - list of strings - A list containing the land names code of each polygon 
	  in the layer.
	aktscenario - string - name of the scenario layout directory to process, i.e. 'SzenarioA'.
	aktlayout - string - name of the layout to process, i.e. 'ScALayout1'.
	nameOfMask - string - name pointing to  the mask raster. This is an ESRI ASCII 
	  raster,  containing 0 for points INSIDE area of interest.
	For all points out side of the domain of interest 1  values. 
	
	outputs:
	ESRI ASCII file - saved in aktscenario + '/' + aktlayout + '.asc'. 
	Each cell contains the land use according to the polygon landuse Shape file
	'aktlayout'.shp 
	nameOfLandUsesRaster - string - pointing to the above created ESRI ASCII raster.
	'''
	#points=[] #possible to make code faster if points is array and values are assigend
	#points2=empty_like(dataArray) #create empty array in the shape of original raster (see numpy.empty_like)
	
	#dataSource = gdal.Open(dataSource)
	dataSource = gdal.Open(nameOfMask)
	
	extent = dataSource.GetGeoTransform()
	width = dataSource.RasterXSize
	height = dataSource.RasterYSize
	xOrigin, yOrigin = extent[0], extent[3]
	print 'origin of raster (x,y):', xOrigin, yOrigin
	print 'dimensions of pollution raster (w, h)', width, height
	print 'extent of pollution raster', extent
	#read mask
	nameOfMask = gdal.Open(nameOfMask)
	mask = gdal_array.DatasetReadAsArray(nameOfMask)
	##find first point which contain the data
	#rowsholder=[]
	#for i in range(mask.shape[0]):
		#if all(mask[i]) != True:
			#rowsholder.append(i)
	##keep only first and last index
	#del rowsholder[1:-1]
	#print rowsholder, 'rowsholder'
	#columnsholder = []
	#for i in range(mask.shape[1]):
		#if all(mask[:,i]) != True:
			#columnsholder.append(i)
	#del columnsholder[1:-1]
	#print columnsholder, 'columnsholder'
	#create empty array with rows as number of 
	#cells in original raster (see numpy.empty), and two culomns
	#points2 = numpy.empty([height * width, 2]) 
	points2 = empty([height * width, 2]) 
	pointsholder = [] #hold all points in mask which are not masked
	c = 0 #just a counter
	maskf = mask.flat
	for item in maskf:
		if item == 0:
			pointsholder.append(c)
		c = c + 1
	#print pointsholder
	cell = 0
	for j in range(height):
		for i in range(width):
			#This bellow line is modified to fit the general case, not just cell size 10
			points2[cell] = [xOrigin + extent[1] / 2 + i * extent[1], yOrigin - extent[1] / 2 - extent[1] * j] 
			cell = cell + 1
	#for every point in the dataSource, check in which polygon it is. If in none, do not modify it...
	#The more efficient way to do it would be to start from the inner places
	landUsesArray = ones((width * height)) * (-9999) #initiate all with nulls
	#print 'type of landUsesArray is ', type(landUsesArray)
	# check wheter a point 
	# here in this loop we are cross checking all polygon vs. all points
	# first scheme : check all boundaries and points
	#for b in range(len(layerBoundaries)):
	#    print 'Rasterizing ', LandUses[b]
	#    for i in range(pointsholder[0],pointsholder[-1]):
	#        if cn_PnPoly(points2[i], layerBoundaries[b])==False:
	#            landUsesArray[i]=LayerCodes[b]
	####second scheme: check point until find boundary and then quit and go to next point
	for point in range(pointsholder[0], pointsholder[-1]):
		for landUse in layerBoundaries:
			if cn_PnPoly(points2[point], landUse) == False:
				landUsesArray[point] = LayerCodes[layerBoundaries.index(landUse)]
				break
	#lun=len(layerBoundaries) #landuses number
	###third scheme point store last boundary and start check from there
	#point = pointsholder[0]
	#endpoint = pointsholder[-1]
	#while point < endpoint:
		#pdb.set_trace()
		#for i in range(lun):
			#loc = cn_PnPoly(points2[point],layerBoundaries[i])
			#if loc == False:
				#luse=LayerCodes[i]
				#landUsesArray[point]=luse
				#boundary=layerBoundaries[i]
				#break
		#point = point+1
		#while cn_PnPoly(points2[point],boundary) == False:
			#landUsesArray[point]=luse
			#point=point+1
	#print 'type of landUsesArray is ', type(landUsesArray)
	landUsesArray = landUsesArray.reshape(height, width)
	#print landUsesArray.shape
	#pdb.set_trace()
	nameOfLandUsesRaster = 'DATA/' + aktlayout + '.asc'
	WriteRaster(aktscenario + '/' + aktlayout + '.asc', landUsesArray, extent)
	#print 'Raster erzeugt... in: ', aktscenario + '/' + aktlayout + '.asc'
	return nameOfLandUsesRaster

### functions originaly in getMask
def getMask(AreaShapeLayer, ContRaster):
	'''This function creates a raster with the extent of the pollution raster. This raster
	has values of 1 if the cell is outside the Area of Interest and 0 if inside. This raster is 
	served as a mask for the masked arrays algebra (see terminology of masked arrays for more details
			input:
			AreaShapeLayer - string - name of shape file containing a polygon bounding the area of interest
			ContRaster - name - A name of contaminant raster
			
			output:
			MaskRaster - file, ESRI ASCII format - The returned mask raster of 0's inside the area of interest,
	1's in all the areas out of the interest.)'''
	dataSource = 'DATA/area_of_interest.shp'
	driver = ogr.GetDriverByName('ESRI Shapefile')
	dataSource = driver.Open(dataSource, 0)
	layer = dataSource.GetLayer()
	extent = layer.GetExtent()
	print 'Extent des Bewertungsgebiets (area_of_interest):', extent
	print 'UL:', extent[0], extent[3]
	print 'LR:', extent[1], extent[2]
	area = layer.GetFeature(0)

	#Useful functions of area :  DumpReadable(*args), GetFieldAsInteger(*args)
	geometry = area.GetGeometryRef()
	boundary_raw = str(geometry.GetBoundary())
	#pdb.set_trace()
	#need some processing to convert boundary raw to a usefull form !
	#remove tail and head
	boundary = boundary_raw[12:-1]
	boundary = list(boundary.split(','))#convert the string to a list of strings
	#print boundary
	lenboundary = len(boundary)
	for i in range(lenboundary): #convert the string to a list of lists (couples of x y coordinates)
		boundary[i] = list(boundary[i].split(' '))
	#print boundary
	#convert each coordinate from string to float
	for i in range(lenboundary):
		boundary[i][0] = float(boundary[i][0])
		boundary[i][1] = float(boundary[i][1])
	#done with preprocessing, now boudary is a list of lists containing processable floats.
	dataSource = ContRaster
	#pdb.set_trace()
	dataSource = gdal.Open(dataSource)
	dataArray = gdal_array.DatasetReadAsArray(dataSource)
	extent = dataSource.GetGeoTransform()
	width = dataSource.RasterXSize
	height = dataSource.RasterYSize
	xOrigin, yOrigin = extent[0], extent[3]
	print 'origin of raster (x,y):', xOrigin, yOrigin
	print 'dimensions of pollution raster (w, h)', width, height
	print 'extent of pollution raster', extent
	##check all nodes - if out than node values 1,     if inside then no mask, node value 0. 
	#points=[] #possible to make code faster if points is array and values are assigend
	#points2=empty_like(dataArray) #create empty array in the shape of original raster (see numpy.empty_like)
	#points2 = numpy.empty([height * width, 2]) #create empty array with rows as number of cells in original raster (see numpy.empty), and two culomns
	points2 = empty([height * width, 2]) #create empty array with rows as number of cells in original raster (see numpy.empty), and two culomns
	###Create the Points array with numpy array instead of list
	cell = 0
	for j in range(height):
		for i in range(width):
			#This below line is modified to fit the general case, not just cell size 10
			points2[cell] = [xOrigin + extent[1] / 2 + i * extent[1], yOrigin - extent[1] / 2 - extent[1] * j] 
			cell = cell + 1
	#mask = numpy.ones(width * height, dtype=int)
	mask = ones(width * height, dtype=int)
	for i in range(len(points2)):
		re = cn_PnPoly(points2[i], boundary)
		mask[i] = re
	mask = mask.reshape(height, width) #rows, columns
	### extract rows with significant data, i.e erase all rows which are only nulls
	### We extract rows from data array - the pollution raster - according to the temp mask we created
	#pdb.set_trace()
	rowsholder = []
	for i in range(mask.shape[0]):
		if all(mask[i]) != True:
			rowsholder.append(i)
	### Extract columns
	columnholder = []
	for i in range(mask.shape[1]):
		if all(mask[:, i]) != True:
			columnholder.append(i)
	# it's only important to know indicies of rows + cols which contain significant data, no
	# real need to go crazy and resize the array
	nameOfMask = 'DATA/' + AreaShapeLayer.replace('shp', 'asc')
	print 'Writing mask raster...', nameOfMask
	WriteRaster(nameOfMask, mask, extent) #destination file, dataArray,extent
	return nameOfMask
### end of getMask 

def CalculateExceedances(nameOfThreasholds, nameOfMask, nameOfClippedRaster, nameOfContaminant, \
		nameOfCompartment, aktscenario, nameOfPollutionRaster):
	'''
	This function reads three rasters and returns one raster.
	It takes the pollution raster, area of interest raster and allowed thresholds raster.
	The area of interest raster is used as a mask for the two others, and the output 
	raster is the result of dividing the contamination raster by the allowed thresholds
	raster. 
	inputs:
	nameOfThreasholds - string - a name pointing to target values raster, 
	  e.g. 'aktscenario/aktlayout/PAK_B_target.asc'.
	nameOfMask - string - name pointing to  the mask raster. This is an ESRI ASCII 
	  raster,  containing 0 for points INSIDE area of interest.
	nameOfClippedRaster - string - path to clipped contaminant raster created by ClipPollution2(), 
	  e.g. 'DATA/TCE_in_GW_clipped.asc'
	nameOfContaminant - string - name of the contaminant as read from Zielwertset.rtv, 
	  e.g. 'PCE'.
	nameOfCompartment - string - where the contaminant is found, e.g 'B' or 'GW'.
	aktscenario - string - name of the scenario layout directory to process, i.e. 'SzenarioA'.
	nameOfPollutionRaster - string - name of the original pollution raster as read
	  from Zielwertset.rtv, e.g. 'PCE_in_GW.aux'.
	nameOfPolutionRaster - string - name pointing to raster of pollution
	
	outputs:
	Exceedance raster to be named 'component'_(contaminant_rastername).asc, 
	i.e B_PAK_B.asc
	Exceedance raster with Boolean values (0 for no exceedance and 1 for exceedance) is 
	saved under TEMP
	'''
	#first step: read raster of pollution
	#2. read raster of targets
	#3. devide pollution by target - result is exceedence raster
	#4. read mask
	#5. clip mask and exceedence
	#6. write file !!!
	# THIS is much more fast than the previous scheme (CreateExidences function) 
	#where each cell was checked !!!
	dataset = gdal.Open(nameOfClippedRaster)
	extent = dataset.GetGeoTransform()
	polutions = gdal_array.DatasetReadAsArray(dataset)
	dataset = gdal.Open(nameOfThreasholds)
	thresholds = gdal_array.DatasetReadAsArray(dataset)
	exceedance = polutions / thresholds
	### Write 0 for False, write 1 for True
	#exidenceOnesZeros = numpy.zeros_like(exidence)
	exceedanceOnesZeros = zeros_like(exceedance)
	#for index, x in numpy.ndenumerate(exidence):
	for index, x in ndenumerate(exceedance):
		if x < 0:
			exceedanceOnesZeros[index[0]][index[1]] = -9999
		elif x < 1:
			exceedanceOnesZeros[index[0]][index[1]] = 0
		elif x > 1:
			exceedanceOnesZeros[index[0]][index[1]] = 1
	#clip unwated 1's if creating exidence raster with excidance values
	dataSetMask = gdal.Open(nameOfMask)
	maska = gdal_array.DatasetReadAsArray(dataSetMask)
	masked_exceedance = ma.MaskedArray(exceedance, mask=maska)
	masked_exceedance = ma.masked_inside(masked_exceedance, -9999, -0.0000001)
	exceedance = ma.filled(masked_exceedance, fill_value= -9999)
	#pdb.set_trace()
	nameOfPollutionRaster = nameOfPollutionRaster.replace('.aux', '')
	###Create Exidance Rasters
	if nameOfCompartment == 'Boden':
		exceedance_raster = aktscenario + '/' + aktlayout + '/' + 'B_' \
				+ nameOfPollutionRaster + '.asc'
	else:
		exceedance_raster = aktscenario + '/' + aktlayout + '/' + 'W_' \
				+ nameOfPollutionRaster + '.asc'
	WriteRaster(exceedance_raster , exceedance, extent)
	#print  'Raster Zielwertueberschreitung(-sfaktor) erzeugt: ', exceedance_raster 
	#pdb.set_trace()
	if os.path.exists(aktscenario + '/' + aktlayout + '/' + 'TEMP/') == False:
		os.mkdir(aktscenario + '/' + aktlayout + '/' + 'TEMP/')
	###Create Excidance Raster of 0's and 1's
	masked_exceedanceOnesZeros = ma.MaskedArray(exceedanceOnesZeros, mask=maska)
	exceedanceOnesZeros = ma.filled(masked_exceedanceOnesZeros, fill_value= -9999)
	exceedance_raster01s = exceedance_raster.replace(aktlayout, aktlayout + '/TEMP')
	exceedance_raster01s = exceedance_raster01s.replace(nameOfContaminant, \
			nameOfContaminant + '01s')
	WriteRaster(exceedance_raster01s , exceedanceOnesZeros, extent)
	print  'True/false Ueberschreitungskarte in: ', exceedance_raster01s
	return exceedance_raster, exceedance_raster01s
	
def CreateTempExceedPolygon(src_filename, dst_filename , dst_layername,
		dst_fieldname='None', src_band_n=1, format='ESRI Shapefile'):
	'''Polygonize the bolean exceedance rasters, which contain 0's for 
	no exceedance and 1's for exceedeance. 
	The polygons are created according to the same classification. Polygons
	containing 1's are for regions with exceedance, polygons with 0's are
	for regions with no exceedance. 
	These exceedance polygons are later intersected with the land used 
	polygons to create the exceedance for each landuse, see function CreatContExceedance.
	This function is coppied from gdal_polygonize utility.
	inputs:
	src_filename - string - path to a ESRI ASCII raster. 
	dst_filename - string - path to ESRI Shape file created as directory 
	which contains all the files of the shape file. 
	dst_layername - string - A string which describes the layer name
	to be created. Should be something informative if left empty the function
	creates a layer named 'out'.
	
	outputs:
	ESRI Shape file - shape file name, created as directory 
	which contains all the files of the shape file. 
	'''
	# =============================================================================
	#    Open source file
	# =============================================================================
	#format = 'ESRI Shapefile'
	options = []
	quiet_flag = 0
	dst_field = -1
	mask = 'default'
	src_ds = gdal.Open(src_filename)
	if src_ds is None:
		print 'Unable to open ', src_filename
		sys.exit(1)
	srcband = src_ds.GetRasterBand(src_band_n)
	if mask is 'default':
		maskband = srcband.GetMaskBand()
	elif mask is 'none':
		maskband = None
	else:
		mask_ds = gdal.Open(mask)
		maskband = mask_ds.GetRasterBand(1)
	# =============================================================================
	#       Try opening the destination file as an existing file.
	# =============================================================================
	try:
		gdal.PushErrorHandler('QuietErrorHandler')
		dst_ds = ogr.Open(dst_filename, update=1)
		gdal.PopErrorHandler()
	except:
		dst_ds = None
	# =============================================================================
	#     Create output file.
	# =============================================================================
	if dst_ds is None:
	# This means that if there's no file you'll see the output. If file exists
	# only the precentage out put will be printed  
	# Creating output HUMUS of format ESRI Shapefile.
		drv = ogr.GetDriverByName(format)
		#print 'Creating output %s of format %s.' % (dst_filename, format)
		dst_ds = drv.CreateDataSource(dst_filename)
	# =============================================================================
	#       Find or create destination layer.
	# =============================================================================
	try:
		dst_layer = dst_ds.GetLayerByName(dst_layername)
	except:
		dst_layer = None

	if dst_layer is None:
		srs = None
		if src_ds.GetProjectionRef() != '':
			srs = osr.SpatialReference()
			srs.ImportFromWkt(src_ds.GetProjectionRef())
		dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
		if dst_fieldname is None:
			dst_fieldname = 'Excidance'
		fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
		dst_layer.CreateField(fd)
		dst_field = 0
	# =============================================================================
	#    Invoke algorithm.
	# =============================================================================
	prog_func = gdal.TermProgress
	result = gdal.Polygonize(srcband, maskband, dst_layer, dst_field, options,
			callback=prog_func)
	srcband = None
	dst_ds = None

def CreatContExceedance(tempShp, aktscenario, aktlayout, nameOfContaminant, nameOfCompartment):
	'''This function creates the final polygon excidance layer. The polygons in this layers are
	only of the exceedence zone interstected with th landuse polygons. The attribute table
	of the Shape File contains three columns: 
	id - serial identifier of the polygons.
	land use c - landuse code of the polygon, according to field 'Kategorie' which was 
	intersected in the layout layer.
	area - the area of the polygon in sq. m. 
	input:
	tempShp - string - path to ESRI Shape file with the exceedence polygons created by CreateTempExceedPolygon()
	aktscenario - string - name of the scenario layout to process, i.e. 'SzenarioA'.
	aktlayout - string - name of the layout to process, i.e. 'ScALayout1'.
	nameOfContaminant - string - name of the contaminant, i.e. 'TCE'
	nameOfCompartment - string - name of the compartment, e.g. 'B' or 'GW'.
	
	output:
	ESRI SHAPE file containing polygons and attribute table as described above.
	 the file will be named according to component name and the contaminant name,
	 e.g GW_TCE or B_PAK inside a direcotry with the same name.
	dst_filename - string - name of the above created file.
	'''
	#pdb.set_trace()
	tempShpString = tempShp
	driver = ogr.GetDriverByName('ESRI Shapefile')
	tempShp = driver.Open(tempShp)
	PAK_B = tempShp.GetLayer()
	PAKB_spatialRef = PAK_B.GetSpatialRef()
	b = PAKB_spatialRef.GetLinearUnitsName()
	#print "Length Units of Layer", tempShpString, b
	dataSource = '../..' + '/' + aktlayout + '.shp' #'SzenarioA/ScALayout1.shp'
	#print "Working on Layout...", aktlayout
	dataSource = driver.Open(dataSource, 0)
	LandUses = dataSource.GetLayer()
	LU_spatialRef = PAK_B.GetSpatialRef()
	b = LU_spatialRef.GetLinearUnitsName()
	LU_SRS = LU_spatialRef.ExportToWkt()
	#print '######################################'
	#print "Number of features in landuses layer... ", LandUses.GetFeatureCount()
	#print "Number of features in contaminant layer.. ", PAK_B.GetFeatureCount()
	#print 'temp shape file: ...', tempShpString
	#os.chdir('DATA')
	#os.chdir(aktscenario+'/'+aktlayout)
	os.chdir('..')
	dst_filename = tempShpString[5:]#+'.shp' #CHANGE THIS TO APPROPRIATE POLLUTANT NAME
	if nameOfCompartment == 'Boden':
		dst_filename = 'B_' + nameOfContaminant#+'.asc'
	else:
		dst_filename = 'W_' + nameOfContaminant#+'.asc'
	#dst_filename=nameOfContaminant
	print 'Kompartiment_Schadstoff: ', dst_filename
	format = 'ESRI Shapefile'
	#LandUses = dataSource.GetLayer()
	src_ds = LandUses
	#import project of the land uses layer srs=Spatial Reference System
	if LU_SRS != '':
		srs = osr.SpatialReference()
		srs.ImportFromWkt(LU_SRS)
	#print 'srs ',srs
	# Create Shapefile named "‏‎YOSEF.shp" for output
	driver = ogr.GetDriverByName(format) #[2] 
	destinationSource = driver.CreateDataSource(dst_filename) #[4]
	#create polygon layer
	dst_layer = destinationSource.CreateLayer(dst_filename, geom_type=ogr.wkbPolygon, srs=srs) #[5]
	# Define attribute table
	field_ID = ogr.FieldDefn()    #[6]
	field_ID.SetName('ID')           #[7]
	field_ID.SetType(ogr.OFTInteger) #[8]
	field_ID.SetWidth(5)             #[9] #damn this width !!!
	dst_layer.CreateField(field_ID)  #[10]
	#field_ID = 0
	#add field named moshe
	field_LUC = ogr.FieldDefn()    #[6]
	#field_LUC.SetName('LAND USE CODE')
	field_LUC.SetName('CATEGORIE')  #Max length of field name is 10 char.
	#see http://en.wikipedia.org/wiki/Shapefile        
	field_LUC.SetType(ogr.OFTInteger) #[8]
	field_LUC.SetWidth(15)             #[9] #damn this width !!!
	dst_layer.CreateField(field_LUC)  #[10]
	
	#dst_fieldname = 'LAND USE CODE'
	#pdb.set_trace()   
	#fd = ogr.FieldDefn(dst_fieldname, ogr.OFTInteger)
	#fd.SetName('LAND USE CODE')
	#fd.SetWidth(15)
	#dst_layer.CreateField(fd)
	#dst_field = 0
	#add field named "area"
	dst_fieldname = 'AREA'
	fd = ogr.FieldDefn(dst_fieldname, ogr.OFTReal)
	dst_layer.CreateField(fd)
	#dst_field = 0
	# Having done all predefinitions, it's time to create features. For this
	# purpose instantiate a new feature.
	ID = 0
	feature = ogr.Feature(dst_layer.GetLayerDefn()) #[12]
	#print 'Iterating over input data ..........'
	# Iterate over input data
	#count = len(data)
	#pdb.set_trace()
	for i in range(LandUses.GetFeatureCount()):
		#ID=ID+1
		Whon = LandUses.GetNextFeature()
		for j in range(PAK_B.GetFeatureCount()):
			#print i,j
			Excidance = PAK_B.GetFeature(j)
			ExcidanceYN = Excidance.GetField(0)
			#print Ex
			if ExcidanceYN == 1:
				a_exi = Excidance.GetGeometryRef()
				#area should be calculated from the elements in the excidance polygons
				a_Whon = Whon.GetGeometryRef()
				#Area = a_Whon.GetArea()
				inter = a_exi.Intersection(a_Whon)
				#print inter
				#print inter #print the boundaries of polygon of intersection
				if inter.GetArea() != 0.0:
					#print 'Area of intersection', inter.GetArea()
					Kategorie = Whon.GetField(1)
					#print Kategorie
					LANDUSECODE = Kategorie
					AREA = inter.GetArea()
					polygon = ogr.CreateGeometryFromWkt(str(inter))
					#polygon = ogr.Geometry(ogr.wkbPolygon)
					# Set geometry
					feature.SetGeometryDirectly(polygon) #[14]
					#print ID
					# Set attributes
					feature.SetField('ID', ID) #[15]
					feature.SetField('AREA', AREA)
					feature.SetField('CATEGORIE', LANDUSECODE)
					dst_layer.CreateFeature(feature)
					#print 'The layer',dst_filename,'has', ID, 'conflict regions'
					#print feature
					ID = ID + 1
					#print ID
	#print 'Creating output %s of format %s.' % (dst_filename, format)
	#print 'The layer', dst_filename, 'has', ID, 'conflict regions'
	# Clean up
	feature.Destroy()#[17]
	destinationSource.Destroy() #[18]
	#pdb.set_trace()
	#print os.getcwd()
	os.chdir('../..')
	### clean TEMP dir from the TEMP polygon's
	#a=os.listdir('TEMP/'+nameOfContaminant)
	#for i in range(len(a)):
	#    os.remove('TEMP/'+nameOfContaminant+'/'+a[i])
	#os.rmdir('TEMP/'+nameOfContaminant)
	return dst_filename
#destinationSource = None

###Max, 2009-12-06: insert correct paths to header.map
###Oz, 19.03.2010: fix unix path names
def ReplaceHeadermapPaths(pathtoheadermap):
	'''This function replaces all project specific paths by path according
	to the actual location of the project data.
	Searches for occurrences of a '\".*DATA' string and puts the actual path in there.
	The header.map file is supposed to be in the /DATA/ directory of the 
	actual project.
	Input: 
	pathtoheadermap - string - containing full filepath to *.map FILE 
	
	output:
	map File - the file is the same input file with all paths updated 
	to fit current system configuration. 
	'''
	import fileinput #re, 
	#try:
	for line in fileinput.FileInput(pathtoheadermap, inplace=1, backup='.bak'):
		phrase = re.search('\\".*DATA', line) # ('\"(.*)DATA', line)
		#replaceterm = phrase.group(0)
		if (phrase) != None:
			replaceterm = phrase.group(0)
			#replaceterm = '/'
			if (replaceterm) in line:
				line = line.replace(replaceterm, '"' + os.getcwd() + '/DATA')
				#line = line.replace(replaceterm, '\\')
		elif (phrase)!= None:
			phrase =  re.search('\\\'.*DATA', line) # ('\"(.*)DATA', line)
			if (phrase) != None:
					replaceterm = phrase.group(0)
					if (replaceterm) in line:
						line = line.replace(replaceterm, '\'' + os.getcwd() + '\DATA')
		print line,
	#except:
		#print "Unexpected error in function ReplaceHeadermapPaths:", sys.exc_info()[0]

###code to create mapfile
def CreateMAP(aktszenario, aktlayout, ExceedanceLayers, \
		mapheader='DATA/header.map', maptail='DATA/tail.map', aoi='DATA/area_of_interest.shp'):
	''' This function creates a map file to be used by the system command shp2img, which produces
		the final jpeg image to be viewed by MMS application. 
		Input:
		aktscenario - string - name of the scenario layout directory to process, i.e. 'SzenarioA'.
		aktlayout - string - name of the layout to process, i.e. 'ScALayout1'.
		ExceedanceLayers - list of strings - containing the names of the contaminant 
		  Exceedances shape files, i.e. ['B_PAK', 'W_PCE', 'W_TCE']
		
		default inputs:
		Basic components of the Scenario.map to be created, to edit these files, 
		consult the documentation of MapServer Mapfile format. 
		
		mapheader='DATA/header.map', maptail='DATA/tail.map' 
		
		area of interest shape file:
		aoi='DATA/area_of_interest.shp'
		
		output:
		map File - a map file describing all the layers and info for the jpeg. The file
		  is created in the scenario direction and is named according to the aktlayout 
		  e.g. ScALayout1.map.
		'''
	#read header write it here
	#add extent
	#read the extent from the area_of_interest.shp
	###                "minx           miny            maxx           maxy"
	#EXTENT         3366307.654280 5813354.151378 3367797.654280 5814955.151378
	###  
	aoiSource = aoi #'DATA/area_of_interest.shp' 
	#pdb.set_trace()
	driver = ogr.GetDriverByName('ESRI Shapefile')
	aoiSource = driver.Open(aoiSource, 0)
	layer = aoiSource.GetLayer()
	extent = layer.GetExtent()
	#print 'Extent of aoi:', extent
	#print 'LOWER LEFT CORNER:', extent[0], extent[2]
	#print 'UPPER RIGHT CORNER:', extent[1], extent[3]
	
	# Now we need convert the extent tuple to the above structure, and add some with space so 
	# we increase maxx, maxy and also decrease miny and minx, this is done by adding 110 m 
	# to max, maxy and removing 110 from minx, miny
	modifiedExtent = [extent[0] - 200, extent[2] - 150, extent[1] + 170, extent[3] + 110]
	# Start of LAYER DEFINITIONS ---------------------------------------------
###
#add layer of land uses
#loop over kategories - from the DBF
#should add in the project ini, description of RGB color for each land use
#but for now this should do
	#print 'Test if colors are here....\n\n\n \t : ', LUTcolours
	#print '\n\n Test if lut-indexex are here....\n\n\n \t : ', LUTName
	#print '\n\n Test if n_Landuses are here....\n\n\n \t : ', n_Landuses
   #each element in the list is line in the file
	AktlayoutLayer = ['LAYER \n', '\tNAME \t\t' + '\'' + aktlayout + '\'', \
			'\n\tTYPE POLYGON \n', \
			'\tDATA \t ' + '\'' + os.getcwd() + '/' + aktscenario + '/'  \
			+ aktlayout + '.shp' + ''' ' ''' '\n', '\tSTATUS DEFAULT \n', '\tTRANSPARENCY 100 \n'] 
	#creation of Conflict layer should be in loop, to add conflict layer for each contaminant
	#ConflictLayer=[]
	#read in the szenario's layout shape file:
	szenarioSource = aktscenario + '/' + aktlayout + '.shp'
	driver = ogr.GetDriverByName('ESRI Shapefile')
	szenarioSource = driver.Open(szenarioSource, 0)
	layer = szenarioSource.GetLayer()
	feature = layer.GetNextFeature()
	# find total number of fields in feature
	numOfFields = feature.GetFieldCount()
	#find all Land Uses and Codes, keep house order so append only non existing values !
	LandUses=[]
	Colors=[]
	Codes=[]
	for i in range(0, len(LUTName)): #define LUTnames and colors for jpeg output
		#print 'i, LUTName: ', i, LUTName[i], '\n'
		LandUses.append(LUTName[i])
		Codes.append(str((i+1)*100))
		tempcolor=''
		for j in range( 0, 3 ):
			tempcolor = tempcolor +' '+ (LUTcolours[i][j])
		Colors.append(tempcolor)
	#get field names
	FieldNames = []
	for i in  range(numOfFields):
		fd = feature.GetFieldDefnRef(i)
		FieldNames.append(fd.GetName())
	#print FieldNames
	AktlayoutLayer.append('CLASSITEM\t\'' + FieldNames[1] + '\'\n')
	#print AktlayoutLayer
	classDef = ('CLASS\n', 'Name \t', '\nEXPRESSION \t', '\nSTYLE \n', 'SYMBOL\t0\n', \
			'OUTLINECOLOR    0 0 0\n', 'COLOR \t', '\nEND\n', 'END\n')
	allClasses = []
	b = range(len(classDef))
	for i in range(len(Codes)):
		TempClass = list(classDef)
		#add the appropriate data to the classDef
		#add name
		TempClass[1] = TempClass[1] + '\'' + LandUses[i] + '\''
		TempClass[2] = TempClass[2] + '\'' + str(Codes[i]) + '\''
		TempClass[6] = TempClass[6] + Colors[i]
		#add the Temp class to allClasses
		allClasses = allClasses + TempClass
	#close layout layer
	allClasses.append('\nEND')
	#Done setting up all classes
	#Now create the layer for Exceedance conflicts
	#do this only if exceedance layers were passed to the function (ExcidanceLayers[])
	if ExceedanceLayers:        
		ExceedanceDef = ('\nLAYER\n', '\tNAME\t', ' \n\tTYPE    POLYGON\n', '\tDATA    ', \
				'\n\tSTATUS off\n', '\tTRANSPARENCY 100\n', '\tCLASS\n', \
				'\tNAME\t' + '\t' + '\'Exceedance  ' , '\tSTYLE\n', \
				'\t\tSYMBOL\t'   '\'hatch\'  \n', \
				'\t\tCOLOR 0 0 0 \n \t\tOUTLINECOLOR     0 0 0 \n', \
				'\t\tANGLE 45 \n \t\tSIZE 10 \n \t\tWIDTH 1\n', \
				'\tEND \n END \n END \n')
		ExceedanceHolder = []
		for i in range(len(ExceedanceLayers)):
			ExceedanceLayout = list(ExceedanceDef)
			ExceedanceLayout[1] = ExceedanceLayout[1] + '\'' + ExceedanceLayers[i] + '\''
			#add the data source
			ExceedanceLayout[3] = ExceedanceLayout[3] + '\t' + '\'' + os.getcwd() \
						 + '/' + aktscenario + '/' + aktlayout + '/' + ExceedanceLayers[i] + '/'\
						 + ExceedanceLayers[i] + '.shp' + '\''
			ExceedanceLayout[7] = ExceedanceLayout[7] + ExceedanceLayers[i][2:] + '\'\n'
			ExceedanceHolder = ExceedanceHolder + ExceedanceLayout
	#collect all pieces together
	mapFileDATA = []
	#adapt path of header.map file in case the project has been relocated:
	ReplaceHeadermapPaths(os.getcwd() + '/' + maptail)
	ReplaceHeadermapPaths(os.getcwd() + '/' + mapheader)
	#read in headerm
	headerfile = open(mapheader, 'r')
	tailfile = open(maptail, 'r')
	headerData = headerfile.readlines()
	tailData = tailfile.readlines()
	headerfile.close()
	tailfile.close()
	#do the same for the file 'tail.map':
	#add layers of land uses
	extentstrList = []
	#print modifiedExtent
	for i in range(len(modifiedExtent)):
		extentstrList.append(str(modifiedExtent[i]))
	extentStr = ''
	for i in range(len(extentstrList)):
		extentStr = extentStr + extentstrList[i] + '\t'
		#print extentStr
	headerData.insert(34, str('\tEXTENT\t' + extentStr))
	#print len(headerData), len(AktlayoutLayer)
	mapFileDATA = headerData + tailData + AktlayoutLayer + allClasses
	#add Exceedance
	if ExceedanceLayers:
		mapFileDATA = mapFileDATA + ExceedanceHolder
	##if no tail.map file
	##don't forget to close the map file with END
	#else:
		#mapFileDATA=mapFileDATA+'\nEND'
	mapFileDATA = mapFileDATA + ['\nEND']
	mapFileName = aktscenario + '/' + aktlayout + '/' + aktlayout + '.map'
	mapFile = open(mapFileName, 'w')
	#mapFile.writelines(AktlayoutLayer)
	mapFile.writelines(mapFileDATA)
	mapFile.close()
	return mapFileName

def convertUnix2DosText(filenames):
	''' Converts a list of filenames from Unix-Text Type
	to DOS-Text type by replacing line endings from
	\n to \r\n.
	input:
	filenames - list - contains the list name of all files in directory. 
	
	output:
	All text files in 'filenames' converted to DOS Text Format.
	'''
	for x in filenames:
		if os.path.isdir(x):
			#print x, "Verzeichnis! Skip..."
			continue
		data = open(x, "rb").read()
		if '\0' in data:
			#print x, "Binärdatei! Skip..."
			continue
		newdata = re.sub("\r?\n", "\r\n",data)
		if newdata != data:
			print x, ' ...konvertiert von UNIX nach DOS Text Format!'
			f= open(x, "wb")
			f.write(newdata)
			f.close()
