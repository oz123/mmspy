# -*- coding: utf-8 -*-  
"""
  This file is the main code of the conflictAnalysis tool developed in 
  the DSITE Group, Center for Applied Geosciences, University of Tuebingen, 
  Germany.
  The algorithm was created by Max Morio, and written in python by Oz Nahum
  and Max Morio.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
    
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
     
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
  MA 02110-1301, USA.
"""

from osgeo import ogr
from osgeo import osr
from osgeo import gdal
from osgeo import gdal_array
from osgeo import gdalnumeric

import os, ConfigParser, sys
import shutil
from matplotlib import nxutils
import numpy as np
import csv
from collections import OrderedDict


def world2Pixel(geoMatrix, x, y):
  """
  Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
  the pixel location of a geospatial coordinate 
  """
  ulX = geoMatrix[0]
  ulY = geoMatrix[3]
  xDist = geoMatrix[1]
  yDist = geoMatrix[5]
  rtnX = geoMatrix[2]
  rtnY = geoMatrix[4]
  pixel = int((x - ulX) / xDist)
  line = int((ulY - y) / xDist)
  return (pixel, line) 


class Project():
    """
    defines objects to hold project properties
    TODO:
    Document all project properties.
    
    LUTcolours - RGB colors for each land use.
    LUTName - Names list of all land uses.
    selcont - List of selected contaminants.
    
    #TEST project in the shell
    
    import mmsca
    proj = mmsca.project()
    proj = getconfig("DATA/Projekt.ini")
    
    #TEST Project in the main script confanalysis.py
    def parseArgs():
        ... omitted ...
        try:
            proj =  project()
            proj.getconfig(os.path.abspath(sys.argv[1])+'/Projekt.ini')
        except IOError:
            print "ERROR: Could not find the Projekt.ini under " \
            + os.path.abspath(sys.argv[1])
            sys.exit(1)
        ... omitted ...

    """
    def __init__(self):
        self.aktscenario = ""
        self.aktlayout = ""
        self.layoutdir = ""
        self.n_contaminants = 0
        self.aoiraster = ""
        self.LUTcolours = ['']
        self.LUTName = ['']
        self.scen_landuseratio = ['']
        self.selcont = ['']
        self.n_landuses = 0

    def getconfig(self, projectini):
        """
        read Project.ini and populat all the properties of a Project instance.
        """
        config = ConfigParser.RawConfigParser()
        config.readfp(open(projectini))
        self.aktscenario = config.get('DSS_project', 'AktSzenario')
        self.aktlayout = config.get('DSS_project', 'AktLayout')
        self.aoiraster = config.get('DSS_project', 'pathstandort')
        self.n_landuses = int(config.get('DSS_project', 'n_Landuses'))
        self.n_contaminants = int(config.get( self.aktscenario, 
                                                'anzahlschadstoffe'))
        self.scen_landuseratio *= self.n_landuses
        self.LUTName *= self.n_landuses
        self.LUTcolours *= self.n_landuses
        for i in range(1, self.n_landuses+1):
            idx = str(i)
            self.scen_landuseratio[i-1] = config.get(self.aktscenario, 
                                    'scen_landuseratio(' + idx + ')') 
            self.LUTName[i-1] = config.get('DSS_project',
                                           'LUTName(' + idx + ')')
            self.LUTcolours[i-1] = [
                config.get('DSS_project','LUTcoloursR(' + idx + ')'), \
                config.get('DSS_project','LUTcoloursG(' + idx + ')'), \
                config.get('DSS_project','LUTcoloursB(' + idx + ')')  ]
        # MAX: note here , len(selcont) = n_contaminants
        #      while in the original readprojectfile 
        #      len(selcont) = n_contaminants - 1
        self.selcont *= self.n_contaminants            
        for i in range(self.n_contaminants):
            idx = str(i+1)
            self.selcont[i] = config.get(self.aktscenario,
                                            'selcont(' + idx + ')') 
    def cleanup(self):
        """
        Clean up files from previous runs ...
        """
        suffixes = ['.cost', '.cosk', '.snh', '.WE', 'opttmep']
        mmsfiles2keep = ['']*len(suffixes)
        for idx, suffix in enumerate(suffixes):
            mmsfiles2keep[idx] = self.aktscenario+'/'+self.aktlayout \
            +'/'+self.aktlayout+suffix
        # make a copy of needed files in mmsfiles2keep
        for idx in mmsfiles2keep:
            if os.path.isfile(idx):
                fname = idx.replace(self.aktlayout + '/' + self.aktlayout, 
                self.aktlayout)
                shutil.copy2(idx, fname) 
            elif os.path.isdir(idx):
                dirname = idx.replace(self.aktlayout + '/opttemp', '/opttemp')
                shutil.copytree(idx, dirname)
        #delete the dir including files and subdirs
        shutil.rmtree(self.aktscenario + '/' + self.aktlayout, 
            ignore_errors=True)
        #create the layout dir, again
        if not os.path.isdir(self.aktscenario + '/' + self.aktlayout):
            os.mkdir(self.aktscenario + '/' + self.aktlayout)
        #move back the *.cost, *cosk , *snh and *WE files to 
        #the aktlayout dirproj=mmsca.project()
        for idx in mmsfiles2keep:
            fname = idx.replace(self.aktlayout + '/' + self.aktlayout, 
                self.aktlayout)
            dirname = idx.replace(self.aktlayout + '/opttemp', '/opttemp')
            if os.path.isfile(fname):
                shutil.copy2(fname, idx) 
            elif os.path.isdir(dirname):
                shutil.copytree(dirname, idx)


class ASCIIRaster():
    """
    Wrapper class to handle ASCII Rasters
    """
    def __init__(self):
        self.ncols      =  100
        self.nrows      =  100
        self.xllcorner  =   0
        self.yllcorner  =   0
        self.xurcorner  =   0
        self.yurcorner  =   0
        self.extent = []
        self.cellsize  =    10
        self.NODATA_value = -9999
        self.mask = np.array((0, 0))
        self.data = np.array((0, 0))
        self.rasterpoints = np.array((0, 0))
        self.Xrange = np.array((0, 0))
        self.Yrange = np.array((0, 0))
        self.boundingvertices = np.array((0, 0))

    def reader(self, filename):
        """
        read an asci file, store a numpy array containing all data
        """
        dataset = gdal.Open(filename)
        self.extent = dataset.GetGeoTransform()
        self.ncols = dataset.RasterXSize #ncols
        self.nrows = dataset.RasterYSize #nrows
        # This extents are correct for shape files
        # note for future, when reading shape files with shapelib
        # In [43]: mask.extent
        #Out[43]: (2.555355548813, 6.4083399774, 49.49721527098, 51.503826015)
        
        #In [44]: r.bbox
        #Out[44]: [2.5553558813, 49.49721527038, 6.4083399744, 51.503826015]
        # the elements 2 and 3 of extent ans bbox are swapped!
        self.xllcorner = self.extent[0]
        self.xurcorner = self.xllcorner + self.ncols * self.cellsize 
        self.yurcorner = self.extent[3]#self.yllcorner - self.nrows * self.cellsize
        self.yllcorner = self.yurcorner - self.nrows * self.cellsize
        self.data = gdal_array.DatasetReadAsArray(dataset)
        # TODO: check if OK to faltten 
        # self.data = self.data.flatten()
        
    def fillrasterpoints(self, xres, yres):
        """
        Create data points of each raster pixel, based on extents
        and X,Y resolution 
        
        If we have a grid with X from 0 to 2 and Y from 0 to 3:
        
         2.5| +     +      
          2 |
         1.5| +     +
          1 |
         .5 | +     +
            +----------        
              .5 1 1.5 2
          
        The pixel centers are: 
        p1x 0.5, p1y 0.5
        p2x 1.5, p2y 0.5
        p3x 0.5, p3y 1.51
        p4x 1.5, p4y 1.5
        p5x 0.5, p5y 2.5
        p6x 1.5, p6y 2.5
        ...
        """
        self.Xrange = np.arange(self.xllcorner, 
                           self.xurcorner-0.5*xres,xres) 
        self.Yrange = np.arange(self.yllcorner, 
                           self.yurcorner-0.5*xres,yres)
        #print self.Xrange.size, self.Yrange.size
        #self.Xrange = np.arange(self.xllcorner+0.5*xres, 
        #                   self.xurcorner+0.5*xres,xres) 
        #self.Yrange = np.arange(self.yllcorner+0.5*yres, 
        #                   self.yurcorner+0.5*yres,yres)
        self.xpts, self.ypts = np.meshgrid(self.Xrange, self.Yrange)
        self.rasterpoints = np.column_stack((self.xpts.flatten(), 
            self.ypts.flatten())) 
            
        
    def writer(self, dst_filename, array, bottomleft, ew_res, ns_res, 
        proj=31468, Flip=True): 
        """
        This is a generic GDAL function to write ASCII Rasters.
        Here it is an aid function called by ClipPollution and others.
        inputs:
        dst_filename - string - relative path to write the ASCII raster
        array - ndarray, numpy array - a data array to be written in to the file. 
        ewRes - East West Resolution
        nsRes - North South Resolution
        proj - projection (separate .prj file specifying the projection. 
               see  http://spatialreference.org/ref/epsg/31468/.
               The default is Gauss Kruger DHDN zone 4, relevant for
               North Eastern Germany.
        
        output:
        file - This function writes an ASCII raster file to the same path as 
        dst_filename.
        """
        if Flip == True:
            array = np.flipud(array)
            
        gformat = "MEM"
        driver = gdal.GetDriverByName(gformat)
        dst_ds = driver.Create(dst_filename, len(array[0]), len(array), \
                1, gdal.GDT_Float32)
        extent = (bottomleft[0], ew_res, 0, bottomleft[1], 0, ns_res )
        dst_ds.SetGeoTransform(extent)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(proj) 
        dst_ds.SetProjection(srs.ExportToWkt())
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.GetRasterBand(1).WriteArray(array)
        aformat = 'AAIGrid'
        driver = gdal.GetDriverByName(aformat)
        dst_ds = driver.CreateCopy(dst_filename, dst_ds)
        dst_ds = None
        
    def polygonize(self, src_raster, dst_layer, dst_field_num):
        """
        ShapeFile(".", dst_layer, fields={"ID":ogr.OFTInteger, "Exceedance":ogr.OFTInteger})
        src_raster - path to ascii raster
        dst_layer - shapefile to  create items in
        dst_field_num - number of field in db.
        example call:
        ras.polygonize("SzenarioA/ScALayout2/TCE_in_gw_exceedance_bool.asc", A.dst_layer, 1)
        """
        src_ds = gdal.Open(src_raster)
        srcband = src_ds.GetRasterBand(1) 
        prog_func = None 
        maskband = None 
        options = []
        gdal.Polygonize( srcband, maskband, dst_layer, dst_field_num, options,
        callback = prog_func )
        dst_layer.SyncToDisk()
        
class MaskRaster(ASCIIRaster):
    """
    Define a Mask raster containing Zero for points outside the area
    of interest, and One for points inside the area of interest.
    
    test usage:
    import time
    a = time.time()
    test = MaskRaster()
    test.getareofinterest("DATA/area_of_interest.shp")
    test.fillRasterPoints(10,10)
    test.getmask()
    test.Writer("test_mask.asc", test.mask, (test.extent[0], test.extent[3]),10,10)
    print "finished in:", time.time() - a , "sec"
    """
    def getareaofinterest(self, aoi_shp_file):
        """
        Read a shape file containing a a single polygon bounding an area of interest.
        """
        datasource = aoi_shp_file
        driver = ogr.GetDriverByName('ESRI Shapefile')
        datasource = driver.Open(datasource, 0)
        layer = datasource.GetLayer()
        self.new_extent = layer.GetExtent()
        print 'Extent des Bewertungsgebiets (area_of_interest):\n', self.new_extent
        print 'UR:', self.new_extent[1], self.new_extent[3]
        print 'LL:', self.new_extent[0], self.new_extent[2]
        self.new_xllcorner = self.new_extent[0]
        self.new_yllcorner = self.new_extent[2]
        self.new_xurcorner = self.new_extent[1] 
        self.new_yurcorner = self.new_extent[3]
        self.corners = np.array([
            [self.new_xllcorner, self.new_yllcorner], #lower left
            [self.new_xllcorner, self.new_yurcorner], #upper left
            [self.new_xurcorner, self.new_yurcorner], #upper right
            [self.new_xurcorner, self.new_yllcorner], #lower right
            ])      
        area = layer.GetFeature(0)
        geometry = area.GetGeometryRef()
        boundary_raw = str(geometry.GetBoundary())
        boundary = boundary_raw[12:-1]
        #some processing to convert boundary raw to a usefull form !
        #remove tail and head
        boundary = boundary.split(',')#convert the string to a list of strings
        #convert the string to a list of lists (couples of x y coordinates)   
        for idx, point in enumerate(boundary):
            #pointx, pointy = point.split()
            #boundary[x] = float(pointx), float(pointy)
            boundary[idx] = point.split()
        
        # THIS CODE is Correct only if AOI has one polygon
        #       We need a fix in case we have multiple polygons
        # np.set_printoptions(precision=18)
        self.boundingvertices = np.asarray(boundary, dtype = np.float64)
   
    def clip2(self, new_extent_polygon=None):
        """
        Select only data points inside the polygon area of interest, 
        leaves the contminant raster in its original size.
        """
        if new_extent_polygon is not None:
            self.getmask(new_extent_polygon)
        else:
            # corners of the raster clock wise
            corners = np.array([
            [self.new_xllcorner, self.new_yllcorner], #lower left
            [self.new_xllcorner, self.new_yurcorner], #upper left
            [self.new_xurcorner, self.new_yurcorner], #upper right
            [self.new_xurcorner, self.new_yllcorner], #lower right
            ])    
            self.getmask(corners)
        # invert 0 to 1
        self.mask=np.logical_not(self.mask)
        #print self.mask
        #import pdb
        #pdb.set_trace()
        self.mask.resize(self.Yrange.size, self.Xrange.size)
        self.mask = np.flipud(self.mask)
        self.data = np.ma.MaskedArray(self.data, mask=self.mask)
        self.data = np.ma.filled(self.data, fill_value=-9999)
        ### added
        self.cdata = self.data[~self.mask]

    #def clip(self, new_extent_polygon=None):
        #"""
        #clip existing raster to new extents, not implemented yet
        #"""
        #if new_extent_polygon:
            #self.getmask(new_extent_polygon)
        #else:
            ## corners of the raster clock wise
            #self.corners = np.array([
            #[self.new_xllcorner, self.new_yllcorner], #lower left
            #[self.new_xllcorner, self.new_yurcorner], #upper left
            #[self.new_xurcorner, self.new_yurcorner], #upper right
            #[self.new_xurcorner, self.new_yllcorner], #lower right
            #])
            #self.getmask(self.corners)    
            
    def clip_to_cutline(self, xres, yres, new_extent_polygon=None):
        """
        Clip raster to fit inside polygon
        """
        minX, maxX = self.new_extent [0], self.new_extent[1]
        minY, maxY = self.new_extent [2], self.new_extent[3]
        ulX, ulY=world2Pixel(self.extent, minX, maxY)
        lrX, lrY=world2Pixel(self.extent, maxX, minY)
        
        self.getmask(self.corners)
        self.mask=np.logical_not(self.mask)
        self.mask.resize(self.Yrange.size,self.Xrange.size)
        # choose all data points inside the square boundaries of the AOI,
        # replace all other points with NULL
        self.cdata= np.choose(np.flipud(self.mask), (self.data, -9999))
        # resise the data set to be the size of the squared polygon
        self.ccdata=self.cdata[ulY:lrY, ulX:lrX]
        # in second step we rechoose all the data points which are inside the
        # bounding vertices of AOI
    
        # need to re-define our raster points
        self.xllcorner, self.yllcorner = minX, minY
        self.xurcorner, self.yurcorner = maxX, maxY
        self.fillrasterpoints(xres,yres)
        self.fillrasterpoints(10,10)
        self.getmask(self.boundingvertices)
        self.data=self.ccdata
        self.clip2(new_extent_polygon=self.boundingvertices)
        self.data = np.ma.MaskedArray(self.data, mask=self.mask)
        self.data = np.ma.filled(self.data, fill_value=-9999)
        
        
        
    def getmask(self, vertices):
        """
        getMask takes a list of points, and a list of vertices, 
        representing the polygon boundaries, and returns a Boolean 
        array with True for points inside the polygon and False for 
        points outside the polygon.
        !make this function use multithreading so we can go
        a bit faster !
        """
        self.mask = nxutils.points_inside_poly(self.rasterpoints, 
                    vertices)
        #self.mask.resize(self.Yrange.size, self.Xrange.size)
        #self.mask = np.flipud(self.mask)


class ShapeFile():
    """
    Wrapper class around gdal to handle ESRI Shapefiles.
    A = mmsca.ShapeFile("SzenarioA/ScALayout2/","TCE_in_gw_exceedance")
    """
    def __init__(self, dst_dir, dst_layername, fields={"ID":ogr.OFTInteger} ,srs=None):
        """
        Create an empty shape file without any features ...
        dst_layername is the collection of files making the Shapefile bundle,
        e.g. example.shp, example.dbf, example.shx.
        srs - spatial reference system.
        """
        # TODO: if updating file, add fields if they don't exist!
        if os.path.exists(os.path.abspath(dst_dir+"/"+dst_layername+".shp")):
            filename = os.path.abspath(dst_dir+"/"+dst_layername+".shp")
            self.dst_ds = ogr.Open( filename, update=1 )
            if self.dst_ds is not None:
                print "Updating ", filename
            else: 
                print "Could not open ", filename
            self.dst_layer = self.dst_ds.GetLayerByName(dst_layername)
        else:
            Format =  'ESRI Shapefile'
            drv = ogr.GetDriverByName(Format)
            self.dst_ds = drv.CreateDataSource(dst_dir)
            self.dst_layer = self.dst_ds.CreateLayer(dst_layername, srs=srs)
            print "created "+dst_dir+"/"+dst_layername
            print fields
            for k,v in fields.iteritems():
                fd = ogr.FieldDefn(k, v)
                self.dst_layer.CreateField(fd)
                self.dst_layer.SyncToDisk()
    
    def select_polygons(self, condi):
        """
        Select polygons in layers based on a certion condition
        """
        pass 
        
    def intersect(self, in_dir, in_layer, condition_func, condition_field ,
        fields, feature_callbacks=[], dst_dir='.', dst_layer="out", srs=None):
        """
        Intersect two shape files creating a new shape file.
        'fields' is i dictionary containing field name and ogr type.
        By default field creates 'ID' with integer type for every interstect
        polygon. Others can be created too.
        A = mmsca.ShapeFile("SzenarioA/ScALayout2/sa","TCE_in_gw_exceedance")
        A.intersect("SzenarioA/ScALayout2/", "ScALayout2_tgl",r'x==1' ,1, fields={"AREA":ogr.OFTReal},feature_callbacks=[r"intersection.GetArea()"] )

        condition - a call back function which decided wheter to check intersection
        f_callbacks - a list of functions to be excuted
        USAGE:
        def incr_ID(fieldUID):
            return fieldUID+1
        
        exceedance_shp.intersect(landUsesShp, landUsesShp,  fields={"ID":ogr.OFTInteger},
        f_callbacks=[incr_ID])
        
        for each field created we map a callback function that populates the right
        value ...
        """
        
        if not isinstance(fields, OrderedDict):
            print "WARNING:, fields is not an OrderedDict, Can't guarantee order"\
            +" of fields..."
                
        # fancy recursion: the method calls the init method of the object
        # it's cool that is possible, but is it O.K.?
        if dst_layer is not None:
             outfile = ShapeFile(dst_dir,dst_layer,fields=fields, srs=srs)
             outfile.dst_layer.SyncToDisk()        
        infile = ShapeFile(in_dir, in_layer)
        
        feature = ogr.Feature(outfile.dst_layer.GetLayerDefn())
        # select all polygons with condition existing
        positives = []
        
        for item in range(self.dst_layer.GetFeatureCount()):
            featureA = self.dst_layer.GetNextFeature()
            #x = featureA.GetField(condition_field)
            if eval(condition_func):  
                positives.append(featureA)
        #print "Positive polygons:", len(positives)
        for x,polygon in enumerate(positives):
            #import pdb; pdb.set_trace()
            for item in range(infile.dst_layer.GetFeatureCount()):
                featureB = infile.dst_layer.GetNextFeature()
                featureB_geom = featureB.GetGeometryRef() 
                featureA_geom = positives[x].GetGeometryRef()          
                intersection = featureA_geom.Intersection(featureB_geom)
                if intersection.GetArea() != 0.0:
                    new_polygon = ogr.CreateGeometryFromWkt(str(intersection))
                    feature.SetGeometryDirectly(new_polygon)
                    
                    for k,v in zip(fields.keys(),feature_callbacks):
                        feature.SetField(k, eval(v))
                    
                    outfile.dst_layer.CreateFeature(feature)
            infile.dst_layer.ResetReading()
        self.dst_layer.ResetReading()
        outfile.dst_layer.SyncToDisk()

        
class LandUseShp():
    """
    Get landuses from shapefile, create landuses raster
    """
    def __init__(self, shapefilepath, copyfile=None):
        self.filepath=shapefilepath
        driver = ogr.GetDriverByName('ESRI Shapefile')
        self.LandUses = []
        self.Codes = [] #Store Categories
        self.NPolygons = 1 # number of polygons in layer
        self.fields = []
        # driver.Open(path,0) -> open read only
        # driver.Open(path,1) -> open with update option
        if copyfile:
            self.copyfilepath=copyfile
            self.createcopy()
            self.dataSource = driver.Open(self.copyfilepath, 1)
        else: 
            self.dataSource = driver.Open(shapefilepath, 0)
        self.layer = self.dataSource.GetLayer()
        self.NPolygons = self.layer.GetFeatureCount()
        # store vertices of each polygon
        self.Boundaries = [""]*self.NPolygons
        feature=self.layer.GetFeature(0)
        self.fields = feature.keys()
        self.layer.ResetReading()
        for polygon in range(self.NPolygons):
            area = self.layer.GetFeature(polygon)
            name = area.GetField(0)
            category = area.GetField(1)
            self.LandUses.append(name)
            self.Codes.append(category)
            geometry = area.GetGeometryRef()
            boundary_raw = str(geometry.GetBoundary())
            #remove tail and head
            boundary = boundary_raw[12:-1]
            boundary = boundary.split(',')
            #convert each coordinate from string to float
            for idx, point in enumerate(boundary):
                boundary[idx] = point.split()
            boundingvertices = np.asarray(boundary, dtype = np.float64)
            self.Boundaries[polygon] = boundingvertices

    def createcopy(self):
        """
        copy the shape file to output file path
        """
        # check if dirname exists, if not create it
        print "Shape file copied to :", self.copyfilepath
        dest_dir=os.path.dirname(self.copyfilepath)
        if not os.path.exists(dest_dir):
            print dest_dir
            os.makedirs(dest_dir)
        suffixes = ['shp', 'shx',  'dbf']
        for suffix in suffixes:
            shutil.copy2(self.filepath.replace("shp", suffix),
                 self.copyfilepath.replace("shp", suffix))
    
    def addfield(self,fieldname):
        """
        Add field to attribute table
        """
        field_TEXT  = ogr.FieldDefn(fieldname, ogr.OFTString) 
        field_TEXT.SetWidth(20)
        self.layer.CreateField(field_TEXT)
        rslt=self.layer.SyncToDisk()
        
        if rslt != 0:
            print "Error: Could not write new feature ", fieldname , \
                " to attribute table of ", self.copyfilepath
            sys.exit(1)
        self.fields.append(fieldname)
        
    def get_value(self, feature, fieldname):
        """
        Get the field value of a polygon.
        """
        idx=feature.GetFieldIndex(fieldname)
        value=feature.GetField(idx)
        return value

    def set_value(self, feature, fieldname ,value):
        """
        Set the property of a polygon.
        This function should  be used inside a loop, since 
        we iterate over the polygons, for example:
        In [42]: A.__class__
        Out[42]: mmsca.LandUseShp
        In [43]: for i in range(A.layer.GetFeatureCount()):
            A.setthresholds("Test", 222+i)
        This will set  the polygon field "Test" to 222+i.
        """
        idx=feature.GetFieldIndex(fieldname)
        #import pdb 
        #pdb.set_trace()
        #print help(feature.SetField2)
        #print help(feature.SetField)
        
        feature.SetField(idx,value)
        self.layer.SetFeature(feature)
        self.layer.SyncToDisk()
        #print "set ", fieldname, "to ", value
        
    def rasterize_field(self, xres, yres, fieldname="Kategorie",rasterfilepath=None):
        """
        This function creates a raster image of the polygons, based
        on one of the database field (the field has to be numeric). 
        """
        raster = MaskRaster()
        raster.extent = self.layer.GetExtent()
        
        raster.xllcorner = raster.extent[0]
        raster.yllcorner = raster.extent[2]
        raster.xurcorner = raster.extent[1] 
        raster.yurcorner = raster.extent[3]
        
        raster.fillrasterpoints(xres, yres)
        raster.data = np.zeros(raster.rasterpoints.shape[0])
        values = [0]*len(self.Codes)
        for polygon in range(self.NPolygons):
            feature = self.layer.GetFeature(polygon)
            values[polygon] = int(self.get_value(feature, fieldname))
        self.layer.ResetReading()
        for boundary, code in zip(self.Boundaries, values):
            vmask = nxutils.points_inside_poly(raster.rasterpoints, boundary)
            raster.mask = np.column_stack([vmask, vmask])
            #deleted_indexes=np.where(raster.mask==True)
            #raster.rasterpoints=np.delete(raster.rasterpoints, 
            #                        np.where(raster.mask==True),0)
            raster.rasterpoints = np.ma.array(raster.rasterpoints, 
                    mask = raster.mask, fill_value = [0, 0])
            vmask = vmask*code
            #raster.mask = np.ma.masked_equal(raster.mask, 0)
            #raster.mask.set_fill_value(-9999)
            raster.data = raster.data + vmask
            raster.rasterpoints = raster.rasterpoints.filled()
            # insert to deleted indecies bogus points so next time
            # raster.mask is the same size as the previous step
        if rasterfilepath:
            raster.data.resize(raster.Yrange.size, raster.Xrange.size)
            #raster.data.reshape(X)
            #raster.data = np.flipud(raster.data)
            np.putmask(raster.data, raster.data == 0, -9999)
            xleftc,yleftc = raster.extent[0]-yres*.5, raster.extent[3]+yres*0.5  
            xleftc, yleftc = round(xleftc,-1),round(yleftc,-1) 
            raster.writer(rasterfilepath, raster.data, 
                    (xleftc,yleftc), xres, yres)
        return raster

class ZielWerte():
    """
    Class to hold all the values from Zielwert.rtv.
    self.landuses describes the columns in the file, 
    these match one possible land uses in the shape-file.
    Feel free to change, but don't complain if your shape-files
    contain "Category" and this class contains "Kategorie".
    """
    def __init__(self,proj):
        n_Landuses=proj.n_landuses
        try:
            infile_target =  'Zielwertset.rtv'
            infile_target = os.path.abspath(proj.aktscenario)+'/'+infile_target
            targetreader = csv.reader(open(infile_target, "r"))
        except IOError:
            print "Error: Could not find 'Zielwertset.rtv' under " \
            , os.path.abspath(aktscenario) 
            sys.exit(1)
        self.contnames = []
        self.compartments = []
        self.targets_LUT={}
        self.contrasternames = []
        self.boolConflictType = []
        self.iContRealisations = []
        
        # landuses describes the columns in the file, 
        # This match one possible land uses in the shape-file.
        # Feel free to change, but don't complain if your shape-files
        # contain "Category" and this class contains "Kategorie".
        
        self.landuses = [ "Wohnen", "Gewerbe (nicht sensibel)", 
        "Gewerbe (sensibel)","Renaturierung", "Akerland/Forst",
        "Freizeit und Erholung", "Stadtpark", "Gruenriegel",
        "Gruenpuffer" ]
        self.no_contaminants = int(targetreader.next()[0])
        print 'Anzahl Kontaminanten (B und GW): ', self.no_contaminants
        # do the actual parsing of the file
        for row in targetreader:
            self.contnames.append(row[0])
            self.compartments.append(row[1])
            self.targets_LUT[row[0]]=row[2:n_Landuses+2]
            self.contrasternames.append(row[n_Landuses+2])
            self.boolConflictType.append(row[n_Landuses+3])
            self.iContRealisations.append(row[n_Landuses +4])
        #print self.contnames 
        #print self.compartments
        #print self.contrasternames 
        #print self.no_contaminants
        #print self.targets_LUT 
        #print self.boolConflictType
        #print self.iContRealisations 
        

if __name__ == "main":
    FILEPATH = "SzenarioA/ScALayout1.shp"
    A = LandUseShp(FILEPATH,"SzenarioA/ScALayout1/ScALayout1_bla.shp")
    A.addfield("Test")

# test with gdal 1.7.3


def RasterClipper():
    craster = MaskRaster()
    contraster2 = 'PCE_in_gw.aux'
    #contraster="PAK30_B.aux"
    craster.reader("DATA/"+contraster2.replace('aux','asc'))
    xres, yres = craster.extent[1], craster.extent[1]
    craster.fillrasterpoints(xres, yres)
    craster.getareaofinterest("DATA/area_of_interest.shp")
    minX, maxX=craster.new_extent [0], craster.new_extent[1]
    minY, maxY= craster.new_extent [2], craster.new_extent[3]
    
    ulX, ulY=world2Pixel(craster.extent, minX, maxY)
    lrX, lrY=world2Pixel(craster.extent, maxX, minY)
    
    
    craster.getmask(craster.corners)
    craster.mask=np.logical_not(craster.mask)
    
    
    craster.mask.resize(craster.Yrange.size,craster.Xrange.size)
    # choose all data points inside the square boundaries of the AOI,
    # replace all other points with NULL
    craster.cdata= np.choose(np.flipud(craster.mask), (craster.data, -9999))
    # resise the data set to be the size of the squared polygon
    craster.ccdata=craster.cdata[ulY:lrY, ulX:lrX]
    
    craster.writer("ccdata2m.asc",craster.ccdata, (minX+xres*.5, maxY+yres*.5), 10,10,Flip=False)
    
    # in second step we rechoose all the data points which are inside the
    # bounding vertices of AOI
    
    # need to re-define our raster points
    craster.xllcorner=minX
    craster.yllcorner=minY
    craster.xurcorner=maxX
    craster.yurcorner=maxY
    
    craster.fillrasterpoints(10,10)
    craster.getmask(craster.boundingvertices)
    craster.data=craster.ccdata
    craster.clip2(new_extent_polygon=craster.boundingvertices)
    craster.data = np.ma.MaskedArray(craster.data, mask=craster.mask)
    craster.data = np.ma.filled(craster.data, fill_value=-9999)
    # write the raster to disk
    print minX , maxY
    minX , maxY= round(minX,-1), round(maxY,-1)
    craster.writer("ccdata2m_clipped.asc",craster.data, (minX-xres, maxY+yres), 10,10,Flip=False)
    
