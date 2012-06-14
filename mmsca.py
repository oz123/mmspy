from osgeo import ogr
from osgeo import osr
from osgeo import gdal
from osgeo import gdal_array
import csv
import os, ConfigParser, sys, re
import shutil
from matplotlib import nxutils
import numpy as np

config = ConfigParser.RawConfigParser()

def readprojectfile(projectinifile): 
    """
    Read the Project.ini file which determines  
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
    LUTColours - 
    LUTName -
    n_Landuses, 
    no_contaminants, 
    selcont
    """
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

class project():
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
        self.n_Landuses = 0
        self.no_contaminants = 0
        self.aoiraster = ""
        self.LUTcolours = ['']
        self.LUTName = ['']
        self.scen_landuseratio = ['']
        self.selcont = ['']
    def getconfig(self, projectini):
        config.readfp(open(projectini))
        self.aktscenario = config.get('DSS_project', 'AktSzenario')
        self.aktlayout = config.get('DSS_project', 'AktLayout')
        self.aoiraster = config.get('DSS_project', 'pathstandort')
        self.n_Landuses= int(config.get('DSS_project', 'n_Landuses'))
        self.no_contaminants = int(config.get( self.aktscenario, 
                                                'anzahlschadstoffe'))
        self.scen_landuseratio *= self.n_Landuses
        self.LUTName *= self.n_Landuses
        self.LUTcolours *= self.n_Landuses
        for i in range(1, self.n_Landuses+1):
            idx=str(i)
            self.scen_landuseratio[i-1] = config.get(self.aktscenario, 
                                    'scen_landuseratio(' + idx + ')') 
            self.LUTName[i-1] = config.get('DSS_project',
                                           'LUTName(' + idx + ')')
            self.LUTcolours[i-1] = [
                config.get('DSS_project','LUTcoloursR(' + idx + ')'), \
                config.get('DSS_project','LUTcoloursG(' + idx + ')'), \
                config.get('DSS_project','LUTcoloursB(' + idx + ')')  ]
        # MAX: note here , len(selcont) = no_contaminants
        #      while in the original readprojectfile 
        #      len(selcont) = no_contaminants - 1
        self.selcont *= self.no_contaminants            
        for i in range(self.no_contaminants):
            idx=str(i+1)
            self.selcont[i]=config.get(self.aktscenario,
                                            'selcont(' + idx + ')') 
    def cleanup(self):
        """
        Clean up files from previous runs ...
        """
        suffixes = ['.cost', '.cosk', '.snh', '.WE', 'opttmep']
        mmsfiles2keep = ['']*len(suffixes)
        for x, suffix in enumerate(suffixes):
            mmsfiles2keep[x]=self.aktscenario+'/'+self.aktlayout \
            +'/'+self.aktlayout+suffix
        # make a copy of needed files in mmsfiles2keep
        for x in mmsfiles2keep:
            if os.path.isfile(x):
                fname=x.replace(self.aktlayout + '/' + self.aktlayout, 
                self.aktlayout)
                shutil.copy2(x, fname) 
            elif os.path.isdir(x):
                dirname = x.replace(self.aktlayout + '/opttemp', '/opttemp')
                shutil.copytree(x, dirname)
        #delete the dir including files and subdirs
        shutil.rmtree(self.aktscenario + '/' + self.aktlayout, ignore_errors=True)
        #create the layout dir, again
        if not os.path.isdir(self.aktscenario + '/' + self.aktlayout):
            os.mkdir(self.aktscenario + '/' + self.aktlayout)
        #move back the *.cost, *cosk , *snh and *WE files to 
        #the aktlayout dirproj=mmsca.project()
        for x in mmsfiles2keep:
            fname = x.replace(self.aktlayout + '/' + self.aktlayout, self.aktlayout)
            dirname=x.replace(self.aktlayout + '/opttemp', '/opttemp')
            if os.path.isfile(fname):
                shutil.copy2(fname, x) 
            elif os.path.isdir(dirname):
                shutil.copytree(dirname, x)
    
    def creatMask(self):
        """
        Create Mask Raster
        """
        print "Not implemented yet"

class ASCIIRaster():
    def __init__(self):
        self.ncols      =  100
        self.nrows      =  100
        self.xllcorner  =   0
        self.yllcorner  =   0
        self.extent = []
        self.cellsize  =    10
        self.NODATA_value = -9999
        self.data = np.zeros((0,0))
    
    def Reader(self, filename):
        """
        read an asci file, store a numpy array containing all data
        """
        dataset = gdal.Open(filename)
        self.extent = dataset.GetGeoTransform()
        self.ncols = dataset.RasterXSize #ncols
        self.nrows = dataset.RasterYSize #nrows
        # TODO: fix me, extent here is wrong!
        # This extents are correct for shape files
        # note for future, when reading shape files with shapelib
        # In [43]: mask.extent
        #Out[43]: (2.555355548858813, 6.40833997726444, 49.49721527099638, 51.50382232666015)
        
        #In [44]: r.bbox
        #Out[44]: [2.555355548858813, 49.49721527099638, 6.40833997726444, 51.50382232666015]
        # the elements 2 and 3 of extent ans bbox are swapped!
        self.xllcorner = self.extent[0]
        self.yllcorner = self.extent[2]
        self.xurcorner = self.extent[1] 
        self.yurcorner = self.extent[3]
        self.data = gdal_array.DatasetReadAsArray(dataset)
        
    def fillRasterPoints(self, Xres, Yres):
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
        self.Xrange = np.arange(self.xllcorner+0.5*Xres, 
                           self.xurcorner+0.5*Xres,Xres) 
        self.Yrange = np.arange(self.yllcorner+0.5*Yres, 
                           self.yurcorner+0.5*Yres,Yres)
        X, Y = np.meshgrid(self.Xrange, self.Yrange)
        self.rasterpoints = np.column_stack((X.flatten(), 
            Y.flatten())) 
    
    def Writer(self, dst_filename, array, topLeftOrigin, ewRes, nsRes, proj=31468): 
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
        format = "MEM"
        driver = gdal.GetDriverByName(format)
        dst_ds = driver.Create(dst_filename, len(array[0]), len(array), \
                1, gdal.GDT_Float32)
        extent = (topLeftOrigin[0], ewRes, 0, topLeftOrigin[1], 0, nsRes )
        dst_ds.SetGeoTransform(extent)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(proj) 
        dst_ds.SetProjection(srs.ExportToWkt())
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.GetRasterBand(1).WriteArray(array)
        format = 'AAIGrid'
        driver = gdal.GetDriverByName(format)
        dst_ds_new = driver.CreateCopy(dst_filename, dst_ds)
        dst_ds = None

class MaskRaster(ASCIIRaster):
    """
    Define a Mask raster containing Zero for points outside the area
    of interest, and One for points inside the area of interest.
    
    test usage:
    import time
    a = time.time()
    test = MaskRaster()
    test.getAreaofInterest("DATA/area_of_interest.shp")
    test.fillRasterPoints(10,10)
    test.getMask()
    test.Writer("test_mask.asc", test.mask, (test.extent[0], test.extent[3]),10,10)
    print "finished in:", time.time() - a , "sec"
    """
    def getAreaofInterest(self,aoi_shp_file):
        dataSource = aoi_shp_file
        driver = ogr.GetDriverByName('ESRI Shapefile')
        dataSource = driver.Open(dataSource, 0)
        layer = dataSource.GetLayer()
        self.extent = layer.GetExtent()
        print 'Extent des Bewertungsgebiets (area_of_interest):', self.extent
        print 'UL:', self.extent[0], self.extent[3]
        print 'LR:', self.extent[1], self.extent[2]
        self.xllcorner = self.extent[0]
        self.yllcorner = self.extent[2]
        self.xurcorner = self.extent[1] 
        self.yurcorner = self.extent[3]
        area = layer.GetFeature(0)
        geometry = area.GetGeometryRef()
        boundary_raw = str(geometry.GetBoundary())
        boundary = boundary_raw[12:-1]
        #some processing to convert boundary raw to a usefull form !
        #remove tail and head
        boundary = boundary.split(',')#convert the string to a list of strings
        #convert the string to a list of lists (couples of x y coordinates)      
        for x, point in enumerate(boundary):
            #pointx, pointy = point.split()
            #boundary[x] = float(pointx), float(pointy)
            boundary[x] = point.split()
        # print boundary
        # TODO: THIS CODE is Correct only if AOI has one polygon
        #       We need a fix in case we have multiple polygons
        np.set_printoptions(precision=18)
        self.boundingvertices=np.asarray(boundary, dtype=np.float64)
        #matplotlib.nxutils.points_inside_poly(mask.gridcenters2, [[0,0],[0,1],[1,0]]
    
    def getMask(self, boundingVertices):
        """
        getMask takes a list of points, and a list of vertices, 
        representing the polygon boundaries, and returns a Boolean 
        array with True for points inside the polygon and False for 
        points outside the polygon.
        TODO: make this function use multithreading so we can go
        a bit faster !
        """
        self.mask=nxutils.points_inside_poly(self.rasterpoints, 
                    boundingVertices)
        #self.mask.resize(self.Yrange.size, self.Xrange.size)
        #self.mask = np.flipud(self.mask)


class LandUseShp():
    """
    Get landuses from shapefile, create landuses raster
    """
    def __init__(self,shapeFilePath):
        driver = ogr.GetDriverByName('ESRI Shapefile')
        self.LandUses = []
        self.Codes = [] #Store Categories
        self.NPolygons = 1 # number of polygons in layer
        self.dataSource = driver.Open(shapeFilePath, 0)
        self.layer = self.dataSource.GetLayer()
        self.NPolygons = self.layer.GetFeatureCount()
        self.Boundaries = [""]*self.NPolygons   # store vertices of each polygon      
        for polygon in range(self.NPolygons):
            area = self.layer.GetFeature(polygon)
            Name = area.GetField(0)
            Kategorie = area.GetField(1)
            self.LandUses.append(Name)
            self.Codes.append(Kategorie)
            geometry = area.GetGeometryRef()
            boundary_raw = str(geometry.GetBoundary())
            #remove tail and head
            boundary = boundary_raw[12:-1]
            boundary = boundary.split(',')
            #print boundary
            #lenboundary = len(boundary)
            #convert each coordinate from string to float
            for x, point in enumerate(boundary):
                boundary[x] = point.split()
            boundingvertices=np.asarray(boundary, dtype=np.float64)
            self.Boundaries[polygon] = boundingvertices

    def Rasterize(self, Xres, Yres, rasterFilePath=None):
        """
        raster object does not needs to be created outside
        """
        raster = MaskRaster()
        raster.extent = self.layer.GetExtent()
        raster.xllcorner = raster.extent[0]
        raster.yllcorner = raster.extent[2]
        raster.xurcorner = raster.extent[1] 
        raster.yurcorner = raster.extent[3]
        raster.fillRasterPoints(Xres, Yres)
        raster.data = np.zeros(raster.rasterpoints.shape[0])
        print raster.data.shape
        for boundary, landuse, code in zip(self.Boundaries, 
                                    self.LandUses, self.Codes):
            vmask=nxutils.points_inside_poly(raster.rasterpoints, boundary)
            raster.mask = np.column_stack([vmask,vmask])
            
            #deleted_indexes=np.where(raster.mask==True)
            #raster.rasterpoints=np.delete(raster.rasterpoints, 
            #                        np.where(raster.mask==True),0)
            
            raster.rasterpoints=np.ma.array(raster.rasterpoints, mask=raster.mask, fill_value=[0,0])
            vmask = vmask*code
            print vmask.shape
            #raster.mask = np.ma.masked_equal(raster.mask, 0)
            #raster.mask.set_fill_value(-9999)
            raster.data = raster.data + vmask
            raster.rasterpoints=raster.rasterpoints.filled()
            # insert to deleted indecies bogus points so next time
            # raster.mask is the same size as the previous step
            
        if rasterFilePath:
            """
            PROBLEM: This still lives 0 intead of no data"
            """
            raster.data.resize(raster.Yrange.size, raster.Xrange.size)
            #raster.data.reshape(X)
            raster.data=np.flipud(raster.data)
            np.putmask(raster.data, raster.data==0, -9999)
            raster.Writer(rasterFilePath, raster.data, (raster.extent[0], raster.extent[3]),Xres,Yres)
        else: return raster



shapeFilePath = "SzenarioA/ScALayout1.shp"
A=LandUseShp(shapeFilePath)
A.Rasterize(10, 10,"scal1.asc")
shapeFilePath = "SzenarioA/ScALayout2.shp"
A=LandUseShp(shapeFilePath)
A.Rasterize(10, 10,"scal2.asc")
print "a"
hapeFilePath = "SzenarioB/ScBLayout1.shp"
A=LandUseShp(shapeFilePath)
A.Rasterize(10, 10,"scbl1.asc")
shapeFilePath = "SzenarioB/ScBLayout2.shp"
A=LandUseShp(shapeFilePath)
A.Rasterize(10, 10,"sccl1.asc")
print "d"
shapeFilePath = "SzenarioC/ScCLayout1.shp"
A=LandUseShp(shapeFilePath)
A.Rasterize(9, 9,"sccl1.asc")
A.Rasterize(5, 5,"sccl15.asc")
shapeFilePath = "SzenarioC/ScCLayout2.shp"
A=LandUseShp(shapeFilePath)
A.Rasterize(10, 10,"sccl2.asc")


import sys
sys.exit()

#MAX TERMIN MONTAG 16:30 25.Juni

test = MaskRaster()
#test.getAreaofInterest("DATA/area_of_interest.shp")
#test.extent=A.layer.GetExtent()

