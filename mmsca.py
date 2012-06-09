from osgeo import ogr
from osgeo import osr
from osgeo import gdal
from osgeo import gdal_array
from numpy import ma, ones,ones_like,  zeros_like, ndenumerate, empty, dtype
import csv
import os, ConfigParser, sys, re
import shutil


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
        self.cellsize  =    10
        self.NODATA_value = -9999
        self.data = np.zeros((0,0))
    def Reader(self, filename):
        """
        read an asci file, store a numpy array containing all data
        """
        dataset = gdal.Open(filename)
        extent = dataset.GetGeoTransform()
        self.ncols = dataset.RasterXSize #ncols
        self.nrows = dataset.RasterYSize #nrows
        self.xllcorner = extent[0]
        self.yllcorner = extent[3]
        self.cellsize  = extent[1]
        self.data = gdal_array.DatasetReadAsArray(dataset)

class MaskRaster(ASCIIRaster):
    """
    Define a Mask raster containing Zero for points outside the area
    of interest, and One for points inside the area of interest.
    """
    #def __init__(self,aoi_shp_file):
    #    print "take in area_of_interest shape file. not implemented yet."
    getAreaofInterest(self,aoi_shp_file):
dataSource = aoi_shp_file
driver = ogr.GetDriverByName('ESRI Shapefile')
dataSource = driver.Open(dataSource, 0)
layer = dataSource.GetLayer()
extent = layer.GetExtent()
print 'Extent des Bewertungsgebiets (area_of_interest):', extent
print 'UL:', extent[0], extent[3]
print 'LR:', extent[1], extent[2]
area = layer.GetFeature(0)
geometry = area.GetGeometryRef()
boundary_raw = str(geometry.GetBoundary())
#need some processing to convert boundary raw to a usefull form !
#remove tail and head
boundary = boundary_raw[12:-1]
boundary = boundary.split(',')#convert the string to a list of strings
#convert the string to a list of lists (couples of x y coordinates)      
for x, point in enumerate(boundary):
    boundary[x] = point.split()
print boundary
b2 = []
#convert each coordinate from string to float
#for i in range(lenboundary):
#    boundary[i][0] = float(boundary[i][0])
#    boundary[i][1] = float(boundary[i][1])
for i in range(len(boundary)):
    b2.append([float(boundary[i][0]), float(boundary[i][1])])
#nicer to use numpy
boundary2=np.asarray(boundary, dtype=float)
