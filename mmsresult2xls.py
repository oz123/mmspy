#  wrapper for conflict analysis modules
# -*- coding: iso-8859-1 -*-
import sys, os, string, glob, shutil

# import the module for reading the project's ini file 'Projekt.ini'
import read_projectini

# import the module to read the MMS-results file from the active
# layout and scenario
import csv # to read and write to csv-files
import readmmsresults

# get the path we are currently in (main.py)
print sys.path[0]



# read the command line argument (path to project ini file)
#os.chdir(sys.path[0])
try:
	if (len(sys.argv) >1):
		projectinifile = sys.argv[1]
		projectinifile.replace(':', '//')
        projectinifile.replace('\\', '/')
        os.chdir(projectinifile)
except:
    print 'FEHLER: Datei zum Projektpfad nicht übergeben oder nicht vorhanden! '
    print 'Es wird angenommen, dass dieses Skript vom Projektordner aus aufgerufen wurde.\n ####################################'
    projectinifile = os.getcwd() + '/DATA/Projekt.ini' 

# # Define the path to the ini-file 'Projekt.ini
#projectinifile = os.getcwd() + '/DATA/Projekt.ini' 

# This function parses the Projekt.ini and returns... an inistring
inistring = read_projectini.readprojectfile(projectinifile) 

# split the inistring with the returned values
aktscenario = inistring[0]
aktlayout = inistring[1]
#path to the current layout of the active scenario
layoutdir = inistring[2] 
print 'Directory of the layout to be used for conflict analisis:',layoutdir
print 'Name of the layout shapefile:', aktlayout + '.shp'
layoutshape = layoutdir.replace(aktlayout,'') + '/' + aktlayout + '.shp'

#land use ratio (quota in percent) for the predefined land use types
scen_landuseratio=[]
scen_landuseratio[0:10] = string.split(inistring[3])
print 'Land use quota:',scen_landuseratio[0:10]

#check for area of interest:
aoiraster = inistring[4]
print '\n area of interest raster: ',aoiraster


try:
	resultsfile = layoutdir.replace('\DATA','') + '/' + aktlayout + '.WE'
	print 'test: ' + resultsfile
	if os.path.isfile(resultsfile):
		#\
		#	'D:\Andere\Max\Vorträge_Poster\optimization\code\conflictanalysis\projdataexample\SzenarioA\ScALayout2\ScALayout2.WE'
		readmmsresults.readMMSResultsWriteXLS(resultsfile)
	else:
		print "Datei " + resultsfile + " nicht gefunden! Bitte MMS-Bewertung speichern und Pfade prüfen."
except:
	print 'Fehler bei Export Wertermittlung!'
try:
	resultsfile =  layoutdir.replace('\DATA','')  + '/' + aktlayout + '.snh'
	if os.path.isfile(resultsfile):
		#	'D:\Andere\Max\Vorträge_Poster\optimization\code\conflictanalysis\projdataexample\SzenarioA\ScALayout2\ScALayout2.snd'
		readmmsresults.readMMSResultsWriteXLS(resultsfile)
	else:
		print "Datei " + resultsfile + " nicht gefunden! Bitte MMS-Bewertung speichern und Pfade prüfen."
except:
	print 'Fehler bei Export Nachhaltigkeitsbewertung!'
try:
	resultsfile =  layoutdir.replace('\DATA','')  + '/' + aktlayout + '.cost'
	if os.path.isfile(resultsfile):
		#	'D:\Andere\Max\Vorträge_Poster\optimization\code\conflictanalysis\projdataexample\SzenarioA\ScALayout2\ScALayout2.cost'
		readmmsresults.readMMSResultsWriteXLS(resultsfile)
	else:
		print "Datei " + resultsfile + " nicht gefunden! Bitte MMS-Bewertung speichern und Pfade prüfen."
except:
	print 'Fehler bei Export Sanierungs- u. Aufbereitungskosten!'

print '\n ERGEBNISEXPORT fertig! \n SHELL-FENSTER schliessen oder mit "exit" um zurück zum MMS zu gelangen!'
