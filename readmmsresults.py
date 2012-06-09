

# module to read MMS results and write them to csv-files
# 
# Max Morio / August 2009
#
# parameters:
# <- resultsfile = 'path-and-filename of MMS-results'
 

import ConfigParser, os, csv
config = ConfigParser.RawConfigParser()
def readMMSResultsWriteXLS(resultsfile):
	''' reads the SAFIRA II Cost Estimation results file '''
	# define the csv writer object:
	# open the csv file to be written, same name as input file
	# ... with suffix *.csv, delimiter is space in order to not
	# confuse delimiters for columns and numbers! 
	fid=open(resultsfile + '.csv','wb') 
	writer=csv.writer(fid, delimiter=' ')
	# open the *.snh file of the active scenario and layout
	config.readfp(open(resultsfile))
	# print the first section's name in the file to the screen 
	# print config.sections()
	# convert to string
	section2read = config.sections()
	test = section2read[0]
	section2read = test.replace('[','')
	section2read = test.replace(']','')
	print "\n" + section2read + ':'
	# get a list of available options (entries) in the section
	#print config.options(section2read)
	# read the values of the options and write them to an CSV file
	# which can be read by xls
	optionslist = config.options(section2read)
	n = 0 # index for number of option in the according section
	while n < len(optionslist) :
		#print optionslist[n] + ' ' + config.get(section2read, optionslist[n])
		# write to csv:
		writer.writerow([optionslist[n], config.get(section2read, optionslist[n])])
		n= n + 1
	fid.close()
	
	print '--------------------------------------------------'
	print 'Results for MMS-Assessment written to *.csv file :'
	print resultsfile + '.csv'
	print '--------------------------------------------------\n'
	section2read = None
	optionslist = None
	test = None
	
	
	
	
