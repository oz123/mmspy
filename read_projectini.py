# module to read MMS project ini 
# 
# Max Morio / July 2009
#
# parameters:
# <- aktscenario = 'scenarioname'
# <- aktlayout = 'layoutname'
# <- layoutdir = 'pathtolayout'
# -> pathtoprojectini 

import ConfigParser, os
config = ConfigParser.RawConfigParser()

def readprojectfile(projectinifile): #(aktscenario, aktlayout, layoutdir, projectinifile):
	#try:
	#open the ini file
	config.readfp(open(projectinifile))
	#print config.items('DSS_project')
	scen_landuseratio = ''	
	#parse for the names
	aktscenario = config.get('DSS_project', 'AktSzenario')
	aktlayout = config.get('DSS_project', 'AktLayout')
	aoiraster = config.get('DSS_project', 'pathstandort')
	layoutdir = os.getcwd() + '/' + aktscenario + '/' + aktlayout
	# read the szenario's land use ratios
	for i in range(1, 10):
		lookfor = 'scen_landuseratio(' + str(i) + ')'
		print lookfor
		scen_landuseratio = str(scen_landuseratio) + ' ' + str(config.get(str(aktscenario), str(lookfor)))
	return aktscenario, aktlayout, layoutdir , scen_landuseratio, aoiraster
	print 'works'
	
	#except:
		
