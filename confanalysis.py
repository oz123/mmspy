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


"""
This is the main script that runs the conflict analysis.
most of the functios are in the module mmsca.py
"""
import stat,sys,os,shutil 


def parseArgs():
    """
    Do some basic checks about the command line arguments.
    - Verify the Projekt.ini exists.
    - Verify that conflict type is specified.  
    """
    if not os.path.isdir(sys.argv[1]):
        print "ERROR: The first argument must be a directory which holds the" \
        + " Projekt.ini"
        sys.exit(1)
    try:
        open(os.path.abspath(sys.argv[1])+'/Projekt.ini')
    except IOError:
        print "ERROR: Could not find the Projekt.ini under " \
        + os.path.abspath(sys.argv[1])
    try: 
        conflicttype = sys.argv[2]
        try: 
            conflicttype=int(sys.argv[2])
            print "conflicttype ", conflicttype
        except ValueError:
            print "Error: Conflict Type must be an integer"
    except IndexError:      
        print "ERROR: You didn't specify conflict type."
        sys.exit(1)

def main():
    print "main."

if __name__ == '__main__':
    parseArgs()
    main()
