#Import the required modules


import psycopg2 #@UnresolvedImport
#import errorcodes 

from netCDF4 import Dataset #@UnresolvedImport
from numpy import *
import time
import os
import sys
from osgeo import ogr #@UnresolvedImport
from osgeo import osr #@UnresolvedImport

#module to allow commands to be sent to the command line (i.e. Terminal)
from subprocess import call
from ctypes import byref, cdll

from datetime import datetime, date, time, timedelta


class ReadOutput:
    
    def __init__(self,run_filename, grid_filename ):
        
        
        
        self.run_filename = run_filename
        self.grid_filename = grid_filename
    
        try: self.run_dataset = Dataset(run_filename,'r',format='NETCDF4')
        except Exception, e:
            print "ERROR: %s" % e
            sys.exit()
    
        try: self.grid_dataset = Dataset(grid_filename,'r',format='NETCDF4')
        except Exception, e:
            print "ERROR: %s" % e
            sys.exit()
        
        run_id = getattr(self.run_dataset,"run_id")          
        grid_dbname = getattr(self.run_dataset,"grid_dbname")          
        grid_db_epsg = getattr(self.run_dataset,"grid_db_epsg")          
        buildings_dbname = getattr(self.run_dataset,"buildings_dbname")          
        name = getattr(self.run_dataset,"name")          
        description = getattr(self.run_dataset,"description")          
        user = getattr(self.run_dataset,"user")    
        
        #Create all the BUILDINGS groups

        self.buildings_grp = self.run_dataset.groups['buildings'] 
        self.building_nodes_grp = self.buildings_grp.groups['nodes']          
        self.building_elements_grp = self.buildings_grp.groups['elements']   
        self.building_sides_grp = self.buildings_grp.groups['sides']
      
        #Create all the Keypoints groups
        self.key_points = self.run_dataset.groups['key_points'] 
        self.key_points_nodes = self.key_points.groups['nodes'] 
        self.key_points_elements = self.key_points.groups['elements']
        
        
        #self.number_of_buildings = getattr(self.buildings_grp,"number_of_buildings")
        self.buildings_dict = {}
        
        self.side_ids_dict = {}
        self.node_ids_dict = {}
        self.element_ids_dict = {}
        
        self.__init_buildings_dictionary__()


        


    def __del__ (self):
        """
        Class deconstructor
    
        """
        self.run_dataset.close()
        self.grid_dataset.close()

    def __init_buildings_dictionary__(self):
        '''
        Initialise the buildings dictionary
        
        '''
        building_id = self.buildings_grp.variables['building_id'][:]
        building_perimeter = self.buildings_grp.variables['building_perimeter'][:]
        building_nodes = self.buildings_grp.variables['building_nodes'][:]
        building_elements = self.buildings_grp.variables['building_elements'][:]
        building_sides = self.buildings_grp.variables['building_sides'][:]
        
        i = 0
        for id in building_id:    
            self.buildings_dict[id] = {'nodes': building_nodes[i].tolist(), 'elements': building_elements[i].tolist(),
                                        'sides': building_sides[i].tolist(), 'perimeter': building_perimeter[i]}
            i+=1
        
        side_ids = self.building_sides_grp.variables['id'][:]
        node_ids = self.building_nodes_grp.variables['id'][:]
        element_ids = self.building_elements_grp.variables['id'][:]
                
        #Create dictionary MAP between id and array location for nodes, sides and elements
        i = 0
        for id in side_ids:
            self.side_ids_dict[id] = i
            i+=1
        i = 0
        for id in node_ids:
            self.node_ids_dict[id] = i
            i+=1
        i = 0
        for id in element_ids:
            self.element_ids_dict[id] = i
            i+=1 


        
    def get_building_output(self,id):
        '''
        Given the building id, return the uv and eta data associated with building
        '''
        #get info from the buildings dictionary
        nodes = self.buildings_dict[id]['nodes']
        elements = self.buildings_dict[id]['elements']
        sides = self.buildings_dict[id]['sides']
        
        slen = self.grid_dataset.variables['slen'][:]
        
        
        speedMax_dict = {}
        speedAvg_dict = {}
        speedAvgMax = zeros(len(self.buildings_dict))
        speedMaxMax = zeros(len(self.buildings_dict)) 
        sides_time = self.building_sides_grp.variables['time'][:]
        i = 3500
        #for step in sides_time:
        while i < 3510:
            #get all the UV results for the current time step            
            sidesU = self.building_sides_grp.variables['uv'][:,i,0]
            sidesV = self.building_sides_grp.variables['uv'][:,i,1]
            j = 0
            print i
            #get the MAX flow speed around each build for the current time step
            for building in self.buildings_dict.iteritems():      
                id = building[0]
                if i == 3500:
                    speedMax_dict[id] = {'sides': []}
                    speedAvg_dict[id] = {'sides': []}
    
                    
                nodes = self.buildings_dict[id]['nodes']
                elements = self.buildings_dict[id]['elements']
                sides = self.buildings_dict[id]['sides']
                
                nsides = len(sides)
                nnodes = len(nodes)
                nelements = len(elements)
                
                speedMax = 0
                speedAvg = 0
                for s in sides:
                    if s != 0:
                        u = sidesU[self.side_ids_dict[s]]
                        v = sidesV[self.side_ids_dict[s]]
    
                        speed = math.sqrt(u*u + v*v)
                        speedAvg = speedAvg + speed
                        if speed > speedMax:
                            speedMax = speed
                        
                speedMax_dict[id]['sides'].append(speedMax)
                speedAvg_dict[id]['sides'].append(speedAvg/nsides)
                
                if (speedAvg/nsides > speedAvgMax[j]):
                    speedAvgMax[j] = speedAvg/nsides


                if (speedMax > speedMaxMax[j]):
                    speedMaxMax[j] = speedMax
                
                j += 1
                
                        
            i+=1
            
        return speedAvgMax, speedMaxMax
        #return array(speedMax_dict[id]['sides']), array(speedAvg_dict[id]['sides']), sides_time
                 
                    
                

                    
                
                
            

            
                       
#            n = row[1]['nodes'] 
#            e = row[1]['elements']
#            s = row[1]['sides']
#            perimeter.append(row[1]['perimeter'])
#            nodesAll.extend(n)
#            elementsAll.extend(e)
#            sidesAll.extend(s)
#            if maxNodes < len(n): maxNodes = len(n)
#            if maxElements < len(e): maxElements = len(e)
#            if maxSides < len(s): maxSides = len(s)



        #output_dict = {'nodes': {}, 'elements': {}, 'sides':{}}

                
        #get the time steps array from the NETCDF file
        #sides_time = self.building_sides_grp.variables['time'][:]

            


#
#        i = 0
#        for s in sides:
#            if s == 0:
#                break
#            else:
#                uv = self.building_sides_grp.variables['uv'][self.side_ids_dict[s]]
#                output_dict['sides'][s] = uv
                

                
                
        #return output_dict
                

            
        
        
        

    def key_point_output_to_GMT(self):
        """
        
        """
        startTime = datetime(2011,03,11,14,46,00)
       
        number_of_key_points = len(self.key_points.dimensions['number_of_key_points'])
        
        if number_of_key_points > 0:
            description = self.key_points.variables['description'][:]
            xy = self.key_points.variables['xy'][:]
            name = self.key_points.variables['name'][:]
            time = self.key_points.variables['time'][:] 
            eta = self.key_points_nodes.variables['eta'][:] 

            i = 0
            outfileNames = []
            outfiles = []

            for id in name:
                #create all the output files
                filename = id + "_RicomTS.d"
                outfileNames.append(filename)
                outfiles.append(open(filename, "w"))
            
            tStep = 0
            for t in time:
                currentDateTime = startTime + timedelta(seconds=t)
                datestring = str(currentDateTime)
                datestring = datestring.replace(' ','T')
                i = 0
                while i < number_of_key_points:
                    outfiles[i].write("%s %s\n" % (datestring, str(eta[i][tStep]*100)))
                    i+=1
                tStep+=1

            #close all the output files
            for file in outfiles:
                file.close()    
    
