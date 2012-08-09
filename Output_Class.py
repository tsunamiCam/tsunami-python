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

from Grid_Class import GridPG


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
  
        self.domain_grp = self.grid_dataset.groups['domain'] 
                
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
        
        self.elevation = self.domain_grp.variables['xyz'][:,2]


        #get the side unit vectors
        self.sdx = self.domain_grp.variables['sdx'][:]
        self.sdy = self.domain_grp.variables['sdx'][:]


        
        self.__init_buildings_dictionary__()

        #connect to the POSTGIS database
        self.grid_pg = GridPG(grid_dbname=grid_dbname, epsg = grid_db_epsg) 
        self.slen = self.grid_pg.get_side_lengths()
        #self.element_areas = self.grid_pg.get_element_areas()

        


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

            #get the average elevation of the building
            zAvg = 0
            for n in building_nodes[i]:
                self.elevation[n-1]
                zAvg = zAvg + self.elevation[n-1]    
            zAvg = zAvg/len(building_nodes[i]) 
            
            self.buildings_dict[id] = {'nodes': building_nodes[i].tolist(), 'elements': building_elements[i].tolist(),
                                        'sides': building_sides[i].tolist(), 'perimeter': building_perimeter[i], 'elevation':zAvg}
            

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

    def get_output_at_building(self,id):
        
        '''
        Given the building id, return the uv and eta data associated with building
        '''
        #Get the elements,sides and nodes that make up the building
        nodes = self.buildings_dict[id]['nodes']
        elements = self.buildings_dict[id]['elements']
        sides = self.buildings_dict[id]['sides']
        perimeter = self.buildings_dict[id]['perimeter']

        nsides = len(sides)
        nnodes = len(nodes)
        nelements = len(elements)
        
        speedMax_dict = {}
        etaMax_dict = {}
        speedMax_dict['nodes'] = []
        speedMax_dict['elements'] = []
        speedMax_dict['sides'] = []
        speedMax_dict['sides2'] = []

        etaMax_dict['nodes'] = []
        etaMax_dict['elements'] = []
        etaMax_dict['sides'] = []


        speedAvg_dict = {}
        etaAvg_dict = {}
        speedAvg_dict['nodes'] = []
        speedAvg_dict['elements'] = []
        speedAvg_dict['sides'] = []
        etaAvg_dict['nodes'] = []
        etaAvg_dict['elements'] = []
        etaAvg_dict['sides'] = []       
        
        etaMaxMax = 0
        speedMaxMax = 0


        time = self.building_nodes_grp.variables['time'][:]

        iMIN = 0
        iMAX = len(time)
        i = iMIN
        #for step in time:
        while i < iMAX:
            speedMax = 0
            etaMax = 0
            speedAvg = 0
            etaAvg = 0
            uAll = 0
            vAll = 0
            uvLIST = []
    
            nodesU = self.building_nodes_grp.variables['uv'][:,i,0]
            nodesV = self.building_nodes_grp.variables['uv'][:,i,1]
            nodesEta = self.building_nodes_grp.variables['eta'][:,i]


            sidesU = self.building_sides_grp.variables['uv'][:,i,0]
            sidesV = self.building_sides_grp.variables['uv'][:,i,1]

    
            #get the node output
            for n in nodes:
                if n != 0:
                    
                    u = nodesU[self.node_ids_dict[n]]
                    v = nodesV[self.node_ids_dict[n]]
                    
                    uvLIST.append([u,v])
                    
                    uAll = uAll + u
                    vAll = vAll + v
    
                    speed = math.sqrt(u*u + v*v)
                    speedAvg = speedAvg + speed                      
    
                    if speed > speedMax:
                        speedMax = speed
                    
                    
                    eta = nodesEta[self.node_ids_dict[n]] - self.elevation[n-1]
                    etaAvg = etaAvg + eta
    
                    if eta > etaMax:
                        etaMax = eta
             
            '''           
            #compute the unit vector of [uAll,vAll]
            ui = uAll/math.sqrt(uAll*uAll + vAll*vAll)
            vi = vAll/math.sqrt(uAll*uAll + vAll*vAll)
            #    print "%s     %s" % (ui,vi)
            
            #print math.acos(ui)*180/math.pi
            
            k = 0
            mUi = 0
            mVi = 0
            kMax1 = 0
            kMax2 = 0

            while k < len(uvLIST):
                u = uvLIST[k][0]
                v = uvLIST[k][1]
                uvLIST[k][0] = u*ui + v*vi
                uvLIST[k][1] = -u*vi + v*ui
                
                if abs(uvLIST[k][0]) > abs(mUi):
                    mUi = uvLIST[k][0]
                    kMax1 = k
                if abs(uvLIST[k][1]) > abs(mVi):
                    mVi = uvLIST[k][1]
                    kMax2 = k

                
                k+=1    
            
            print "%s     %s   %s   %s     %s" % (math.acos(ui)*180/math.pi, mUi,mVi, nodes[kMax1],nodes[kMax2])
            '''   
                
                
                    
                
                

            speedMax_dict['nodes'].append(speedMax)
            if speedMax > speedMaxMax:
                speedMaxMax = speedMax
                
            etaMax_dict['nodes'].append(etaMax)
            if etaMax > etaMaxMax:
                etaMaxMax = etaMax
           
            speedAvg = speedAvg/nnodes
            speedAvg_dict['nodes'].append(speedAvg)
            
            etaAvg = etaAvg/nnodes
            etaAvg_dict['nodes'].append(etaAvg)            
            
            

            speedMax = 0
            speedAvg = 0
            speedAvgWeighted = 0
            uWeighted = 0
            vWeighted = 0
            
            etaAvg = 0    
            
            for s in sides:
                if s != 0:
                    u = sidesU[self.side_ids_dict[s]]
                    v = sidesV[self.side_ids_dict[s]]

                    speed = math.sqrt(u*u + v*v)
                    speedAvg = speedAvg + speed
                    #slen is degrees!!! must fix
                    if speed > speedMax:    
                        speedMax = speed


                    speedAvgWeighted = speedAvgWeighted + speed*(self.slen[s-1])/perimeter
                    #uWeighted = uWeighted + u*(self.slen[s-1])/perimeter
                    #vWeighted = vWeighted + v*(self.slen[s-1])/perimeter
                    
                    #Calculate the velocity components tangential and normal to the side
                    ut = self.sdx[s-1]*u - self.sdy[s-1]*v
                    un = self.sdy[s-1]*u + self.sdx[s-1]*v
                    uWeighted = uWeighted + abs(ut*(self.slen[s-1])/perimeter)
                    vWeighted = vWeighted + abs(un*(self.slen[s-1])/perimeter)  
        

            speedAvg = speedAvg/nsides
            speedAvgWeighted = math.sqrt(uWeighted*uWeighted + vWeighted*vWeighted)


            speedWeighted2 = uWeighted
            speedMax_dict['sides'].append(uWeighted)

            speedAvg = speedAvg/nsides
            speedMax_dict['sides2'].append(speedAvgWeighted)


                    


            i+=1
        
        return speedMax_dict, etaMax_dict,time, speedAvg_dict, etaAvg_dict            


    def get_building_output(self,id):
        '''
        Given the building id, return the uv and eta data associated with building
        '''
        #get info from the buildings dictionary
        nodes = self.buildings_dict[id]['nodes']
        elements = self.buildings_dict[id]['elements']
        sides = self.buildings_dict[id]['sides']
                
        
        speedMax_dict = {}
        speedAvg_dict = {}
        speedAvgWeighted_dict = {}
        etaAvg_dict = {}
        etaMax_dict = {}

        speedAvg_nodes_dict = {}


        speedAvgMax = {}
        speedAvgMax['nodes'] = zeros(len(self.buildings_dict))
        speedAvgMax['sides'] = zeros(len(self.buildings_dict))
        #speedAvgMax['elements'] = zeros(len(self.buildings_dict))
        
        speedMaxMax = {}
        speedMaxMax['nodes'] = zeros(len(self.buildings_dict))
        speedMaxMax['sides'] = zeros(len(self.buildings_dict))
        
        print len(speedMaxMax['sides'])

        speedAvgMaxWeighted = {}
        speedAvgMaxWeighted['nodes'] = zeros(len(self.buildings_dict))
        speedAvgMaxWeighted['sides'] = zeros(len(self.buildings_dict))

        etaAvgMax = {}
        etaAvgMax['elements'] = zeros(len(self.buildings_dict))
        etaAvgMax['nodes'] = zeros(len(self.buildings_dict))        

        etaMaxMax = {}
        etaMaxMax['elements'] = zeros(len(self.buildings_dict))
        etaMaxMax['nodes'] = zeros(len(self.buildings_dict)) 
        
        time = self.building_sides_grp.variables['time'][:]

        iMIN = 0
        iMAX = len(time)
        i = iMIN
        #for step in time:
        while i < iMAX:
            #get all the UV results for the current time step            
            sidesU = self.building_sides_grp.variables['uv'][:,i,0]
            sidesV = self.building_sides_grp.variables['uv'][:,i,1]
            
            nodesEta = self.building_nodes_grp.variables['eta'][:,i]
            nodesU = self.building_nodes_grp.variables['uv'][:,i,0]
            nodesV = self.building_nodes_grp.variables['uv'][:,i,1]

            elementsETA = self.building_elements_grp.variables['eta'][:,i]

            j = 0   #the building index
            
            print i
            #get the MAX flow speed around each build for the current time step
            for building in self.buildings_dict.iteritems():      
                id = building[0]
                if i == iMIN:
                    speedMax_dict[id] = {'nodes': [],'elements': [],'sides': []}
                    speedAvg_dict[id] = {'nodes': [],'elements': [],'sides': []}
                    speedAvgWeighted_dict[id] = {'nodes': [],'elements': [],'sides': []}    
                    etaAvg_dict[id] = {'nodes': [],'elements': [],'sides': []}
                    etaMax_dict[id] = {'nodes': [],'elements': [],'sides': []}
                    
                nodes = self.buildings_dict[id]['nodes']
                elements = self.buildings_dict[id]['elements']
                sides = self.buildings_dict[id]['sides']
                perimeter = self.buildings_dict[id]['perimeter']
                
                nsides = len(sides)
                nnodes = len(nodes)
                nelements = len(elements)
                
                
                #-----------------------------------------------------------------------------------
                #NODES OUTPUT
                #-----------------------------------------------------------------------------------
                #
                #
                #
                #-----------------------------------------------------------------------------------
                speedMax = 0
                speedAvg = 0
                speedAvgWeighted = 0
                etaAvg = 0
                etaMax = 0
                
                for n in nodes:
                    if n != 0:
                        
                        u = nodesU[self.node_ids_dict[n]]
                        v = nodesV[self.node_ids_dict[n]]

                        speed = math.sqrt(u*u + v*v)
                        speedAvg = speedAvg + speed                      

                        if speed > speedMax:
                            speedMax = speed

                        eta = nodesEta[self.node_ids_dict[n]]
                        etaAvg = etaAvg + eta
        
                        if eta > etaMax:
                            etaMax = eta
                        
                speedAvg = speedAvg/nnodes
                speedAvg_dict[id]['nodes'].append(speedAvg)
                if (speedAvg > speedAvgMax['nodes'][j]):
                    speedAvgMax['nodes'][j] = speedAvg

                speedMax_dict[id]['nodes'].append(speedMax)
                if (speedMax > speedMaxMax['nodes'][j]):
                    speedMaxMax['nodes'][j] = speedMax   
                
                etaAvg = etaAvg/nnodes
                etaAvg_dict[id]['nodes'].append(etaAvg)
                if (etaAvg > etaAvgMax['nodes'][j]):
                    etaAvgMax['nodes'][j] = etaAvg

                etaMax_dict[id]['nodes'].append(etaMax)                        
                if (etaMax > etaMaxMax['nodes'][j]):
                    etaMaxMax['nodes'][j] = etaMax


             
                #-----------------------------------------------------------------------------------
                #SIDES OUTPUT
                #-----------------------------------------------------------------------------------
                #
                #
                #
                #-----------------------------------------------------------------------------------
                
                
                speedMax = 0
                speedAvg = 0
                speedAvgWeighted = 0
                etaAvg = 0    
                
                for s in sides:
                    if s != 0:
                        u = sidesU[self.side_ids_dict[s]]
                        v = sidesV[self.side_ids_dict[s]]
    
                        speed = math.sqrt(u*u + v*v)
                        speedAvg = speedAvg + speed
                        #slen is degrees!!! must fix
                        speedAvgWeighted = speedAvgWeighted + speed*(self.slen[s-1]/perimeter)
                        if speed > speedMax:
                            speedMax = speed
                            
                speedAvg = speedAvg/nsides
                speedAvgWeighted = (speedAvgWeighted/nsides)*perimeter

                speedMax_dict[id]['sides'].append(speedMax)
                speedAvg_dict[id]['sides'].append(speedAvg)                
                speedAvgWeighted_dict[id]['sides'].append(speedAvgWeighted)
                
                if (speedAvg > speedAvgMax['sides'][j]): speedAvgMax['sides'][j] = speedAvg
                if (speedAvgWeighted > speedAvgMaxWeighted['sides'][j]): speedAvgMaxWeighted['sides'][j] = speedAvgWeighted
                if (speedMax > speedMaxMax['sides'][j]): speedMaxMax['sides'][j] = speedMax
                    
                    
                #-----------------------------------------------------------------------------------
                #ELEMENTS OUTPUT
                #-----------------------------------------------------------------------------------
                #
                #
                #
                #-----------------------------------------------------------------------------------
                
                if nelements > 0:
                    etaAvg = 0
                    etaMax = 0  
                    for e in elements:
                        if e != 0:
                            
                            eta = elementsETA[self.element_ids_dict[e]]
                            etaAvg = etaAvg + eta

                            if eta > etaMax:
                                etaMax = eta
    
                    etaAvg = etaAvg/nelements
                    etaAvg_dict[id]['elements'].append(etaAvg)
                    if (etaAvg > etaAvgMax['elements'][j]):
                        etaAvgMax['elements'][j] = etaAvg

                    etaMax_dict[id]['elements'].append(etaMax)                        
                    if (etaMax > etaMaxMax['elements'][j]):
                        etaMaxMax['elements'][j] = etaMax
    

                j += 1
                
                        
            i+=1
        return speedAvgMax,speedMaxMax,etaAvgMax,etaMaxMax,etaMax_dict[500]['elements'],etaMax_dict[500]['nodes'], time[iMIN:iMAX]
 
 
 
 
    def output_row_histogram(self,row_number):
        '''
        Shielding TEST output
        
        Given a row number, get the max speed and eta values for each building in the row
        
        AREA ID: (in order of building row)
        24
        27
        28
        29
        30
        31
        32
        33
        34
        37
        36
        '''
        
        id1 = 0
        id2 = 0
        etaMaxMaxLIST = []
        speedMaxMaxLIST = []
        if row_number == 1:
            id1 = 24
            id2 = 27   
        elif row_number == 2:
            id1 = 27
            id2 = 28
        elif row_number == 3:
            id1 = 28
            id2 = 29
        elif row_number == 4:
            id1 = 29
            id2 = 30        
        elif row_number == 5:
            id1 = 30
            id2 = 31
        elif row_number == 6:
            id1 = 31
            id2 = 32
        elif row_number == 7:
            id1 = 32
            id2 = 33            
        elif row_number == 8:
            id1 = 33
            id2 = 34
        elif row_number == 9:
            id1 = 34
            id2 = 37        
        elif row_number == 10:
            id1 = 37
            id2 = 36        
        else:
            print "Invalid Row number..."    
            return
        
        self.grid_pg.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (id1))        
        a1 = self.grid_pg.cur.fetchall()[0][0]
        
        
        self.grid_pg.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (id2))        
        a2 = self.grid_pg.cur.fetchall()[0][0]
        
                
        self.grid_pg.cur.execute("SELECT b.id FROM buildings AS b \
                        WHERE \
                        (ST_Contains('%s',b.geom) AND NOT ST_Contains('%s',b.geom)) \
                        ORDER BY b.id;" % (a1, a2))
        buildings = self.grid_pg.cur.fetchall()
        
        print "Number of buildings in row = %s" % len(buildings)
        for b in buildings:
            id = b[0]
            speedMax_dict, etaMax_dict,time, etaMaxMax, speedMaxMax = self.get_output_at_building(id)
            etaMaxMaxLIST.append(etaMaxMax)
            speedMaxMaxLIST.append(speedMaxMax)
        
        spdAvg = 0
        for spd in speedMaxMaxLIST:
            spdAvg = spdAvg + spd

        spdAvg = spdAvg/len(speedMaxMaxLIST)

        
        etaAvg = 0
        for eta in etaMaxMaxLIST:
            etaAvg = etaAvg + eta
            
        etaAvg = etaAvg/len(etaMaxMaxLIST)
    
        
        print "Row ID = %s, Speed Avg = %s, ETA Avg = %s " % (row_number, etaAvg, spdAvg)
        
        return spdAvg, etaAvg

        
        
        

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
    
