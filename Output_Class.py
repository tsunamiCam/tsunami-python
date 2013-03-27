#Import the required modules


import psycopg2 #@UnresolvedImport
#import errorcodes 

from netCDF4 import Dataset #@UnresolvedImport
import time
import os
import sys
from osgeo import ogr #@UnresolvedImport
from osgeo import osr #@UnresolvedImport



from Grid_Class import GridPG


from numpy import *
import numpy as np
from scipy.special import erf,erfinv
from numpy import linspace, exp, sqrt,log
import matplotlib.pyplot as plt
from scipy import polyval,polyfit



from numpy import histogram

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
        
        if grid_dbname == "Yuriage1a":
            grid_dbname = "Yuriage1"
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


        
        #self.__init_buildings_dictionary__()

        #connect to the POSTGIS database
        self.grid_pg = GridPG(grid_dbname=grid_dbname, epsg = grid_db_epsg) 
        self.slen = self.grid_pg.get_side_lengths()
        #self.element_areas = self.grid_pg.get_element_areas()
        
        
        self.output_tablename = ""

        


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
        building_type = self.buildings_grp.variables['building_type'][:]
        
        i = 0
        for id in building_id:    

            #get the average elevation of the building
            zAvg = 0
            for n in building_nodes[i]:
                self.elevation[n-1]
                zAvg = zAvg + self.elevation[n-1]    
            zAvg = zAvg/len(building_nodes[i])
            
            nodes = building_nodes[i].tolist()
            j = len(nodes) - 1
            while j >= 0:
                if nodes[j] == 0:
                    nodes.pop(j)
                j =j-1


            elements = building_elements[i].tolist()
            j = len(elements) - 1
            while j >= 0:
                if elements[j] == 0:
                    elements.pop(j)
                j =j-1
                
            sides = building_sides[i].tolist()
            j = len(sides) - 1
            while j >= 0:
                if sides[j] == 0:
                    sides.pop(j)
                j =j-1                
            
            self.buildings_dict[id] = {'nodes': nodes, 'elements': elements,
                                        'sides': sides, 'perimeter': building_perimeter[i], 'elevation':zAvg, 'type': building_type[i]}
            
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






    def create_fragility_probit(self,damage_class, bv_low = -1, bv_high = 1,damage_low = 1, type=1,bin_size = 8,weight=False):

        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename
    
            damage_level = 7
            bin_number = 20
            if type == 1:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and bv >= %s and bv < %s;" % (t, bv_low,bv_high))        
            
            if type == 2:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and s >= %s and s <= %s;" % (t, bv_low,bv_high))        
    
            if type == 3:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and m >= %s and m <= %s;" % (t, bv_low,bv_high))        
    
    
            if type == 4:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and f >= %s and f <= %s;" % (t, bv_low,bv_high))        
    
   
            if type == 5:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage < 7 and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t))        
 
            if type == 6:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and s >= 1 and m >= 1;" % (t))        

            if type == 7:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and s <= 0.5 and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t))        


            if type == 8:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage < 7 and s = 1 and m >= 1 and f >= 0.75" % (t))        

            if type == 9:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage < 7  and m >= 0.5" % (t))        
   
   
   

            
            damaged_bldgs = self.grid_pg.cur.fetchall()

            
            print "PROBIT: Num. of bldgs with damage %s = %s" % (damage_class,len(damaged_bldgs))
            
            #print damaged_bldgs
            depth = []
            phiinv = []
            histogram = []
            for b in damaged_bldgs:
                depth.append([b[1],b[2]])
                histogram.append(b[2])  

                '''
                if b[2] == 5:
                    histogram.append(6)            
                else:
                    histogram.append(b[2])  
                '''
           
            '''
            n, bins, patches = plt.hist(histogram, [1,2,3,4,5,6])
            
            plt.xticks((1.5,2.5, 3.5,4.5,5.5), ('Flood', 'Minor','Moderate','Major','Collapse'), size = 12)
            plt.yticks(size = 12)
           
            plt.xlabel('Damage Level',size=15,weight='bold')
            plt.ylabel('No. of Buildings', size=15,weight='bold' )

            plt.show()
            plt.close()
            '''
            depth.sort()
            
            num_bldgs = len(depth)
            depth_class1 = []
            i = 0
            start = depth[0][0]
            end = depth[0][0] + 1.0
            count7 = 0
            p7 = []
            
            
            
            p = []
            depth_bins = []
            depths_in_bin = []
            damage_in_bin = []
            #bin_size = 8
            count = 0
            count_class = 0
            
            binary_damage = []
            binary_depths = []

            k = 0
            p = 0
            while i < num_bldgs:
                damage = depth[i][1]
                fd = depth[i][0]
                binary_depths.append(fd)
                if damage >= damage_class:
                    k+=1
                    binary_damage.append(1)
                else:
                    p+=1
                    binary_damage.append(0)
                i+=1
            
            #print k,p
            #print damage_class
            
            #np.ones(5)
            #print binary_damage
            #print binary_depths
            #probit_model = Probit(endog=binary_damage,exog=binary_depths)
            #probit_res = probit_model.fit()
            #print probit_model.cdf(0)
            #print probit_res.params[0]
            #print probit_model.predict(probit_res.params[0],5)
            import pysal
            from pysal.spreg.probit import Probit
            if weight == True:
                #Weight the 0m depth and 0 damage part of the curve
                #This ensures that the curve passes through 0 probability at 0m water depth
                                
                depths_add = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                binary_add = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                binary_depths = depths_add+depths_add+binary_depths
                binary_damage = binary_add+binary_add+binary_damage
                
            x = np.array(binary_depths, ndmin=2).T
            y = np.array(binary_damage, ndmin=2).T
            
            #print binary_damage

            model = Probit(y,x)
            betas = np.around(model.betas, decimals=6)
            #print betas

            #The CDF for PROBIT
            def cdf(x):    
                Prob = 0.5*(1 + erf(x/sqrt(2)))
                return Prob   

            '''
            dep = linspace(0,8,1000)
            dep_beta = betas[1][0]*dep + betas[0][0]
            
            
            plt.plot(dep,cdf(dep_beta))
            
            plt.show()
            plt.close()
            '''
            

            
            b0 = betas[0][0]
            b1 = betas[1][0]
            
            #plt.plot(pinv, depth_bins, 'o', label='Original data', markersize=10)
            #x = linspace(-3,3,51)
            x = linspace(0,5,1000)
            
            i=0
            while i< len(binary_damage):
                
                
                if binary_damage[i] == 0:
                    binary_damage[i] = 0.02*100
                else:
                    binary_damage[i] = 0.98*100
                i+=1
            
            '''
            plt.plot(binary_depths, binary_damage,'o',markersize=10,color='b')
            plt.plot(x,(self.cdf_probit(x,b0,b1))*100,color='r',linewidth=3)
            
            plt.xticks(size = 12,weight='bold')
            plt.yticks((0,50, 100), ('0', '50%',  '100%'), size = 12,weight='bold')
            
            plt.xlabel('Inundation Depth (m)',size=15,weight='bold')
            plt.ylabel('Damage Probability (%)', size=15,weight='bold' )
            plt.ylim([0,100])
            plt.xlim([0,5])
            
            #plt.legend()
            plt.show()
            
            plt.close()
            '''
         
            return b0,b1,binary_depths,binary_damage



    def create_fragility_porter(self,damage_class, bv_low = -1, bv_high = 1,damage_low = 1, type=1,bin_size = 8):

        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename
    
            damage_level = 7
            bin_number = 20
            if type == 1:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage <= %s and bv >= %s and bv < %s;" % (t,damage_class, bv_low,bv_high))        
            
            if type == 2:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage <= %s and s >= %s and s <= %s;" % (t,damage_class, bv_low,bv_high))        
    
            if type == 3:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage <= %s and m >= %s and m <= %s;" % (t,damage_class, bv_low,bv_high))        
    
    
            if type == 4:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage <= %s and f >= %s and f <= %s;" % (t,damage_class, bv_low,bv_high))        
    
    
            
            damaged_bldgs = self.grid_pg.cur.fetchall()
            
            print "Number of buildings with damage %s = %s" % (damage_class,len(damaged_bldgs))
            
            
            depth = []
            phiinv = []
            for b in damaged_bldgs:
                depth.append([b[1],b[2]])
            
            
            depth.sort()
            
            num_bldgs = len(depth)
            depth_class1 = []
            i = 0
            start = depth[0][0]
            end = depth[0][0] + 1.0
            count7 = 0
            p7 = []
            
            
            
            p = []
            depth_bins = []
            depths_in_bin = []
            damage_in_bin = []
            #bin_size = 8
            count = 0
            count_class = 0
            
            #according to Porter 2007 - bin size should be SQRT(num_bldgs) rounded up to the nearest integer value
            bin_size = int(sqrt(num_bldgs))
            if np.mod(float(sqrt(float(num_bldgs))),1) > 0:
                bin_size = bin_size + 1
            
            while i < num_bldgs:
                damage = depth[i][1]
                fd = depth[i][0]
                count += 1
                depths_in_bin.append(fd)
                damage_in_bin.append(damage)
                
                if damage <= damage_class and damage >= damage_low:
                    count_class += 1
                
                if count == bin_size:
                     
                    print damage_in_bin
                    prob = float(count_class)/float(count)
                    
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep         
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    
                    if prob == 1:     
                        #p.append(float(count_class)/float(count+1))
                        #p.append(0.999)
                        CAM=1
                    
                    elif prob == 0:
                        CAM = 1
                    else:
                        depth_bins.append(log(average_depth_bin)) 
                        p.append(prob)
                        
                        
                    depths_in_bin = []
                    damage_in_bin = []
                    
                    count = 0
                    count_class = 0
                
                i+=1

            def phiinv(p):    
                y = sqrt(2)*erfinv(2*p - 1)
                return y
            
            def fragility(x,mu,sig):    
                Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
                return Prob    
            
            pinv = []
            
            for prob in p:
                pinv.append(phiinv(prob))
            
            print p
            print pinv

            x_bar = 0.0
            for dep in depth_bins:
                x_bar = x_bar + log(dep)
            x_bar = x_bar / float(len(depth_bins))
            
            y_bar = 0.0
            for prob in pinv:
                y_bar = y_bar + prob
            y_bar = y_bar / float(len(pinv))
                

            print x_bar
            print y_bar


            sig,mu =  polyfit(pinv,depth_bins,1,full=True)[0]
            print num_bldgs,bin_size,sig,mu
            
            print exp(sig),exp(mu)
    
    
            
            return sig,mu


    def create_fragility_suppasri(self,damage_class, bv_low = -1, bv_high = 1,damage_low = 1, type=1,bin_size = 8):

        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename
    
            damage_level = 7
            bin_number = 20
            if type == 1:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and bv >= %s and bv < %s;" % (t, bv_low,bv_high))        
            
            if type == 2:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and s >= %s and s <= %s;" % (t, bv_low,bv_high))        
    
            if type == 3:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and m >= %s and m <= %s;" % (t, bv_low,bv_high))        
    
    
            if type == 4:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and f >= %s and f <= %s;" % (t, bv_low,bv_high))        
    
   
            if type == 5:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage < 7 and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t))        
 
            if type == 6:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and s >= 1 and m >= 1;" % (t))        

            if type == 7:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and s <= 0.5 and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t))        


            if type == 8:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage < 7 and s = 1 and m >= 1 and f >= 0.75" % (t))        

            if type == 9:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= 1 and damage < 7  and m >= 0.5" % (t))        

      

            
            damaged_bldgs = self.grid_pg.cur.fetchall()
            
            print "SUPPASRI: Num. of bldgs. with damage %s = %s" % (damage_class,len(damaged_bldgs))
            
            
            depth = []
            phiinv = []
            for b in damaged_bldgs:
                depth.append([b[1],b[2]])
            
            
            depth.sort()
    
            num_bldgs = len(depth)
            depth_class1 = []
            i = 0
            start = depth[0][0]
            end = depth[0][0] + 1.0
            count7 = 0
            p7 = []
            
            
            
            p = []
            depth_bins = []
            depths_in_bin = []
            damage_in_bin = []
            #bin_size = 8
            count = 0
            count_class = 0

            #according to Porter 2007 - bin size should be SQRT(num_bldgs) rounded up to the nearest integer value
            bin_size = int(sqrt(num_bldgs))
            if np.mod(float(sqrt(float(num_bldgs))),1) > 0:
                bin_size = bin_size + 1
            j = 0
            start = 0
            
            print "SUPPASRI: Bin Size = %s" % bin_size
            while i < num_bldgs:
                damage = depth[i][1]
                fd = depth[i][0]
                count += 1
                depths_in_bin.append(fd)
                damage_in_bin.append(damage)
                
                if damage >= damage_class:
                    count_class += 1
                
                if count == bin_size:
                    
                    j += 1
                    '''
                    prob = float(count_class)/float(count)
                    if prob > 0 and prob < 1:
                        average_depth_bin = 0
                        for dep in depths_in_bin:
                            average_depth_bin = average_depth_bin + dep    
                        average_depth_bin = average_depth_bin/len(depths_in_bin)
                        depth_bins.append(average_depth_bin)             
                        p.append(float(count_class)/float(count))
                    '''
                    prob = float(count_class)/float(count)
                    start = i+1
                    if prob == 0 and j == 1:
                    #if prob == 0:
                        average_depth_bin = 0
                        for dep in depths_in_bin:
                            average_depth_bin = average_depth_bin + dep    
                        average_depth_bin = average_depth_bin/len(depths_in_bin)
                        depth_bins.append(average_depth_bin)             
                        p.append(0.01)
                                            
                    elif prob == 1:
                        prob = 0.99
                        average_depth_bin = 0
                        for dep in depths_in_bin:
                            average_depth_bin = average_depth_bin + dep    
                        average_depth_bin = average_depth_bin/len(depths_in_bin)
                        depth_bins.append(average_depth_bin)             
                        p.append(0.99)
                    
                    elif prob < 1 and prob > 0:
                        average_depth_bin = 0
                        for dep in depths_in_bin:
                            average_depth_bin = average_depth_bin + dep    
                        average_depth_bin = average_depth_bin/len(depths_in_bin)
                        depth_bins.append(average_depth_bin)             
                        p.append(float(count_class)/float(count))                   
                    
                    depths_in_bin = []
                    damage_in_bin = []
                    count = 0
                    count_class = 0
                
                i+=1
            #Set probability for the remaining bin
            if num_bldgs - start >= 5:
                i = start
                count = 0
                count_class = 0
                depths_in_bin = []
                damage_in_bin = []
                while i < num_bldgs:
                    damage = depth[i][1]
                    fd = depth[i][0]
                    count += 1
                    depths_in_bin.append(fd)
                    damage_in_bin.append(damage)

                    if damage >= damage_class:
                        count_class += 1
                    i+=1
                        
                prob = float(count_class)/float(count)
                if prob == 0 and j == 1:
                #if prob == 0:
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep    
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)             
                    p.append(0.01)
                                        
                elif prob == 1:
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep    
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)             
                    p.append(0.99)
                
                elif prob < 1 and prob > 0:
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep    
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)             
                    p.append(float(count_class)/float(count))  
                    
                    
            

            
            '''
            while i < num_bldgs:
                damage = depth[i][1]
                fd = depth[i][0]
                
                if fd >= start and fd < end:
                    depth_class1.append(damage)
                    count += 1
                    if damage == 7:
                        count7+=1
                    depths_in_bin.append(fd)
                
                elif fd >= end:
                    print "Number in class = %s" % len(depth_class1)
                    print depth_class1
#                    n, bins, patches = plt.hist(depth_class1, [1,2,3,4,5,6,7])
#                    plt.show()
#                    plt.close()
                

                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep
                        
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)
                    
                    print "average depth in bin = %s" % average_depth_bin
                    
                    p7.append(float(count7)/float(count))

                    count = 0
                    count7 = 0
                
                    depths_in_bin = []
                    depth_class1 = []
                    depth_class1.append(damage)

                    count += 1
                    if damage == 7:
                        count7+=1
                    depths_in_bin.append(fd)

                    start = end
                    end = start + 0.25
                    
                i+=1
            '''
                
            
            
            
#            print "Number in class = %s" % len(depth_class1)
#            print depth_class1
#
#            average_depth_bin = 0
#            for dep in depths_in_bin:
#                average_depth_bin = average_depth_bin + dep
#                
#            average_depth_bin = average_depth_bin/len(depths_in_bin)
#            depth_bins.append(average_depth_bin)
#            
#            p7.append(float(count7)/float(count))
#
#
#            print depth_bins
#            print p7
#            


            def phiinv(p):    
                y = sqrt(2)*erfinv(2*p - 1)
                return y
            
            def fragility(x,mu,sig):    
                Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
                return Prob    
            
            pinv = []

            p100 = []
            for prob in p:
                pinv.append(phiinv(prob))
                p100.append(prob*100)
            sig,mu =  polyfit(pinv,depth_bins,1,full=True)[0]


            '''


            i = bin_number -1
            j = 0
            bin_edges = []
            bin_edges.append(depth[0][0])
            bins_damage_23 = []
            total_count = 0
            while i < num_bldgs:
                edge = depth[i][0]
                bin_edges.append(edge)
            
                count = 0
                while j <= i:
            
                    if depth[j][1] == damage_level:
                        count += 1
                        total_count+=1
                    j+=1
                #bins_damage_23.append(float(count)/float(bin_number))
                bins_damage_23.append(float(count))
          
                i+=bin_number
            
            print bins_damage_23
            

            c23 = []
            c_cumul = []
            
            for b in bins_damage_23:
                c23.append(float(b)/float(bin_number))

            
#            r = 0
#            for b in bins_damage_23:
#                c_cumul.append(r+b)
#                r = r+b
#            
#            for b in c_cumul:
#                #c23.append(float(b)/float(bin_number))
#                c23.append(float(b)/float(r))
#            
            print c23
            def phiinv(p):    
                y = sqrt(2)*erfinv(2*p - 1)
                return y
            
            def fragility(x,mu,sig):    
                Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
                return Prob    
            
            print c23

            pinv23 = []
            for p in c23:
                pinv23.append(phiinv(p))
                        
            pinv23_NEW = []
            bins_depth = []
            
            i = 0
            while i < len(pinv23):
                p = pinv23[i]
                
                b = (bin_edges[i+1] + bin_edges[i]) / 2
                
                if p > -50 and p < 50:
                    pinv23_NEW.append(p)
                    bins_depth.append(b)
                i+=1
            
            print pinv23_NEW
            print bins_depth
            
            sig,mu =  polyfit(pinv23_NEW,bins_depth,1,full=True)[0]
            print sig,mu
            
            
            plt.plot(pinv23_NEW, bins_depth, 'o', label='Original data', markersize=10)
            x = linspace(-2,2,51)
            plt.plot(x, sig*x + mu, 'r', label='Fitted line')
            plt.legend()
            plt.show()
            
            plt.close()
            
            x = linspace(0,5,100)
            plt.plot(x,fragility(x,mu,sig),label='Minor',linestyle ='--',color='b')

            
            

            #Fragility from Suprassi
            
            sig4 = 1.0159        #Complete
            mu4 = 4.2243
            
            sig3 = 0.8516        #Major
            mu3 = 3.8458
            
            sig2 = 0.6777        #Moderate
            mu2 = 2.9028
            
            sig1 = 0.6409        #Minor
            mu1 = 2.4409
            
            
            plt.plot(x,fragility(x,mu1,sig1),label='SUP - Minor',linestyle ='--',color='r')
            plt.plot(x,fragility(x,mu2,sig2),label='SUP - Moderate',linestyle ='-.',color='r')
            plt.plot(x,fragility(x,mu3,sig3),label='SUP - Major',linestyle =':',color='r')
            plt.plot(x,fragility(x,mu4,sig4),label='SUP - Complete',color='r')
            
            plt.legend()

            
            plt.show()

            plt.close()

            '''

            
            
            bin = depth_bins
            prob = p100
            
            
            print bin
            print prob
            print sig, mu
            return sig,mu,bin,prob

    def create_fragility_cumfreq2(self,damage_class, bv_low = -1, bv_high = 1,damage_low = 1, type=1,bin_size = 8):
        
        

        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename
    
            if type == 1:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and bv >= %s and bv < %s;" % (t,damage_class, bv_low,bv_high))        
            
            if type == 2:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s >= %s and s <= %s;" % (t,damage_class, bv_low,bv_high))        
    
            if type == 3:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and m >= %s and m <= %s;" % (t,damage_class, bv_low,bv_high))        
    
    
            if type == 4:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and f >= %s and f <= %s;" % (t,damage_class, bv_low,bv_high))        
    
   
            if type == 5:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t,damage_class))        
 
            if type == 6:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s >= 1 and m >= 1;" % (t,damage_class))        

            if type == 7:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s <= 0.5 and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t,damage_class))        


            if type == 8:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s = 1 and m >= 1 and f >= 0.75" % (t,damage_class))        

            if type == 9:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and m >= 0.5" % (t,damage_class))        


            damaged_bldgs = self.grid_pg.cur.fetchall()
            
            print "CUMULATIVE: Num bldgs. with damage %s = %s" % (damage_class,len(damaged_bldgs))

            depth = []
            phiinv = []
            for b in damaged_bldgs:
                depth.append(b[1])
            
            
            depth.sort()
    
            num_bldgs = len(depth)
            
            def phiinv(p):    
                y = sqrt(2)*erfinv(2*p - 1)
                return y

            p = []
            p100 = []
            depths_in_bin = []
            n_bins = []
            bin_avg = []
            bin_avgLN = []
            pinv = []
            depth_avg = 0
            freq = 0
            i = 0
            count = 0

            #according to Porter 2007 - bin size should be SQRT(num_bldgs) rounded up to the nearest integer value
            bin_size = int(sqrt(num_bldgs))
            if np.mod(float(sqrt(float(num_bldgs))),1) > 0:
                bin_size = bin_size + 1
            
            print "CUMFREQ_2: Bin Size = %s" % bin_size
            while i < num_bldgs:
                
                count += 1
                depths_in_bin = depth[i]
                depth_avg = depth_avg + depth[i]
                i+=1
                
                if count == bin_size: 
                    bin_avg.append(float(depth_avg)/float(bin_size))
                    bin_avgLN.append(log(float(depth_avg)/float(bin_size)))
                    freq = freq + float(bin_size)/float(num_bldgs)
                    pinv.append(phiinv(freq))
                    p100.append(freq*100)
                    p.append(freq)
                    count = 0
                    depth_avg = 0
                    
 
            sig,mu =  polyfit(pinv,bin_avgLN,1,full=True)[0]    
#            plt.plot(pinv, bin_avg, 'o', label='Original data', markersize=10)
#            x = linspace(-5,5,51)
#            plt.plot(x, sig*x + mu, 'r', label='Fitted line')
#            plt.legend()
#            plt.show()        
#            plt.close()  
            
                    
            def fragilityln(x,mu,sig):    
                """
                Fragility with depth logged (i.e. ln(x))
                """ 
                Prob = 0.5*(1 + erf((log(x) - mu)/(sqrt(2)*sig)))
                return Prob 
			

#            x = linspace(0,10,1000)
#            plt.plot(bin_avg, p100,'o',markersize=10,color='b')
#            plt.plot(x,(fragilityln(x,mu,sig))*100,color='b',linewidth=3)
#
#            plt.xticks(size = 12,weight='bold')
#            plt.yticks((0,10,20,30,40,50,60,70,80,90,100), size = 12,weight='bold')
#            plt.xlabel('Inundation Depth (m)',size=15,weight='bold')
#            plt.ylabel('Damage Probability (%)', size=15,weight='bold' )
#            plt.ylim([0,100])
#            plt.xlim([0,10])
#            plt.show()
#            plt.close()

            return sig,mu,p100,bin_avg
            

    def create_fragility_cumfreq(self,damage_class, bv_low = -1, bv_high = 1,damage_low = 1, type=1,bin_size = 8):

        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename
    

            if type == 1:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and bv >= %s and bv < %s;" % (t,damage_class, bv_low,bv_high))        
            
            if type == 2:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s >= %s and s <= %s;" % (t,damage_class, bv_low,bv_high))        
    
            if type == 3:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and m >= %s and m <= %s;" % (t,damage_class, bv_low,bv_high))        
    
    
            if type == 4:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and f >= %s and f <= %s;" % (t,damage_class, bv_low,bv_high))        
    
   
            if type == 5:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t,damage_class))        
 
            if type == 6:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s >= 1 and m >= 1;" % (t,damage_class))        

            if type == 7:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s <= 0.5 and m <= 1 and m >= 0.75 and f <= 0.25 and s <= 0.75;" % (t,damage_class))        


            if type == 8:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and s = 1 and m >= 1 and f >= 0.75" % (t,damage_class))        

            if type == 9:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and m >= 0.5" % (t,damage_class))        


            damaged_bldgs = self.grid_pg.cur.fetchall()
            
            print "CUMULATIVE: Num bldgs. with damage %s = %s" % (damage_class,len(damaged_bldgs))

            depth = []
            phiinv = []
            for b in damaged_bldgs:
                depth.append(b[1])
            
            
            depth.sort()
    
            num_bldgs = len(depth)
            depth_class1 = []
            i = 0
            start = depth[0]
            end = depth[0] + 1.0
            count7 = 0
            p7 = []

            def phiinv(p):    
                y = sqrt(2)*erfinv(2*p - 1)
                return y
            
            def fragility(x,mu,sig):    
                Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
                return Prob   
            
            
            p = []
            depth_bins = []
            depths_in_bin = []
            damage_in_bin = []
            count = 0
            count_class = 0

            bin_size = int(sqrt(num_bldgs))
            if np.mod(float(sqrt(float(num_bldgs))),1) > 0:
                bin_size = bin_size + 1
            j = 0
            start = 0
            
            print "BIN size = %s" % bin_size
            
            cumfreq = []
            cumfreq.append([0.01,0])
           
            p = []
            depth_bins = []
            freq = 0
            bins_out = []
            bins_outLN = []
            pinv = []
            p100 = []

            for d in depth:
                freq = freq + 1.0/float(num_bldgs)
                if freq < 1:
                    cumfreq.append([d, freq])
                    depth_bins.append(d)
                    p.append(freq)
                    bins_out.append(d)
                    bins_outLN.append(log(d))
                    pinv.append(phiinv(freq))
                    p100.append(freq*100)
                                
            
            
            bin_size = 0.25
            i = 1
            range_min = float(0)
            range_max = float(10)
            current = range_min
            bin_edges = []
            bin_mid = []
            bin_edges.append(range_min)
            bin_mid.append(0.0)

            while current < range_max:
                current = range_min+i*bin_size
                bin_edges.append(current)
                mid = current
                bin_mid.append(mid)
                i+=1
                
            
            n_bins, bins, patches = plt.hist(depth, bins = bin_edges)
            plt.close()
            p_BIN = []
            freq = 0
            p100_BIN = []
            pinv_BIN = []
            previous = 0
            bin_mid2 = []
            bin_midLN = []
            i=0
            bin_mid = 0

            
            for n in n_bins:                
                freq = previous + float(n)/float(num_bldgs)
                bin_mid = (bins[i+1] + bins[i])/2
                print bin_mid, n
                if freq > 0 and freq <= 1.0:
                    if bin_mid > 2 and bin_mid < 4.875:
                        if freq == 1.0:
                            freq = 0.999   
                            p_BIN.append(freq)
                            p100_BIN.append(freq*100)
                            pinv_BIN.append(phiinv(freq))
                            bin_mid2.append(bin_mid)
                            bin_midLN.append(log(bin_mid))
                            previous = 10.0
                        else:
                            p_BIN.append(freq)
                            p100_BIN.append(freq*100)
                            pinv_BIN.append(phiinv(freq))
                            bin_mid2.append(bin_mid)
                            bin_midLN.append(log(bin_mid))
                            previous = freq
                        

                i+=1

            #sig,mu =  polyfit(pinv,bins_outLN,1,full=True)[0]    
            sig,mu =  polyfit(pinv_BIN,bin_midLN,1,full=True)[0]    

#            plt.plot(pinv, bins_outLN, 'o', label='Original data', markersize=10)
#            x = linspace(-5,5,51)
#            plt.plot(x, sig*x + mu, 'r', label='Fitted line')
#            plt.legend()
#            plt.show()        
#            plt.close()  
            

            def fragilityln(x,mu,sig):    
                """
                Fragility with depth logged (i.e. ln(x))
                """ 
                Prob = 0.5*(1 + erf((log(x) - mu)/(sqrt(2)*sig)))
                return Prob 
            
            '''
            x = linspace(0,10,1000)
            #plt.plot(bins_out, p100,'o',markersize=10,color='b')
            #plt.plot(x,(fragilityln(x,mu,sig))*100,color='b',linewidth=3)

            plt.plot(bin_mid2, p100_BIN,'o',markersize=10,color='b')
            plt.plot(x,(fragilityln(x,mu,sig))*100,color='b',linewidth=3)


            plt.xticks(size = 12,weight='bold')
            plt.yticks((0,10,20,30,40,50,60,70,80,90,100), size = 12,weight='bold')
            plt.xlabel('Inundation Depth (m)',size=15,weight='bold')
            plt.ylabel('Damage Probability (%)', size=15,weight='bold' )
            plt.ylim([0,100])
            plt.xlim([0,10])
            plt.show()
            plt.close()

            bin = bins_out
            print size(bin), size(prob), size(bin_mid2)
            '''
            prob = p100

            return sig,mu,bin,prob

    def cumulative_valencia(self,bins,bin_count, num_bldgs):
        def phiinv(p):    
            y = sqrt(2)*erfinv(2*p - 1)
            return y
        
        def fragility(x,mu,sig):    
            Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
            return Prob 

        def fragilityln(x,mu,sig):    
            """
            Fragility with depth logged (i.e. ln(x))
            """ 
            Prob = 0.5*(1 + erf((log(x) - mu)/(sqrt(2)*sig)))
            return Prob    
    
        p = []
        p100 = []
        pinv = []
        freq = 0
        previous = 0
        bins_out = []
        bins_outLN = []
        #bins_out.append(0)
        #p.append(0.01)
        #p100.append(0.01*100)
        #pinv.append(phiinv(0.01))
        i=0
        for n in bin_count:
            freq = previous + float(n)/float(num_bldgs)
            if freq < 0.00000001:
                freq = 0.00000001
            if freq >= 1:
                freq = 0.9999999
            if bins[i] > 2.54 and bins[i] <= 9.4:
                p.append(freq)
                pinv.append(phiinv(freq))
                p100.append(freq*100)
                bins_out.append(bins[i])
                bins_outLN.append(log(bins[i]))
                previous = freq
            i+=1
        


        sig,mu =  polyfit(pinv,bins_outLN,1,full=True)[0]
        print p
        print sig,mu


        plt.plot(pinv, bins_outLN, 'o', label='Original data', markersize=10)
        x = linspace(-5,5,51)
        plt.plot(x, sig*x + mu, 'r', label='Fitted line')
        plt.legend()
        plt.show()        
        plt.close()
        
        x = linspace(0,10,1000)
        plt.plot(bins_out, p100,'o',markersize=10,color='b')
        plt.plot(x,(fragilityln(x,mu,sig))*100,color='b',linewidth=3)
        plt.plot(x,(fragilityln(x,log(5.26),0.30))*100,color='r',linewidth=3)

        plt.xticks(size = 12,weight='bold')
        plt.yticks((0,10,20,30,40,50,60,70,80,90,100), size = 12,weight='bold')
        plt.xlabel('Inundation Depth (m)',size=15,weight='bold')
        plt.ylabel('Damage Probability (%)', size=15,weight='bold' )
        plt.ylim([0,100])
        plt.xlim([0,10])
        plt.show()
        plt.close()
        return p,bins

    def create_fraglility_curve(self,damage_level, type = 1, bv_low = -1, bv_high = 1, damage_low = 1, damage_high = 7):

        '''
        '''
            
        
        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename
            

            '''
            self.grid_pg.cur.execute("SELECT bv from %s WHERE damage > 0 ;" % (t))        

            bv_all = self.grid_pg.cur.fetchall()

            bv_hist = []
            for b in bv_all:
                bv_hist.append(b[0])
                
            n, bins, patches = plt.hist(bv_hist,21)
            plt.show()
            '''
            if type == 1:
                self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s;" % (t, damage_level))        

            if type == 2:
                if (bv_low >= -1 and bv_high <= 1):
                    self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s and bv > %s and  bv <= %s;" % (t, damage_level, bv_low, bv_high))
                else:
                    print "BV range incorrect - must be within 0 and 1"
         
                    sys.exit()
            
            if type == 3:
                if (bv_low >= -1 and bv_high <= 1):
                    self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage > 0 and damage <= %s and bv > %s and  bv <= %s;" % (t, damage_level, bv_low, bv_high))
                
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()

            if type == 4:
                if (bv_low >= -1 and bv_high <= 1):
                    self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= %s and bv > %s and  bv <= %s;" % (t, damage_level, bv_low, bv_high))
                
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()
            

            if type == 5:
                if (bv_low >= -1 and bv_high <= 1):
                    if (damage_low >= 1 and damage_high <= 7):
                        self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= %s and damage <= %s and bv > %s and  bv <= %s;" % (t, damage_low, damage_high, bv_low, bv_high))
                    else:
                        print "Damage range incorrect - must be within 1 and y"
                        sys.exit()                        
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()

            if type == 6:
                if (bv_low >= -1 and bv_high <= 1):
                    if (damage_low >= 1 and damage_high <= 7):
                        self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage >= %s and damage <= %s and bv > %s and  bv <= %s;" % (t, damage_low, damage_high, bv_low, bv_high))
                    else:
                        print "Damage range incorrect - must be within 1 and y"
                        sys.exit()                        
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()
  

            if type == 7:
                if (bv_low >= -1 and bv_high <= 1):
                    if (damage_low >= 1 and damage_high <= 7):
                        self.grid_pg.cur.execute("SELECT id,fdmax_avg_nodes,damage FROM %s WHERE damage >= %s and damage <= %s and bv > %s and  bv <= %s ;" % (t, damage_low, damage_high, bv_low, bv_high))
                    else:
                        print "Damage range incorrect - must be within 1 and y"
                        sys.exit()                        
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()
    
            
            damaged_bldgs = self.grid_pg.cur.fetchall()

            print "Number of buildings with damage %s = %s" % (damage_level,len(damaged_bldgs))
            depth = []
            phiinv = []
            for b in damaged_bldgs:
                depth.append(b[1])

            #n, bins, patches = plt.hist(depth, 50,normed=1,cumulative=True)
            #plt.show()
            #plt.close()
            
            n, bins, patches = plt.hist(depth, 50,normed=1,cumulative=True)
            plt.close()

            x = []
            i = 0
            while i < (len(bins)-1):
                avg = (bins[i+1]+bins[i])/2
                x.append(avg)
                i+=1

            def phiinv(p):    
                y = sqrt(2)*erfinv(2*p - 1)
                return y
            
            pinv = []
            for p in n:
                pinv.append(phiinv(p))
            
            pinv.pop()
            x.pop()
#            plt.scatter(pinv,x)
#            plt.scatter(x,n)
#            plt.show()
            
            #fit the fragility data
            sig,mu =  polyfit(pinv,x,1)
            
            print sig,mu
 

               
            #A = np.vstack([x, np.ones(len(x))]).T   
                  
            #m, c = np.linalg.lstsq(A, pinv)[0]
            
            #print m,c
            
            '''
            
            xr=polyval([sig,mu],pinv)
            print sig,mu
            plt.scatter(pinv,x)
            plt.plot(pinv,xr)
            plt.show()    
            n, bins, patches = plt.hist(depth, 20,normed=1,cumulative=True)
  
            x1 = []
            i = 0
            while i < (len(bins)-1):
                avg = (bins[i+1]+bins[i])/2
                x1.append(avg)
                i+=1  
  
            x = linspace(0,5,51)
            
            plt.plot(x,self.fragility(x,mu,sig))
            plt.scatter(x1,n)
            
            plt.show()
            '''
            return sig,mu

            
              
        
        '''    

        def I(p,mu,sig):    
            x  = mu + sig*sqrt(2)*erfinv(2*p - 1)
            return x
        
        def I2(phi,mu,sig):    
            x  = mu + sig*sqrt(2)*erfinv(2*p - 1)
            return x



        p = linspace(0, 1, 101)
        plt.plot(p, I(p,2.99,1.12))
        '''


    def create_fraglility_curve_avg(self,damage_level, type = 1, bv_low = -1, bv_high = 1, damage_low = 1, damage_high = 7):

        '''
        '''
            
        
        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename

            '''
            self.grid_pg.cur.execute("SELECT bv from %s WHERE damage > 0 ;" % (t))        

            bv_all = self.grid_pg.cur.fetchall()

            bv_hist = []
            for b in bv_all:
                bv_hist.append(b[0])
                
            n, bins, patches = plt.hist(bv_hist,21)
            plt.show()
            '''
            if type == 1:
                self.grid_pg.cur.execute("SELECT id,fdmax_avg_nodes,damage FROM %s WHERE damage = %s;" % (t, damage_level))        

            if type == 2:
                if (bv_low >= -1 and bv_high <= 1):
                    self.grid_pg.cur.execute("SELECT id,fdmax_avg_nodes,damage FROM %s WHERE damage = %s and bv > %s and  bv <= %s;" % (t, damage_level, bv_low, bv_high))
                else:
                    print "BV range incorrect - must be within 0 and 1"
         
                    sys.exit()
            
            if type == 3:
                if (bv_low >= -1 and bv_high <= 1):
                    self.grid_pg.cur.execute("SELECT id,fdmax_avg_nodes,damage FROM %s WHERE damage > 0 and damage <= %s and bv > %s and  bv <= %s;" % (t, damage_level, bv_low, bv_high))
                
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()

            if type == 4:
                if (bv_low >= -1 and bv_high <= 1):
                    self.grid_pg.cur.execute("SELECT id,fdmax_avg_nodes,damage FROM %s WHERE damage >= %s and bv > %s and  bv <= %s;" % (t, damage_level, bv_low, bv_high))
                
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()
            

            if type == 5:
                if (bv_low >= -1 and bv_high <= 1):
                    if (damage_low >= 1 and damage_high <= 7):
                        self.grid_pg.cur.execute("SELECT id,fdmax_avg_nodes,damage FROM %s WHERE damage >= %s and damage <= %s and bv > %s and  bv <= %s;" % (t, damage_low, damage_high, bv_low, bv_high))
                    else:
                        print "Damage range incorrect - must be within 1 and y"
                        sys.exit()                        
                else:
                    print "BV range incorrect - must be within 0 and 1"
                    sys.exit()
            
  
            damaged_bldgs = self.grid_pg.cur.fetchall()

            print "Number of buildings with damage %s = %s" % (damage_level,len(damaged_bldgs))
            speed = []
            phiinv = []
            for b in damaged_bldgs:
                speed.append(b[1])

            n, bins, patches = plt.hist(speed, 50,normed=1,cumulative=True)
            plt.close()
            x = []
            i = 0
            while i < (len(bins)-1):
                avg = (bins[i+1]+bins[i])/2
                x.append(avg)
                i+=1

            def phiinv(p):    
                y = sqrt(2)*erfinv(2*p - 1)
                return y
            
            pinv = []
            for p in n:
                pinv.append(phiinv(p))
            
            pinv.pop()
            x.pop()

            sig,mu = polyfit(pinv,x,1)
            

            return sig,mu
      

    def add_max_building_values_to_database(self,new_tablename):
        '''
        Add the max values at each building to the PostGIS database for the grid

    
        '''       
        
        self.output_tablename = new_tablename
        
        self.grid_pg.add_building_results_table(new_tablename)
        
        etamaxNodesALL = self.building_nodes_grp.variables['eta_max'][:]            
        speedmaxNodesALL = self.building_nodes_grp.variables['speed_max'][:]    
        
        for building in self.buildings_dict.iteritems():      
            id = building[0]
            nodes = building[1]['nodes']
            sides = building[1]['sides']
            elements = building[1]['elements']
            perimeter = building[1]['perimeter']
            type = building[1]['type']
            
            #get the node output
            etamax = 0
            speedmax = 0
            exposure = 0
    
            speedNodes = []
            etaNodes = []
            flooddepthNodes = []
            elevation = 0
            
            for n in nodes:
                z = self.elevation[n-1]
                
                '''
                if etamax < eta:
                    etamax = eta
                
                '''
                speed = speedmaxNodesALL[self.node_ids_dict[n]]
                fd = etamaxNodesALL[self.node_ids_dict[n]] - z
                eta = etamaxNodesALL[self.node_ids_dict[n]]
                elevation = elevation + z 

               # if speedmax < speed:
               #     speedmax = speed
               
                speedNodes.append(speed)
                etaNodes.append(eta) 
                flooddepthNodes.append(fd)
                #etamax = etamax + (etamaxNodesALL[self.node_ids_dict[n]] - z)
                #speedmax = speedmax + speedmaxNodesALL[self.node_ids_dict[n]]

            '''
            if id == 295:
                speedNodes.sort()
                print speedNodes
                
                nodeEdge = self.grid_pg.get_nodes_at_building_edge(id)
                speedEdge295 = []
                for n in nodeEdge:
                    speed = speedmaxNodesALL[self.node_ids_dict[n]]
                    speedEdge295.append(speed)
                speedEdge295.sort()
                print speedEdge295
                sys.exit()
                    
            '''
            etaNodes.sort()
            etaNodes.reverse()

            speedNodes.sort()
            speedNodes.reverse()    

            flooddepthNodes.sort()
            flooddepthNodes.reverse()
            
            fdMaxAverage = 0
            for max in flooddepthNodes:
                fdMaxAverage = fdMaxAverage + max
            
                    
            if len(nodes) > 0:
                fdMaxAverage = fdMaxAverage/len(nodes)
                
                
                #calculate the average values
#                speedmaxAVG = 0
#                i = 0
#                while i < len(nodes)/2:
#                    speedmaxAVG = speedmaxAVG + speedNodes[i]
#                    i+=1
                    
                #speedmax = speedmaxAVG/(i-1)
                
                #speedmax = speedmax/len(nodes)
                #etamax = etamax/len(nodes)
                etamax = etaNodes[0]
                fdmax = flooddepthNodes[0]
                speedmax = speedNodes[0]
                exposure = (0.6666*etamax + 0.3333*speedmax)
                #mu and sig from supprasi et al
                #Try - Wooden house complete damage
                #mu = 4.2243
                #sig = 1.0159
                #Wooden house - minor damage

                self.grid_pg.cur.execute("SELECT bv, damage from %s WHERE id = %s ;" % (new_tablename,id)) 
                building_return = self.grid_pg.cur.fetchall()

                if building_return[0][1] > 0:        #damage is greater than 0 (i.e. building has been surveyed)

                    bv = building_return[0][0]
                    
                    if bv > -1 and bv <= 0.25:
                        mu = 3.45374067692
                        sig = 0.601532235676    
                    if bv > 0.25 and bv <= 0.75:
                        mu = 3.09628891017
                        sig = 0.576207872969    
                    if bv > 0.75:
                        mu = 2.88133508574
                        sig = 0.662402196488    
                
                    fragility = self.fragility(fdmax,mu,sig)

                    #self.grid_pg.cur.execute("UPDATE %s SET fragility = %s WHERE id = %s;" % (new_tablename,fragility,id))

                #else:
                    #self.grid_pg.cur.execute("UPDATE %s SET fragility = -1 WHERE id = %s;" % (new_tablename,id))

    
                elevation = elevation / len(nodes)
                self.grid_pg.cur.execute("UPDATE %s SET etamax_nodes = %s WHERE id = %s;" % (new_tablename,etamax,id))
                self.grid_pg.cur.execute("UPDATE %s SET fdmax_nodes = %s WHERE id = %s;" % (new_tablename,fdmax,id))

                self.grid_pg.cur.execute("UPDATE %s SET speedmax_nodes = %s WHERE id = %s;" % (new_tablename,speedmax,id))
                self.grid_pg.cur.execute("UPDATE %s SET exposure = %s WHERE id = %s;" % (new_tablename,exposure,id))
                self.grid_pg.cur.execute("UPDATE %s SET z = %s WHERE id = %s;" % (new_tablename,z,id))
                
                #add the etaMaxAverage
                self.grid_pg.cur.execute("UPDATE %s SET fdmax_avg_nodes = %s WHERE id = %s;" % (new_tablename,fdMaxAverage,id))
        
        self.grid_pg.conn.commit()



    
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
                    uWeighted = uWeighted + u*(self.slen[s-1])/perimeter
                    vWeighted = vWeighted + v*(self.slen[s-1])/perimeter
                    
                    #Calculate the velocity components tangential and normal to the side
                    ut = self.sdx[s-1]*u - self.sdy[s-1]*v
                    un = self.sdy[s-1]*u + self.sdx[s-1]*v


            speedAvg = speedAvg/nsides


            speedAvgWeighted2 = math.sqrt(uWeighted*uWeighted+vWeighted*vWeighted)

            speedWeighted2 = uWeighted
            speedMax_dict['sides'].append(speedAvgWeighted2)

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
        startTime = datetime(2011,03,11,14,46,24)
       
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
            print time
            print eta[6]
            print len(eta[6])
            print len(time)
            for t in time:
                currentDateTime = startTime + timedelta(seconds=t)
                datestring = str(currentDateTime)
                datestring = datestring.replace(' ','T')
                i = 0
                while i < number_of_key_points:
                    outfiles[i].write("%s %s\n" % (datestring, str(eta[i][tStep]*100)))
                    i+=1
                tStep+=1
            
            print tStep
            #close all the output files
            for file in outfiles:
                file.close()
                

    def fragility(self, x,mu,sig):	
        """
        Fragility function - Supprasi et al.
        
        
        """ 
        Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
        return Prob    

    
    def fragilityln(self, x,mu,sig):    
        """
        Fragility with depth logged (i.e. ln(x))
        """ 
        Prob = 0.5*(1 + erf((log(x) - mu)/(sqrt(2)*sig)))
        return Prob  


    #The CDF for PROBIT
    def cdf_probit(self,x,b0,b1):    
        
        x = b1*x + b0
        Prob = 0.5*(1 + erf(x/sqrt(2)))
        return Prob   
    

    #Given the fitted PROBIT parameters (b0 and b1), solve for mu (i.e. the standard mean)    
    def get_mu_from_cdf_probit(self,b0,b1):    
        
        #The standard mean of a CDF is were probability = 0.5
        mu = (sqrt(2)*erfinv(0) - b0)/b1
        x = 1
        sig = (x - mu) / ( sqrt(2) * erfinv(erf((b1*x + b0)/sqrt(2))) )

        return sig,mu
    
#    def create_fragility_suppasri_norm(self,damage_class=2,mu1=1,sig1=0.4, mu2=2.5,sig2=0.6, mu3=3.25,sig3=0.6, c1=100, c2=100, c3=100):
    def create_fragility_suppasri_norm(self,damage_class,s1,s2,s3):

        '''
        Create a binned (i.e. method used in suppasri/koshimura) fragility function using a set of
        normal distributed buildings
        
        This will be used as an example to show how the different technique work
        
        This will be compared against other techniques (probit and cumulative) to illustrate the differences
        '''
        
        #CREATE THE DATA
#        s1 = np.random.normal(mu1, sig1, c1).tolist()
#        s2 = np.random.normal(mu2, sig2, c2).tolist()
#        s3 = np.random.normal(mu3, sig3, c3).tolist()

        depth = []

        for s in s1:
            depth.append([s,1])

        for s in s2:
            depth.append([s,2])
            
        for s in s3:
            depth.append([s,3])
                
        depth.sort()
        
        
        
        num_bldgs = len(depth)
        depth_class1 = []
        i = 0
        start = depth[0][0]
        end = depth[0][0] + 1.0
        count7 = 0
        p7 = []
        
        
        
        p = []
        depth_bins = []
        depths_in_bin = []
        damage_in_bin = []
        #bin_size = 8
        count = 0
        count_class = 0
        
        n_out = []              #the specimens in each bin
        bins_out = []         #the edges of the bins
        n_d1 = 0
        n_d2 = 0
        n_d3 = 0
        n1 = []
        n2 = []
        n3 = []
        
        
        #according to Porter 2007 - bin size should be SQRT(num_bldgs) rounded up to the nearest integer value
        bin_size = int(sqrt(num_bldgs))
        if np.mod(float(sqrt(float(num_bldgs))),1) > 0:
            bin_size = bin_size + 1
        j = 0
        start = 0
        
        prob_one_FLAG = False
        
        print "SUPPASRI: Bin Size = %s" % bin_size
        while i < num_bldgs:
            damage = depth[i][1]
            fd = depth[i][0]
            count += 1
            depths_in_bin.append(fd)
            damage_in_bin.append(damage)
            if damage == 1:
                n_d1+=1
            if damage == 2:
                n_d2+=1
            if damage == 3:
                n_d3+=1
                
            
            
            if damage >= damage_class:
                count_class += 1
            
            if count == bin_size:
                
                j += 1
                '''
                prob = float(count_class)/float(count)
                if prob > 0 and prob < 1:
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep    
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)             
                    p.append(float(count_class)/float(count))
                '''
                prob = float(count_class)/float(count)
                start = i+1
                if prob == 0 and j == 1:
                #if prob == 0:
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep    
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)             
                    p.append(0.01)
                                        
                elif prob == 1:
                    prob = 0.99
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep    
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    if prob_one_FLAG == False:
                        depth_bins.append(average_depth_bin)             
                        p.append(0.99)
                        prob_one_FLAG = True
                
                elif prob < 1 and prob > 0:
                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep    
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)             
                    p.append(float(count_class)/float(count))                   
                
                
                bins_out.append(average_depth_bin)
                n_out.append([n_d1,n_d2,n_d3])
                n1.append(n_d1)
                n2.append(n_d2)
                n3.append(n_d3)
        
    
                n_d1 = 0
                n_d2 = 0
                n_d3 = 0
                
                depths_in_bin = []
                damage_in_bin = []
                count = 0
                count_class = 0
            
            i+=1
        #Set probability for the remaining bin
        if num_bldgs - start >= 5:
            i = start
            count = 0
            count_class = 0
            depths_in_bin = []
            damage_in_bin = []
            while i < num_bldgs:
                damage = depth[i][1]
                fd = depth[i][0]
                count += 1
                depths_in_bin.append(fd)
                damage_in_bin.append(damage)
    
                if damage == 1:
                    n_d1+=1
                if damage == 2:
                    n_d2+=1
                if damage == 3:
                    n_d3+=1
                    
                if damage >= damage_class:
                    count_class += 1
                i+=1
                    
            prob = float(count_class)/float(count)
            if prob == 0 and j == 1:
            #if prob == 0:
                average_depth_bin = 0
                for dep in depths_in_bin:
                    average_depth_bin = average_depth_bin + dep    
                average_depth_bin = average_depth_bin/len(depths_in_bin)
                depth_bins.append(average_depth_bin)             
                p.append(0.01)
                                    
            elif prob == 1:
                average_depth_bin = 0
                for dep in depths_in_bin:
                    average_depth_bin = average_depth_bin + dep    
                average_depth_bin = average_depth_bin/len(depths_in_bin)
                if prob_one_FLAG == False:
                    depth_bins.append(average_depth_bin)             
                    p.append(0.99)
                    prob_one_FLAG = True
            
            elif prob < 1 and prob > 0:
                average_depth_bin = 0
                for dep in depths_in_bin:
                    average_depth_bin = average_depth_bin + dep    
                average_depth_bin = average_depth_bin/len(depths_in_bin)
                depth_bins.append(average_depth_bin)             
                p.append(float(count_class)/float(count))  

            bins_out.append(average_depth_bin)
            n_out.append([n_d1,n_d2,n_d3])
            n1.append(n_d1)
            n2.append(n_d2)
            n3.append(n_d3)
        
        def phiinv(p):    
            y = sqrt(2)*erfinv(2*p - 1)
            return y
        
        def fragility(x,mu,sig):    
            Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
            return Prob    
        
        pinv = []
        
        p100 = []
        for prob in p:
            pinv.append(phiinv(prob))
            p100.append(prob*100)
        sig,mu =  polyfit(pinv,depth_bins,1,full=True)[0]
        
        print p 
        print depth_bins
        
        bin = depth_bins
        prob = p100
        
        
        
        i = len(prob)
        while i<len(bins_out):
            prob.append(100)
            i+=1

                
        
        
        #print n_out
        #print depth_bins
        
        return sig,mu,bins_out,prob,n1,n2,n3


 #   def create_fragility_probit_norm(self, damage_class,mu1=1,sig1=0.4, mu2=2.5,sig2=0.6, mu3=3.25,sig3=0.6, c1=100, c2=100, c3=100, weight=False):
    def create_fragility_probit_norm(self, damage_class,s1,s2,s3,weight=False):
  
        '''
        Create a probit (i.e. method used in reese) fragility function using a set of
        normal distributed buildings
        
        This will be used as an example to show how the different technique work
        
        This will be compared against other techniques (binned and cumulative) to illustrate the differences
        '''
                
        #CREATE THE DATA
#        s1 = np.random.normal(mu1, sig1, c1).tolist()
#        s2 = np.random.normal(mu2, sig2, c2).tolist()
#        s3 = np.random.normal(mu3, sig3, c3).tolist()

        depth = []

        for s in s1:
            depth.append([s,1])

        for s in s2:
            depth.append([s,2])
            
        for s in s3:
            depth.append([s,3])
                
        depth.sort()

        
        #print damaged_bldgs
        phiinv = []
        histogram = []

        num_bldgs = len(depth)
        depth_class1 = []
        i = 0
        start = depth[0][0]
        end = depth[0][0] + 1.0
        count7 = 0
        p7 = []
        
        
        
        p = []
        depth_bins = []
        depths_in_bin = []
        damage_in_bin = []
        #bin_size = 8
        count = 0
        count_class = 0
        
        binary_damage = []
        binary_depths = []

        k = 0
        p = 0
        while i < num_bldgs:
            damage = depth[i][1]
            fd = depth[i][0]
            binary_depths.append(fd)
            if damage >= damage_class:
                k+=1
                binary_damage.append(1)
            else:
                p+=1
                binary_damage.append(0)
            i+=1
        
        #print k,p
        #print damage_class
        
        #np.ones(5)
        #print binary_damage
        #print binary_depths
        #probit_model = Probit(endog=binary_damage,exog=binary_depths)
        #probit_res = probit_model.fit()
        #print probit_model.cdf(0)
        #print probit_res.params[0]
        #print probit_model.predict(probit_res.params[0],5)
        import pysal
        from pysal.spreg.probit import Probit
        if weight == True:
            #Weight the 0m depth and 0 damage part of the curve
            #This ensures that the curve passes through 0 probability at 0m water depth
                            
            depths_add = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            binary_add = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            binary_depths = depths_add+depths_add+binary_depths
            binary_damage = binary_add+binary_add+binary_damage
            
        x = np.array(binary_depths, ndmin=2).T
        y = np.array(binary_damage, ndmin=2).T
        
        #print binary_damage

        model = Probit(y,x)
        betas = np.around(model.betas, decimals=6)
        #print betas

        #The CDF for PROBIT
        def cdf(x):    
            Prob = 0.5*(1 + erf(x/sqrt(2)))
            return Prob   

        
        b0 = betas[0][0]
        b1 = betas[1][0]
        
        #plt.plot(pinv, depth_bins, 'o', label='Original data', markersize=10)
        #x = linspace(-3,3,51)
        x = linspace(0,5,1000)
        
        i=0
        while i< len(binary_damage):
            
            
            if binary_damage[i] == 0:
                binary_damage[i] = 0.02*100
            else:
                binary_damage[i] = 0.98*100
            i+=1
        
        return b0,b1,binary_depths,binary_damage


#    def create_fragility_cumfreq2_norm(self, damage_class, mu1=1,sig1=0.4, mu2=2.5,sig2=0.6, mu3=3.25,sig3=0.6, c1=100, c2=100, c3=100):
    def create_fragility_cumfreq2_norm(self, damage_class,s1,s2,s3 ):
        
        '''
        Create a cumulative frequency (i.e. method used in valencia) fragility function using a set of
        normal distributed buildings
        
        This will be used as an example to show how the different technique work
        
        This will be compared against other techniques (probit and cumulative) to illustrate the differences
        '''
        
        
        #Damage Levels and associated mu and sig distribution parameters
#        txt1 =  "D1"
#        mu1 = 1
#        sig1 = 0.4 
#        
#        
#        txt2 =  "D2"
#        mu2 = 2.5
#        sig2 = 0.6
#        
#        txt3 =  "D3"
#        mu3 = 3.25
#        sig3 = 0.6
        
        #CREATE THE DATA
#        s1 = np.random.normal(mu1, sig1, c1).tolist()
#        s2 = np.random.normal(mu2, sig2, c2).tolist()
#        s3 = np.random.normal(mu3, sig3, c3).tolist()

        depth = []
               
        if damage_class == 1:
            depth = s1
        elif damage_class == 2:
            depth = s2
        elif damage_class == 3:
            depth = s3
        else:
            sys.exit()

        depth.sort()
        
        #count, bins, ignored = plt.hist(depth, cumulative=True)
        #plt.show()
        #plt.close()

 
        phiinv = []
        
        
        num_bldgs = len(depth)
        print "NUMBER OF BUILDINGS = %s" % num_bldgs
        def phiinv(p):    
            y = sqrt(2)*erfinv(2*p - 1)
            return y
        
        p = []
        p100 = []
        depths_in_bin = []
        n_bins = []
        bin_avg = []
        bin_avgLN = []
        pinv = []
        depth_avg = 0
        freq = 0
        i = 0
        count = 0
        
        #according to Porter 2007 - bin size should be SQRT(num_bldgs) rounded up to the nearest integer value
        bin_size = int(sqrt(num_bldgs))
        print "bins size 1 = %s" % bin_size
        print "MOD thing = %s" % np.mod(float(sqrt(float(num_bldgs))),1)
        if np.mod(float(sqrt(float(num_bldgs))),1) > 0:
            print "INSIDE!!!"
            bin_size = bin_size + 1
        
        print "CUMFREQ_2: Bin Size = %s" % bin_size
        while i < num_bldgs:
            
            count += 1
            depths_in_bin = depth[i]
            depth_avg = depth_avg + depth[i]
            i+=1
            
            if count == bin_size: 
                bin_avg.append(float(depth_avg)/float(bin_size))
                bin_avgLN.append(log(float(depth_avg)/float(bin_size)))
                freq = freq + float(bin_size)/float(num_bldgs)
                pinv.append(phiinv(freq))
                p100.append(freq*100)
                p.append(freq)
                count = 0
                depth_avg = 0
                
        
        sig,mu =  polyfit(pinv,bin_avgLN,1,full=True)[0]    


        plt.plot(pinv, bin_avgLN, 'o', label='Original data', markersize=10)
        x = linspace(-10,10,200)
        plt.plot(x, sig*x + mu, 'r', label='Fitted line')
        plt.legend()
        plt.show()        
        plt.close() 
                
        def fragilityln(x,mu,sig):    
            """
            Fragility with depth logged (i.e. ln(x))
            """ 
            Prob = 0.5*(1 + erf((log(x) - mu)/(sqrt(2)*sig)))
            return Prob 
        
        return sig,mu,bin_avg,p100

    #def create_fragility_cumfreq_no_bin_norm(self, damage_class, mu1=1,sig1=0.4, mu2=2.5,sig2=0.6, mu3=3.25,sig3=0.6, c1=100, c2=100, c3=100):
        
    def create_fragility_cumfreq_no_bin_norm(self, damage_class, s1,s2,s3):
        '''
        Create a cumulative frequency (i.e. method used in valencia) fragility function using a set of
        normal distributed buildings
        
        This will be used as an example to show how the different technique work
        
        This will be compared against other techniques (probit and cumulative) to illustrate the differences
        '''
        
        
        #Damage Levels and associated mu and sig distribution parameters
#        txt1 =  "D1"
#        mu1 = 1
#        sig1 = 0.4 
#        
#        
#        txt2 =  "D2"
#        mu2 = 2.5
#        sig2 = 0.6
#        
#        txt3 =  "D3"
#        mu3 = 3.25
#        sig3 = 0.6
        
        #CREATE THE DATA
#        s1 = np.random.normal(mu1, sig1, c1).tolist()
#        s2 = np.random.normal(mu2, sig2, c2).tolist()
#        s3 = np.random.normal(mu3, sig3, c3).tolist()

        depth = []
               
        if damage_class == 1:
            depth = s1
        elif damage_class == 2:
            depth = s2
        elif damage_class == 3:
            depth = s3
        else:
            sys.exit()

        depth.sort()
        
        #count, bins, ignored = plt.hist(depth, cumulative=True)
        #plt.show()
        #plt.close()

 
        phiinv = []
        
        
        num_bldgs = len(depth)
        print "NUMBER OF BUILDINGS = %s" % num_bldgs
        def phiinv(p):    
            y = sqrt(2)*erfinv(2*p - 1)
            return y
        
        p = []
        p100 = []
        depths_in_bin = []
        n_bins = []
        bin_avg = []
        bin_avgLN = []
        pinv = []
        depth_avg = 0
        freq = 0
        i = 0
        count = 0
        depth_ln = []
        depth_out = []
        
        while i < num_bldgs:
            
            freq = float(i+1)/float(num_bldgs)
            
            if freq < 0.98 and freq > 0.02:
                pinv.append(phiinv(freq))
                p100.append(freq*100)
                p.append(freq)
                depth_out.append(depth[i])
                depth_ln.append(log(depth[i]))
            i+=1

        
        
        sig,mu =  polyfit(pinv,depth_ln,1,full=True)[0]    
        
        
#        print sig,mu
#
#        plt.plot(pinv, depth_ln, 'o', label='Original data', markersize=10)
#        x = linspace(-10,10,200)
#        plt.plot(x, sig*x + mu, 'r', label='Fitted line')
#        plt.legend()
#        plt.show()        
#        plt.close() 
                
        def fragilityln(x,mu,sig):    
            """
            Fragility with depth logged (i.e. ln(x))
            """ 
            Prob = 0.5*(1 + erf((log(x) - mu)/(sqrt(2)*sig)))
            return Prob 
        
        return sig,mu,depth_out,p100

