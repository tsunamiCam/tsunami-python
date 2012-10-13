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

import numpy as np


from Grid_Class import GridPG



from scipy.special import erf,erfinv
from numpy import linspace, exp, sqrt
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


        
        self.__init_buildings_dictionary__()

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



    def create_fraglility_curve3(self):

        '''
        '''
        
        
        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            
            t = self.output_tablename
            self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage>0" % (t))        
            damaged_bldgs = self.grid_pg.cur.fetchall()

            depth = []
            damaged  = []
            for b in damaged_bldgs:
                depth.append(b[1])
                damaged.append(b[2])
                
            n, bins, patches = plt.hist(damaged, 7)
            plt.show()
  




    def create_fraglility_curve2(self,damage_level,stories):

        '''
        '''
        
        
        if self.output_tablename == "":
            print "Please add output data to the buildings table for this run..."
            return
        else:
            

            
            t = self.output_tablename

            
            

            self.grid_pg.cur.execute("SELECT id,fdmax_nodes,damage FROM %s WHERE damage = %s AND s=%s;" % (t, damage_level,stories))        
            damaged_bldgs = self.grid_pg.cur.fetchall()

            print "Number of buildings with damage %s = %s" % (damage_level,len(damaged_bldgs))
            depth = []
            phiinv = []
            for b in damaged_bldgs:
                depth.append(b[1])
                
                
            n, bins, patches = plt.hist(depth, 20,normed=1,cumulative=True)
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

    def create_fraglility_suppasri(self,damage_class, bv_low = -1, bv_high = 1,damage_low = 1, type=1,bin_size = 8):

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

            while i < num_bldgs:
                damage = depth[i][1]
                fd = depth[i][0]
                count += 1
                depths_in_bin.append(fd)
                damage_in_bin.append(damage)
                
                if damage <= damage_class and damage >= damage_low:
                    count_class += 1
                
                if count == bin_size:
                     
                    prob = float(count_class)/float(count)
                    if prob > 0 and prob < 1:
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

            for prob in p:
                pinv.append(phiinv(prob))
            
            
            sig,mu =  polyfit(pinv,depth_bins,1,full=True)[0]
            print sig,mu


#            plt.plot(pinv, depth_bins, 'o', label='Original data', markersize=10)
#            x = linspace(-3,3,51)
#            plt.plot(x, sig*x + mu, 'r', label='Fitted line')
#            plt.legend()
#            plt.show()
#            
#            plt.close()

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

            
            

            return sig,mu
            

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
                

    def fragility(self, x,mu,sig):	
        """
        Fragility function - Supprasi et al.
        
        
        """ 
        Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
        return Prob    
    
    

"""

def I(p,mu,sig):    
    x  = mu + sig*sqrt(2)*erfinv(2*p - 1)
    return x

def I2(phi,mu,sig):    
    x  = mu + sig*sqrt(2)*erfinv(2*p - 1)
    return x



p = linspace(0, 1, 101)
plt.plot(p, I(p,2.99,1.12))


"""

