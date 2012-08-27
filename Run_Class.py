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


class RunClass:
    
    def __init__(   self,
                    run_id,                         #unique id of the run (INTEGER)
                    name,
                    run_filename,                            #relative directory location
                    run_dir,
                    description, 
                    grid,                           #GridClass object that defines the grid
                    number_of_iterations,           #number of iterations for the run
                    delta_time,                      #Time step in seconds between iterations
                    lat0,
                    long0,
                    courantCheck = 0,				#default - courant Check is OFF
                    omega0 = 0,
                    sea_level_offset = 0,           #offset in meters above or below mean sea level
                    user='tbone'):
        '''
        Initialization function for the RunClass
        '''
        self.run_id = run_id
        self.grid = grid
        self.name = name
        self.description = description
        self.run_filename = run_filename
        self.run_dir = run_dir

        self.grid_filename = self.grid.grid_filename
        self.buildings_dbname = grid.buildings_dbname

        self.run_nc = RunNC(run_dir+"/"+self.run_filename, self.grid.grid_dir+"/"+self.grid.grid_filename,
                            run_id, name, description,grid.grid_dbname, grid.epsg, grid.buildings_dbname,grid.user)      
        
        
        '''
        FRICTION
        -------------------------------------------------
        Get the elements and element_codes in the domain.
        These are to be saved the run file before a model run occurs.
        '''
        
        self.elements = self.grid.get_elements_in_area(0)
        i = 0
        length = len(self.elements)
        while i < length:
            self.elements[i][1] = 1
            i+=1	
        
        self.ntype = 1
        self.friction_parameters = []
        		
		#set default friction parameters for element code = 1
        self.friction_parameters.append([1, 0.0,   0.1,   0.0,   0.0025,   0.0,   0.0,   0.0])

      
        '''
        Form DRAG
        -------------------------------------------------
        All form drag elements drag parameters are saved to FD2dInit.dat in self.run()
        '''  
        self.form_drag_elements = []				#list of form drag elements and corresponding element codes
        self.number_of_drag_types = 0
        self.drag_parameters = []				#Form Drag parameters list - each row in the list corresponds to a different element type
        self.ifdrag = 0								#Form drag type - 0 = no drag
        
        
        #TEST: Alter element adjacency matrix
        #self.ieadj = self.grid.grid_nc.get_grid_ieadj()
        #print "SIZE if ieadj = %s" % len(self.ieadj)

        
		
        self.ricom = Ricom(run_dir)
        self.ricom.title = description
        self.ricom.gridfilename = self.grid.grid_dir+"/"+self.grid.grid_filename        
        self.ricom.run_filename = self.run_dir + "/" + self.run_filename  
        self.ricom.nitn = number_of_iterations
        self.ricom.delt = delta_time
        self.ricom.elev = sea_level_offset
        self.ricom.lat0 = lat0
        self.ricom.long0 = long0
        self.ricom.omega0 = omega0
        
        self.courantCheck = courantCheck
        self.keyPointsAdded = False
        self.buildingsAdded = False
        
        #subdomains
        self.number_of_subdomains = 0
        
    def __del__ (self):
        """
        Class deconstructor
    
        """
        cam = 1 #TODO

    def run(self):
        """
        Run RiCOM
        
        """


		#Write the requested element codes for the run to the Run File
        self.run_nc.set_element_codes(self.elements,self.ntype)
        self.ricom.ntype = self.ntype
        
        #write the element adjacency array if change:
        
        #if self.ieadj != []:
        #    self.run_nc.set_ieadj(self.ieadj)


        #write the drag elements and parameters to FD2dInit.dat
        dragfile = open(self.run_dir + "/FD2dInit.dat","w")
        
        dragfile.write("%s\n" % self.number_of_drag_types)
        for fd in self.drag_parameters:
            dragfile.write("%s %s %s %s\n" % (fd[0],fd[1],fd[2],fd[3]))
        
        dragfile.write("%s\n" % len(self.form_drag_elements))
        for ele in self.form_drag_elements:
            dragfile.write("%s %s\n" % (ele[0],ele[1])) 
        
        dragfile.close()       

        
            
        

        self.write_rcm(self.run_dir + "/tsunami.rcm")
        #self.write_rcm("tsunami.rcm")
        self.ricom.write_fault_file(self.run_dir + "/fault.param")
        self.grid.add_run(self.run_id, self.name, self.description,
                          run_filename = self.run_filename,grid_filename = self.grid_filename)
        
        if (self.courantCheck == 1):
        	ricomLIB = cdll.LoadLibrary("/opt/local/bin/PhD/Ricom/ricom11.5.9_nc_dcCrcheck.dylib")
        else:
        	ricomLIB = cdll.LoadLibrary("/opt/local/bin/PhD/Ricom/ricom12.1.x.nc.omp.dylib") 
        	#ricomLIB = cdll.LoadLibrary("/opt/local/bin/PhD/Ricom/ricom11.5.9_nc.dylib") 
        
        self.run_nc.close()
        
        pwd = os.getcwd()               #change the directory to the Run directory
        os.chdir(self.run_dir)

        ricomLIB.ricom_()
        os.chdir(pwd)

        
        
        
        
        self.run_nc.open('r')               #open output dataset for reading
        #self._write_output_()

        #print self.run_nc.grid_filename
        #print "finished"

    def add_fault_rupture_from_file(self,fault_filename):
        """
        Add fault rupture from a file (.param format)
        """
        
        self.ricom.add_fault_rupture_from_file(fault_filename)


    def add_restart(self,restart_filename):
        """
        Add fault rupture from a file (.param format)
        """
        
        self.ricom.add_restart(restart_filename)
    
    
    def add_friction_from_file(self, fric_filename = "friction.par"):
        """
        Add friction definitions from a file (.par format)
        
        """
        self.ricom.add_friction_from_file(fric_filename)




    def add_building_output_locations(self,area_id,start,end,step,type='BUILDINGS_AS_HOLES'):
        """
        For each of node,element and side corresponding to the buildings in the grid,
        request output
        
        Input:
            area_id - id of the area polygon that you want to get the info inside
            start: iteration # where output should start
            end: iteration # where output should end
            type:   BUILDINGS_AS_HOLES
                    BUILDINGS_AS_POINTS
                    BUILDINGS_GRIDDED
        
        """ 
        print "Getting buildings locations..."
        
        dictionary = self.grid.get_building_output_locations(area_id,type)
        if (dictionary != {}):
            self.run_nc.add_building_output_locations(dictionary, start, end,step)


    def add_building_output_locations2(self,areasList,start,end,step):
        """
        For each of node,element and side corresponding to the buildings in the grid,
        request output
        
        Input:
            areasList - list of areas and corresponding building types:
                    type:   BUILDINGS_AS_HOLES
                            BUILDINGS_AS_POINTS
                            BUILDINGS_GRIDDED
            start: iteration # where output should start
            end: iteration # where output should end
        
        
        """ 
        print "Getting buildings locations..."
        
        dictionaries = []
        dictionary = {}
        
        for a in areasList:
                
            dictionaries.append(self.grid.get_building_output_locations(a[0],a[1]))
       
        for dict in dictionaries:
            for row in dict.iteritems(): 
                dictionary[row[0]] = row[1]              

        print "Number of buildings = %s" % (len(dictionary))

        if (dictionary != {}):
            self.run_nc.add_building_output_locations(dictionary, start, end,step)


    def add_element_output_locations(self, xy, epsgIN,start,end,step):
        """
        Given a list of points [x,y], return a list of ids corresponding to the element that 
        is closest to each point 

        
        """ 
        elementIds = self.grid.get_element_output_locations(xy,epsgIN)
        if(elementIds != []):
            self.run_nc.add_element_output_locations(elementIds,start,end,step)

    def add_node_output_locations(self, xy,epsgIN,start,end,step):
        """
        Given a list of points [x,y], return a list of ids corresponding to the node that 
        is closest to each point 
        
        """         
        nodeIds = self.grid.get_node_output_locations(xy,epsgIN)
        if(elementIds != []):
            self.run_nc.add_node_output_locations(nodeIds,start,end,step)


    def add_full_grid_output(self,output_filename,output_type, start, step):
        """
        Setup full domain output to the run (i.e. tecplot etc.)
        
        output_type:  output type, integer corresponding to 'nopt' in RiCOM
        start: iteration # where ouput
        
        
        """
        self.ricom.nopt = output_type        
        self.ricom.noptstart = start
        self.ricom.nskip = step
        self.ricom.outputFileFull = output_filename
        
    def add_key_points_output_locations(self, key_points, epsgIN):
        """
        Given a list of key points [ x  ,  y, "name", "description"] where output is required
        
        Locate the closest element and node in the grid and request them for output

        
        """ 
        
        if (key_points != []) and (len(key_points[0]) == 4):    #check if the key point array is valid

            if self.keyPointsAdded == False:                    #check if key points have already been added
                
                elementIds, distanceElements = self.grid.get_closest_elements(key_points,epsgIN)
                nodeIds, distanceNodes = self.grid.get_closest_nodes(key_points,epsgIN)

                
                i = 0 
                while i < len(key_points):
                    if (elementIds[i] == 0) or (nodeIds[i] == 0):
                        elementIds.pop(i)
                        nodeIds.pop(i)
                        key_points.pop(i)
                    i += 1
                
                if len(key_points) <= 0:
                    print "WARNING: Key points NOT added. Pts not valid in domain"
                    return -1
                
                
                print nodeIds
                print elementIds
                
                self.run_nc.add_key_points_output_locations(key_points,nodeIds,distanceNodes,elementIds,distanceElements)
                self.keyPointsAdded = True
            else:
                print "Warning: Key Points already added to Run."
                return -1   
        
        else:
            print "Error: Invalid key point array."


    def add_subdomain_output_by_area_id(self,filename,area_id, start,stop,step):
        """
        Given an area id defined for the grid create an output subdomain that is equal to the bounds of the area geometry
        
        """
        
        #check to see if the start and stop are in the time of the run
        if (stop > self.ricom.nitn): stop = self.ricom.nitn        
        if (start > self.ricom.nitn): return -1             #start is after run finishes - no output required    
        
        bounds = self.grid.get_area_bounds_by_id(area_id)
        points = []
        if(bounds == -1):
            print "ERROR: Invalid area!"
        else:
            i = bounds.find('(')
            j = bounds.find(')')
            bounds = bounds[i+1:j]
            bounds = bounds.replace(',',' ').split()
            
            points.append([float(bounds[0]),float(bounds[1])])
            points.append([float(bounds[2]),float(bounds[3])])

            bounds = convert_points(points,self.grid.epsg,4326)         #convert bounds to Lat Long
            bounds = LL2LocalRicom(bounds, self.ricom.lat0, self.ricom.long0, self.ricom.latoff, self.ricom.longoff)                      #convert bounds to Local Ricom Coordinated
            self.run_nc.add_subdomain_output(filename,bounds[0][0],bounds[0][1], bounds[1][0], bounds[1][1],start,stop,step,area_id)
    
    def add_subdomain_output(self,filename,ll_x,ll_y, ur_x, ur_y, epsgIN,start,stop,step):
        """
        Add a subdomain to be outputted.
        
        ll_x - lower left corner X
        ll_y - lower left corner Y
        ur_x - upper right corner X
        ur_y - upper right corner Y        
        
        """
        bounds = []
        bounds.append([ll_x,ll_y])
        bounds.append([ur_x,ur_y])
        bounds = convert_points(bounds,epsgIN,4326)         #convert bounds to Lat Long

        bounds = LL2LocalRicom(bounds, self.ricom.lat0, self.ricom.long0, self.ricom.latoff, self.ricom.longoff)                      #convert bounds to Local Ricom Coordinated
        self.run_nc.add_subdomain_output(filename,bounds[0][0],bounds[0][1], bounds[1][0], bounds[1][1],start,stop,step)
        
        #self.grid.get_grid_bounds()
        
        

    def write_rcm(self,rcm_filename):
        """
        Write the .RCM file
        
        """
        
        if(self.buildingsAdded != True):
            self.run_nc.add_building_output_locations({}, 0, 0,0)       #Set building locations # 0 in NETCDF file
    
        if(self.keyPointsAdded != True):
            self.run_nc.add_key_points_output_locations([], 0, 0, 0, 0)                                #Set key points to 0 in netcdf file
        
        self.ricom.write_rcm(rcm_filename)
        #self.run_nc.close()


    def _write_output_(self):
        """
        """
        if(self.grid.is_run_defined(self.run_id)):
            print "Adding node output to the database..."
            
            nodeIds, eta, uv, time = self.run_nc.get_output_nodes()
            
            i = 0
            if nodeIds != []:
                if (eta != [] and uv != []):
                    for node_id in nodeIds: 
                        self.grid.add_output_at_node(self.run_id, node_id, eta[i], uv[i], time)
                elif (eta != [] and uv == []):
                    for node_id in nodeIds: 
                        self.grid.add_output_at_node(self.run_id, node_id, eta[i], uv, time)
                elif (eta == [] and uv != []):
                    for node_id in nodeIds: 
                        self.grid.add_output_at_node(self.run_id, node_id, eta, uv, time)


            i = 0
            print "Adding element output to the database..."


            elementIds, eta, uv, time = self.run_nc.get_output_elements()
            if elementIds != []:
                if (eta != [] and uv != []):
                    for element_id in elementIds: 
                        self.grid.add_output_at_element(self.run_id, element_id, eta[i], uv[i], time)
                elif (eta != [] and uv == []):
                    for element_id in elementIds: 
                        self.grid.add_output_at_element(self.run_id, element_id, eta[i], uv, time)
                elif (eta == [] and uv != []):
                    for element_id in elementIds: 
                        self.grid.add_output_at_element(self.run_id, element_id, eta, uv, time)

            i = 0
            print "Adding side output to the database..."

            sideIds, eta, uv, time = self.run_nc.get_output_sides()
            
            if sideIds != []:
                if (eta != [] and uv != []):
                    for side_id in sideIds: 
                        self.grid.add_output_at_side(self.run_id, side_id, eta[i], uv[i], time)
                elif (eta != [] and uv == []):
                    for side_id in sideIds: 
                        self.grid.add_output_at_side(self.run_id, side_id, eta[i], uv, time)
                elif (eta == [] and uv != []):
                    for side_id in sideIds: 
                        self.grid.add_output_at_side(self.run_id, side_id, eta, uv[i], time)
            
        else:
            print "RUN not defined in database"
     
        
    def init_friction(self,ntype = 1, friction_parameters = []):
    	'''
    	Initialise friction for the model run
    	
    	'''
    	
    	#ERROR checking
    	if(len(friction_parameters) != ntype):
    		print "Number of friction types does not equal ntype.  Exiting..."
    		sys.exit()
    	
        i = 1
        for frc in friction_parameters:
        	if frc[0] != i:
        		"Friction parameters list is not sequential. Please revise. Exiting...."
        		sys.exit()
        	if len(frc) != 8:
        		"Incorrect number of parameter in friction type row. Please Revise. Exiting..."
        		sys.exit()
        	i+=1
        
        #Error checks PASSED - Continue
        self.friction_parameters = friction_parameters
        #friction = FrictionType()
        #friction.fricTypeList = friction_parameters
        #friction.ntype = self.ntype
        
        self.ricom.friction = friction_parameters
        self.ricom.ntype = ntype
        self.ntype = ntype
            
    def set_friction_in_area(self,area_id, element_code):
        """
        Set element code of all elements inside a given area to passed element_code
        
        area_id - the id (in the postgis table) of the area
        element_code - target element code of areas
            
        NOTE: if area_id = 0 all elements in the domain are set
        
        """
        
        #ERROR Checking - check if the element code is defined yet
        if element_code > self.ntype:
            print "The element code is not defined in the friction parameters list. Please Revise. Exiting..."
            sys.exit()
        
        elements_in_area = self.grid.get_elements_in_area(area_id) 
        for el in elements_in_area:
            self.elements[el[0]-1][1] = element_code

		
    def init_form_drag(self, ifdrag = -1, number_of_drag_types = 1, drag_parameters = []):
        '''
        Initialise form drag for the model run
        
        '''
        #ERROR checking
        if(len(drag_parameters) != number_of_drag_types):
            print "Number of drag types in drag_parameter list does not number_of_drag_types.  Exiting..."
            sys.exit()
    
        if(ifdrag < -2 or ifdrag > 0): 
            print "Undefined ifdrag.  Choose: 0,-1,-2.  Exiting..."
            sys.exit()
        
        for fd in drag_parameters:
            if fd[0] > self.ntype or fd[0] < 1:
                "Element code not defined in friction parameter. Please Revise. Exiting...."
                sys.exit()
        
        #Error checks PASSED - Continue
        self.drag_parameters = drag_parameters
        self.number_of_drag_types = number_of_drag_types
        self.ricom.ifdrag = ifdrag 
        self.ifdrag = ifdrag 
        
    def set_form_drag_in_area(self, area_id, element_code):
        """
        Set the drag for elements in a given area
        
        area_id - the id (in the postgis table) of the target area
        element_code - target element code of elements in the areas
        
        NOTE: 	The element_code must correspond to an element code number in the drag_parameters list
        		and the friction parameters list
        		
        		Drag definitions can be assigned multiple time to an element, however, only
        		the last defined code will be used for computations in RiCOM
        		
        """
        #ERROR Checking
        dragDefined = False
        for fd in self.drag_parameters:
            if fd[0] == element_code:
                dragDefined = True
        
        #if the element code is defined in the drag_parameters list
        if dragDefined:
            elements_in_area = self.grid.get_elements_in_area(area_id)              #get all elements in the area
            for el in elements_in_area:                                             #save elements to the form_drag_elements array
                self.form_drag_elements.append([el[0],element_code])         
        else:
            print "Given element_code is is not defined in the drag_parameters list. Exiting..."
            sys.exit()
        


    def set_form_drag_inside_buildings(self, area_id, element_code):
        """
        Set the drag for elements inside the building areas
        
        area_id - the id (in the postgis table) of the target area
        element_code - target element code of elements in the areas
        
        NOTE:     The element_code must correspond to an element code number in the drag_parameters list
                and the friction parameters list
                
                Drag definitions can be assigned multiple time to an element, however, only
                the last defined code will be used for computations in RiCOM
                
        """
        
        elements_in_area = []
        #ERROR Checking
        dragDefined = False
        for fd in self.drag_parameters:
            if fd[0] == element_code:
                dragDefined = True
        
        #if the element code is defined in the drag_parameters list
        if dragDefined:
            elements_in_area = self.grid.get_elements_inside_buildings(area_id)              #get all elements in the area
            for el in elements_in_area:                                             #save elements to the form_drag_elements array
                self.form_drag_elements.append([el[0],element_code])         
        else:
            print "Given element_code is is not defined in the drag_parameters list. Exiting..."
            sys.exit()

        return len(elements_in_area)
    

    def get_nodes_at_building_edges(self, area_id, node_code):
        """

        """
        
        
        #outfile = open('NodesAtBuildings_code2.dat',"w")

        nodes_at_buildings = []
        nodes_all = []
        
        nodes_all, nodes_at_buildings = self.grid.get_nodes_at_building_edges(area_id)
        print "Number of Nodes at building edges = %s" % len(nodes_at_buildings)
        
        '''

        for node in nodes_at_buildings:
            nodes_all[node[0]-1][1] = 2
        
        for n in nodes_all:
            if n[1] == 10 or n[1] == 100:
                n[1] = 1
            outfile.write("%s %s\n" % (str(n[0]).ljust(12),n[1]))
            
        
        outfile.close()
        '''
    def adjust_ieadj(self):
        '''
        
        
        
        
        '''
        
        self.ieadj = self.grid.adjust_ieadj()
        

class RunNC:
    
    """
    class GridNC
    
    This class reads and writes existing RiCOM netcdf4 grid file
    INPUT:  filename - netcdf4 grid file from PreMOD
            grid - a NCGrid class object containing the grid data
            flag - new output    
    
   
    *****************************************************************************
    NOTE:  Netcdf does not provide functionality for deleting variables or groups!
    *****************************************************************************
   
   
   
    """

    
    def __init__(self, run_filename,grid_filename,run_id,
                 name, description, grid_dbname,grid_db_epsg,buildings_dbname, user = "tbone"):
        '''


        Initializing function for RunNC object


        '''

        self.grid_filename = grid_filename
        self.run_filename = run_filename
        self.run_id = run_id
        self.grid_dbname = grid_dbname
        self.grid_db_epsg = grid_db_epsg
        self.user = user   
        self.buildings_dbname = buildings_dbname   
        self.name = name   
        self.description = description   
                
        self.buildingsAdded = False
        self.keyPointsAdded = False
        
        self.number_of_subdomains = 0
        
        #Open the Netcdf GRID file output from PreMOD
        try: self.dataset = Dataset(self.run_filename,'w',format='NETCDF4')
        except Exception, e:
            print "ERROR: %s" % e            

        if(self.run_filename == ""):
            print "ERROR: No Run Filename provided"
            #return -1
        
       
        #Create all the BUILDINGS groups
        try: self.buildings = self.dataset.createGroup('buildings')
        except Exception, e:
            print "WARNING: %s" % e
            self.buildings = self.dataset.groups['buildings']
        
        
        try: self.building_nodes = self.buildings.createGroup('nodes')
        except Exception, e:
            print "WARNING: %s" % e           
            self.building_nodes = self.buildings.groups['nodes']
            
        try: self.building_elements = self.buildings.createGroup('elements')
        except Exception, e:
            print "WARNING: %s" % e            
            self.building_elements = self.buildings.groups['elements']   
      
        try: self.building_sides = self.buildings.createGroup('sides')
        except Exception, e:
            print "WARNING: %s" % e            
            self.building_sides = self.buildings.groups['sides']
        
        
        
        self.building_nodes.eta_output_added = 0
        self.building_nodes.uv_output_added = 0
        self.building_elements.eta_output_added = 0
        self.building_elements.uv_output_added = 0        
        self.building_sides.eta_output_added = 0
        self.building_sides.uv_output_added = 0

        
        
        #Create the KEY POINT output groups.  These are the nodes and elements
        # that are not part of a building where output is required but  
        
        try: self.key_points = self.dataset.createGroup('key_points')
        except Exception, e:
            print "WARNING: %s" % e            
            self.key_points = self.dataset.groups['key_points']    

        try: self.key_nodes = self.key_points.createGroup('nodes')
        except Exception, e:
            print "WARNING: %s" % e           
            self.key_nodes = self.key_points.groups['nodes']
       
        try: self.key_elements = self.key_points.createGroup('elements')
        except Exception, e:
            print "WARNING: %s" % e            
            self.key_elements = self.key_points.groups['elements']         
        


        #Create the KEY POINT output groups.  These are the nodes and elements
        # that are not part of a building where output is required but  
        
        try: self.subdomains = self.dataset.createGroup('subdomains')
        except Exception, e:
            print "WARNING: %s" % e            
            self.subdomains = self.dataset.groups['subdomains']    
        
        
        #ATTRIBUTE
        self.subdomains.number_of_subdomains = 0
        self.subdomainGroups = []

        self.dataset.grid_filename = self.grid_filename
        self.dataset.run_filename = self.run_filename
        self.dataset.run_id = self.run_id
        self.dataset.name = self.name
        self.dataset.description = self.description        
        self.dataset.grid_dbname = self.grid_dbname
        self.dataset.grid_db_epsg = self.grid_db_epsg
        self.dataset.buildings_dbname = self.buildings_dbname
        self.dataset.user = self.user   
		
	
        
    def __del__ (self):
        """
        Class deconstructor
    
        """
        self.dataset.close()

    def close(self):
        self.dataset.close()

    def open(self,type = 'r'):
        """
        open the Netcdf dataset
        
        type: 'r' = read, 'r+' = read/write, 'w' = write (overwrites exisitng)
        
        """
        
        self.dataset = Dataset(self.run_filename,'r',format='NETCDF4')      
    
    def add_building_output_locations(self,dictionary, start,end,step):
        """
        Get id's of the nodes, elements and sides that are inside the buildings polygons
        
        Input:
            area_id - id of the area polygon that you want to get the info inside
        
        
        """ 
        """
        Given a dictionary of building footprints and associated nodes,element and sides, add the values 
        to the netcdf grid file.
        
        The nodes, elements and sides associated with each footprint correspond to the there index in the RiCOM grid file
        
        Dictionary format:
        {id1: {'nodes': [n1, n2,...nn] }, {'elements': [e1,e2,...,en] },{'sides': [s1,s2,...,sn]}, id2: {}, id3 {}, ...., idn {} } 
        
        idn = the id of the building footprint that the node, elements and sides belong to
        
        """
        
        if (dictionary != {}):
            maxNodes = 0
            maxElements = 0
            maxSides = 0
            nodesAll = []
            elementsAll = []
            sidesAll = []
            id = []
            perimeter = []
            type = []
            for row in dictionary.iteritems(): 
                id.append(row[0])              
                n = row[1]['nodes'] 
                e = row[1]['elements']
                s = row[1]['sides']
                perimeter.append(row[1]['perimeter'])
                
                if row[1]['type'] == "BUILDINGS_AS_HOLES":
                    typeNUM = 1
                elif row[1]['type'] == "BUILDINGS_GRIDDED":
                    typeNUM = 2

                elif row[1]['type'] == "BUILDINGS_AS_POINTS":
                    typeNUM = 3
                else:
                    typeNUM = 0
                type.append(typeNUM)
                
                nodesAll.extend(n)
                elementsAll.extend(e)
                sidesAll.extend(s)
                if maxNodes < len(n): maxNodes = len(n)
                if maxElements < len(e): maxElements = len(e)
                if maxSides < len(s): maxSides = len(s)
            
            
            #remove repeated elements, sides and nodes
            nodesAll = list(set(nodesAll))
            elementsAll = list(set(elementsAll))
            sidesAll = list(set(sidesAll))
                        
            print "# elements = %s" % len(elementsAll)
            print "# sides = %s" % len(sidesAll)
            print "# nodes = %s" % len(nodesAll)

            
            #initialise arrays for entry into netcdf file
            nodes = zeros((len(dictionary),maxNodes))
            elements = zeros((len(dictionary),maxElements))
            sides = zeros((len(dictionary),maxSides)) 
            
            i = 0
            for row in dictionary.iteritems():    
                nodes[i,0:(len(row[1]['nodes']))] = row[1]['nodes']
                elements[i,0:(len(row[1]['elements']))] = row[1]['elements']
                sides[i,0:(len(row[1]['sides']))] = row[1]['sides']
                i+=1  
            
            #create dimensions
            try: self.buildings.createDimension('max_number_nodes',maxNodes)
            except Exception, e: print "WARNING: %s" % e
            try: self.buildings.createDimension('max_number_elements',maxElements)
            except Exception, e: print "WARNING: %s" % e
            try: self.buildings.createDimension('max_number_sides',maxSides)
            except Exception, e: print "WARNING: %s" % e
            try: self.buildings.createDimension('number_of_buildings',len(dictionary))
            except Exception, e: print "WARNING: %s" % e        
            try: self.building_nodes.createDimension('number_of_nodes',len(nodesAll))
            except Exception, e: print "WARNING: %s" % e
            try: self.building_elements.createDimension('number_of_elements',len(elementsAll))
            except Exception, e: print "WARNING: %s" % e
            try: self.building_sides.createDimension('number_of_sides',len(sidesAll))
            except Exception, e: print "WARNING: %s" % e
            
            
            #create variables
            try: building_id = self.buildings.createVariable(varname = 'building_id',datatype = 'i', dimensions=('number_of_buildings',)) 
            except Exception, e:
                building_id = self.buildings.variables['building_id']
                print "WARNING: %s" % e
    
            try: building_wkt = self.buildings.createVariable(varname = 'building_wkt',datatype = str, dimensions=('number_of_buildings',)) 
            except Exception, e:
                building_wkt = self.buildings.variables['building_wkt']            
                print "WARNING: %s" % e

            try: building_perimeter = self.buildings.createVariable(varname = 'building_perimeter',datatype = 'd', dimensions=('number_of_buildings',)) 
            except Exception, e:
                building_perimeter = self.buildings.variables['building_perimeter']            
                print "WARNING: %s" % e


            try: building_type = self.buildings.createVariable(varname = 'building_type',datatype = 'i', dimensions=('number_of_buildings',)) 
            except Exception, e:
                building_type = self.buildings.variables['building_type']            
                print "WARNING: %s" % e

            try: building_nodes = self.buildings.createVariable(varname = 'building_nodes',datatype = 'i', dimensions=('number_of_buildings','max_number_nodes',)) 
            except Exception, e:
                building_nodes = self.buildings.variables['building_nodes']            
                print "WARNING: %s" % e
    
            try: building_elements = self.buildings.createVariable(varname = 'building_elements',datatype = 'i', dimensions=('number_of_buildings','max_number_elements',)) 
            except Exception, e:
                building_elements = self.buildings.variables['building_elements']
                print "WARNING: %s" % e
            
            try: building_sides = self.buildings.createVariable(varname = 'building_sides',datatype = 'i', dimensions=('number_of_buildings','max_number_sides',)) 
            except Exception, e:
                building_sides = self.buildings.variables['building_sides']
                print "WARNING: %s" % e
            
            building_nodes[:] = nodes
            building_elements[:] = elements
            building_sides[:] = sides
            building_id[:] = array(id) 
            building_perimeter[:] = array(perimeter)
            building_type[:] = array(type)
            #Set the attributes
            self.building_nodes.start = start
            self.building_nodes.finish = end
            self.building_nodes.step = step
            self.building_elements.start = start
            self.building_elements.finish = end
            self.building_elements.step = step
            self.building_sides.start = start
            self.building_sides.finish = end
            self.building_sides.step = step
            
            #assign the data
            output_ids = {'nodes': [], 'elements': [], 'sides': []}
            try: output_ids['nodes'] = self.building_nodes.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_nodes',))
            except Exception, e:
                output_ids['nodes'] = self.building_nodes.variables['id']
                print "WARNING: %s" % e
            try: output_ids['elements'] = self.building_elements.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_elements',))
            except Exception, e:
                output_ids['elements'] = self.building_elements.variables['id']
                print "WARNING: %s" % e
            try: output_ids['sides'] = self.building_sides.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_sides',))
            except Exception, e:
                output_ids['sides'] = self.building_sides.variables['id']
                print "WARNING: %s" % e
    
    
            output_ids['nodes'][:] = array(nodesAll)
            output_ids['elements'][:] = array(elementsAll)
            output_ids['sides'][:] =  array(sidesAll)
            
            
            self.buildingsAdded = True
        else:
            #create dimensions
            try: self.buildings.createDimension('number_of_buildings',0)
            except Exception, e: print "WARNING: %s" % e        
            try: self.building_nodes.createDimension('number_of_nodes',0)
            except Exception, e: print "WARNING: %s" % e
            try: self.building_elements.createDimension('number_of_elements',0)
            except Exception, e: print "WARNING: %s" % e
            try: self.building_sides.createDimension('number_of_sides',0)
            except Exception, e: print "WARNING: %s" % e  
            self.buildingsAdded = True
				  
    
    def add_key_points_output_locations(self, key_points,nodeIds,distanceNodes, elementIds, distanceElements):
        """
        Given a list of key points [ x  ,  y, "description"] where output is required
        
        Locate the closest element and node in the grid and request them for output

        
        """ 
       
       
        if (key_points != []):
            try: self.key_points.createDimension('number_of_key_points',len(key_points))
            except Exception, e: print "WARNING: %s" % e
    
            try: self.key_points.createDimension('dimensions',2)
            except Exception, e: print "WARNING: %s" % e
                  
            try: nodesVar = self.key_nodes.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_key_points',))
            except Exception, e:
                nodesVar = self.key_nodes.variables['id']
                
            try: elementsVar = self.key_elements.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_key_points',))
            except Exception, e:
                elementsVar = self.key_elements.variables['id']             
    
    
    
            try: distanceNodesVar = self.key_nodes.createVariable(varname = 'distance',datatype = 'f', dimensions=('number_of_key_points',))
            except Exception, e:
                distanceNodesVar = self.key_nodes.variables['distance']
            
            
            try: distanceElementsVar = self.key_elements.createVariable(varname = 'distance',datatype = 'f', dimensions=('number_of_key_points',))
            except Exception, e:
                distanceElementsVar = self.key_elements.variables['distance'] 
            
            nodesVar[:] = array(nodeIds)
            elementsVar[:] = array(elementIds)    
            distanceNodesVar[:] = array(distanceNodes)
            distanceElementsVar[:] = array(distanceElements)        
            
            xy = []
            description = []
            name = []

            for p in key_points:
                xy.append([p[0],p[1]])
                name.append(p[2])
                description.append(p[3])
            
            try: xyVar = self.key_points.createVariable(varname = 'xy',datatype = 'f', dimensions=('number_of_key_points','dimensions'))
            except Exception, e:
                xyVar = self.key_points.variables['xy']
                
            try: desciptionVar = self.key_points.createVariable(varname = 'description',datatype = str, dimensions=('number_of_key_points',))
            except Exception, e:
                desciptionVar = self.key_points.variables['description'] 

            try: nameVar = self.key_points.createVariable(varname = 'name',datatype = str, dimensions=('number_of_key_points',))
            except Exception, e:
                nameVar = self.key_points.variables['name'] 
                
            xyVar[:] = array(xy)
    
            #add the description data
            data = empty(len(description),'O')
            i = 0
            for d in description:
                data[i] = d
                i+=1        
            desciptionVar[:] = data 
           
            data = empty(len(name),'O')
            i = 0
            for d in name:
                data[i] = d
                i+=1        
            nameVar[:] = data
            self.keyPointsAdded = True  
        
        else:    

            try: self.key_points.createDimension('number_of_key_points',0)
            except Exception, e: print "WARNING: %s" % e
        
            self.keyPointsAdded = True  


    def add_subdomain_output(self,filename,ll_x,ll_y, ur_x, ur_y,start,stop,step,area_id = 0):
        
        """
        Add output at a subdomain
        
        """        
        
        self.number_of_subdomains += 1
        self.subdomains.number_of_subdomains = self.number_of_subdomains            #set the 'number_of_subdomains' attribute 
        name = 'subdomain' + str(self.number_of_subdomains) 
        self.subdomainGroups.append(self.subdomains.createGroup(name) )             #great a new subdomain Group
        
        self.subdomainGroups[self.number_of_subdomains-1].filename = filename       #set the bounds attributes for the subdomain

        self.subdomainGroups[self.number_of_subdomains-1].ll_x = ll_x               #set the bounds attributes for the subdomain
        self.subdomainGroups[self.number_of_subdomains-1].ll_y = ll_y
        self.subdomainGroups[self.number_of_subdomains-1].ur_x = ur_x
        self.subdomainGroups[self.number_of_subdomains-1].ur_y = ur_y
        self.subdomainGroups[self.number_of_subdomains-1].start = start
        self.subdomainGroups[self.number_of_subdomains-1].stop = stop
        self.subdomainGroups[self.number_of_subdomains-1].step = step
        self.subdomainGroups[self.number_of_subdomains-1].area_id = area_id


    def get_output_sides(self):
        """
        Reads the netcdf grid file and returns the side outputs 
    
        """

        self.buildings = self.dataset.groups['buildings']
        self.building_sides = self.buildings.groups['sides']

        eta_output_added = getattr(self.building_sides,'eta_output_added')
        uv_output_added = getattr(self.building_sides,'uv_output_added')
        eta = []
        uv = []
        sideIds = []
        time = []
        if(eta_output_added or uv_output_added ):
            time = self.building_sides.variables['time'][:].tolist()
            sideIds = self.building_sides.variables['id'][:].tolist()
            if eta_output_added: eta = self.building_sides.variables['eta'][:].tolist()
            if uv_output_added: uv = self.building_sides.variables['uv'][:].tolist()
        
        return sideIds,eta, uv, time

    def get_output_elements(self):
        """
        Reads the netcdf grid file and returns the element outputs 
    
        """

        self.buildings = self.dataset.groups['buildings']
        self.building_elements = self.buildings.groups['elements']

        eta_output_added = getattr(self.building_elements,'eta_output_added')
        uv_output_added = getattr(self.building_elements,'uv_output_added')
        eta = []
        uv = []
        elementIds = []
        time = []
                
        if(eta_output_added or uv_output_added ):
            time = self.building_elements.variables['time'][:].tolist()
            elementIds = self.building_elements.variables['id'][:].tolist()
            if eta_output_added: eta = self.building_elements.variables['eta'][:].tolist()
            if uv_output_added: uv = self.building_elements.variables['uv'][:].tolist()
        
        return elementIds,eta, uv, time
        

    def get_output_nodes(self):
        """
        Reads the netcdf grid file and returns the nodal outputs at the requested building id
    
        """
  

        self.buildings = self.dataset.groups['buildings']
        self.building_nodes = self.buildings.groups['nodes']

        eta_output_added = getattr(self.building_nodes,'eta_output_added')
        uv_output_added = getattr(self.building_nodes,'uv_output_added')

        eta = []
        uv = []
        nodeIds = []
        time = []
        
        if(eta_output_added or uv_output_added ):
            time = self.building_nodes.variables['time'][:].tolist()
            nodeIds = self.building_nodes.variables['id'][:].tolist()
            if eta_output_added: eta = self.building_nodes.variables['eta'][:].tolist()
            if uv_output_added: uv = self.building_nodes.variables['uv'][:].tolist()

        
        return nodeIds,eta, uv, time


    def set_element_codes(self,elements,number_of_element_types = 1):
        """
        Set the element codes within a given area. 
        
        USE: To set the element codes for the grid at runtime.
        
        """         
        
        #create dimensions
        try: self.dataset.createDimension('number_of_elements_in_domain',len(elements))
        except Exception, e: print "WARNING: %s" % e
        
        #create dimensions
        try: self.dataset.createDimension('number_of_element_types',number_of_element_types)
        except Exception, e: print "WARNING: %s" % e
        
        
        try: element_codes = self.dataset.createVariable(varname = 'element_codes',datatype = 'i', dimensions=('number_of_elements_in_domain',)) 
        except Exception, e:
        	element_codes = self.dataset.variables['element_codes']
        	print "WARNING: %s" % e
        
        codes = []
        for el in elements:
        	codes.append(el[1])
        	
        element_codes[:] = array(codes) 


    def set_ieadj(self,ieadj_in):
 
 

        #create dimensions
        try: self.dataset.createDimension('ieadj_dimensions',5)
        except Exception, e: print "WARNING: %s" % e

        try: ieadj = self.dataset.createVariable(varname = 'ieadj',datatype = 'i', dimensions=('number_of_elements_in_domain','ieadj_dimensions')) 
        except Exception, e:
            ieadj = self.dataset.variables['ieadj']
            print "WARNING: %s" % e
        

        ieadj[:] = array(ieadj_in) 



        '''
        
        dict_id = {}        #dictionary maps the the node id to the location in the id array
        n = 0
        for id in node_id:
            dict_id[id] = n 
            n += 1
                    
        dict_nodes = {}     #dictionary maps the footprint id to the nodal output from the netcdf file
        dict_nodes['time'] = time.tolist()
        n=0
        for id in building_id:
            dict_nodes[id] = {'node_ids': [], 'eta': {}, 'uv': {}}
            nodesList = []
            for node in building_nodes[n]:
                if node != 0:
                    dict_nodes[id]['eta'][node] = eta[dict_id[node]]
                    dict_nodes[id]['uv'][node] = uv[dict_id[node]] 
                    nodesList.append(node)
                else: break            
            dict_nodes[id]['node_ids'] = nodesList
            n+=1

        
        id,xyz,neigh,code = self.get_grid_nodes()
        neigh = []
        code = []
        fp_xyz = []
        fp_uv = []
        fp_eta = []

        #xyz = array(xyz)

        #Get footprint nodes
        fp_ids = dict_nodes[1]['node_ids']
        for id in fp_ids:
            fp_xyz.append(xyz[id-1])
            fp_uv.append(dict_nodes[1]['uv'][id])
            
       
        fp_xyz = array(fp_xyz) 
        fp_uv = array(fp_uv)
        
        x = fp_xyz[:,0]
        y = fp_xyz[:,1]
        uv = fp_uv
        
        #q = plt.quiver(x,y,u,v,angles='xy',scale=1000,color='r')
        #plt.show()
            
        
        return dict_nodes,x,y,uv
        
#        return dict_nodes
        '''


        
    '''
    def get_building_nodes(self, building_id = 0):
        """
        Reads the netcdf grid file and returns the nodal outputs at the requested building id
    
        """
        
        import matplotlib.pyplot as plt
        
        eta = self.building_nodes.variables['eta'][:]
        time = self.building_nodes.variables['time'][:]
        uv = self.building_nodes.variables['uv'][:]
        node_id = self.building_nodes.variables['id'][:]
        building_id = self.buildings.variables['building_id'][:]
        building_nodes = self.buildings.variables['building_nodes'][:]
        
        dict_id = {}        #dictionary maps the the node id to the location in the id array
        n = 0
        for id in node_id:
            dict_id[id] = n 
            n += 1
                    
        dict_nodes = {}     #dictionary maps the footprint id to the nodal output from the netcdf file
        dict_nodes['time'] = time.tolist()
        n=0
        for id in building_id:
            dict_nodes[id] = {'node_ids': [], 'eta': {}, 'uv': {}}
            nodesList = []
            for node in building_nodes[n]:
                if node != 0:
                    dict_nodes[id]['eta'][node] = eta[dict_id[node]]
                    dict_nodes[id]['uv'][node] = uv[dict_id[node]] 
                    nodesList.append(node)
                else: break            
            dict_nodes[id]['node_ids'] = nodesList
            n+=1

        
        id,xyz,neigh,code = self.get_grid_nodes()
        neigh = []
        code = []
        fp_xyz = []
        fp_uv = []
        fp_eta = []

        #xyz = array(xyz)

        #Get footprint nodes
        fp_ids = dict_nodes[1]['node_ids']
        for id in fp_ids:
            fp_xyz.append(xyz[id-1])
            fp_uv.append(dict_nodes[1]['uv'][id])
            
       
        fp_xyz = array(fp_xyz) 
        fp_uv = array(fp_uv)
        
        x = fp_xyz[:,0]
        y = fp_xyz[:,1]
        uv = fp_uv
        
        #q = plt.quiver(x,y,u,v,angles='xy',scale=1000,color='r')
        #plt.show()
            
        
        return dict_nodes,x,y,uv
        
#        return dict_nodes
    
    def get_building_elements(self, building_id = 0):
    
        """
        Reads the netcdf grid file and returns the element outputs at the requested building id
    
        """
        eta = []
        uv = []
        
        return eta, uv
    
    def get_building_sides(self, building_id = 0):
    
        """
        Reads the netcdf grid file and returns the element outputs at the requested building id
    
        """
        eta = []
        uv = []
        
        return eta, uv
    
    
    
    
    def get_output_footprints_max(self):
        
        import math
        footprints_max = {}
        dict_ids = {'nodes':{}, 'elements': {}, 'sides': {}, 'footprints': {}, 'time': {}}        #dictionary maps the the node id to the location in the id array
        building_ids = self.buildings.variables['building_id'][:]
        node_id = self.building_nodes.variables['id'][:]
        time = self.building_nodes.variables['time'][:]
        building_nodes = self.buildings.variables['building_nodes'][:]
        xyz = self.domain.variables['xyz'][:]
        eta_all = self.building_nodes.variables['eta'][:]
        uv_all = self.building_nodes.variables['uv'][:]
        
        n = 0
        for id in node_id:
            dict_ids['nodes'][id] = n 
            n += 1
            
        n = 0
        for id in building_ids:
            dict_ids['footprints'][id] = n 
            n += 1
        
       
        i_footprint = 0
        for id in building_ids:            #iterate throught the footprints
            print "Footprint Id = %s\n" % str(id)
            nodes = building_nodes[i_footprint]
            i_time = 0
            maxData = []
            for step in time:
                fd_max = 0
                speed_max = 0
                for node in nodes:          #For current footprint: get the max node data for the time step
                    if node != 0:
                        #eta = self.building_nodes.variables['eta'][dict_ids['nodes'][node]][i_time]
                        eta = eta_all[dict_ids['nodes'][node]][i_time]
                        uv = uv_all[dict_ids['nodes'][node]][i_time]
                        #uv = self.building_nodes.variables['uv'][dict_ids['nodes'][node]][i_time]
                        speed  = math.sqrt(uv[0]*uv[0] + uv[1]*uv[1])
                        if eta - xyz[node-1][2] < 0.01: fd = 0
                        else: fd = eta - xyz[node-1][2] 
                        
                        if fd > fd_max: fd_max = fd
                        if speed > speed_max: speed_max = speed
    
                    else: break
        
                i_time += 1
                maxData.append([step, fd_max, speed_max])

            i_footprint += 1
            maxData = array(maxData)
            fd_absMax = maxData[:,1].max()
            speed_absMax = maxData[:,2].max()
            footprints_max[id] = {'time_series': maxData, 'fd_max': fd_absMax, 'speed_max': speed_absMax }
            
#            print "id = %s, fd_max = %s, speed_max = %s\n" % (str(id), str(footprints_max[id]['fd_max']), str(footprints_max[id]['speed_max']))
        return footprints_max
        
                    
            
    def get_output_footprints(self, building_id = 2):     
           
  
        footprint = {}
        dict_ids = {'nodes':{}, 'elements': {}, 'sides': {}, 'footprints': {}}        #dictionary maps the the node id to the location in the id array
        building_ids = self.buildings.variables['building_id'][:]
        node_id = self.building_nodes.variables['id'][:]
        element_id = self.building_elements.variables['id'][:]
        side_id = self.building_sides.variables['id'][:]

        n = 0
        for id in node_id:
            dict_ids['nodes'][id] = n 
            n += 1

        n = 0
        for id in element_id:
            dict_ids['elements'][id] = n 
            n += 1
            
        n = 0
        for id in side_id:
            dict_ids['sides'][id] = n 
            n += 1
                    
        n = 0
        for id in building_ids:
            dict_ids['footprints'][id] = n 
            n += 1
        
                
        nodes = self.buildings.variables['building_nodes'][dict_ids['footprints'][building_id]]
        elements = self.buildings.variables['building_elements'][dict_ids['footprints'][building_id]]
        sides = self.buildings.variables['building_sides'][dict_ids['footprints'][building_id]]
            
        for n in nodes:
            if n != 0:
                xyz = self.domain.variables['xyz'][n-1]
                eta = self.building_nodes.variables['eta'][dict_ids['nodes'][n]]
                i = 0
                while (i < len(eta)):    
                    if eta[i] - xyz[2] < 0.01: eta[i] = 0
                    else: eta[i] = eta[i] - xyz[2] 
                    i+=1
                    
                uv = self.building_nodes.variables['uv'][dict_ids['nodes'][n]]
            else: break            

        for e in elements:
            if e != 0:
                element_nodes = self.domain.variables['nen'][e-1]
                eta = self.building_elements.variables['eta'][dict_ids['elements'][e]]
                uv = self.building_elements.variables['uv'][dict_ids['elements'][e]]
            else: break 

        for s in sides:
            if s != 0:
                side_nodes = self.domain.variables['sides'][s-1]        #get the node numbers of the side
                xyz0 = self.domain.variables['xyz'][side_nodes[0]-1]    #get the node coordinates for node0 of the side
                xyz1 = self.domain.variables['xyz'][side_nodes[1]-1]    #get the node coordinates for node1 of the side
                eta = self.building_elements.variables['eta'][dict_ids['sides'][s]]
                uv = self.building_elements.variables['uv'][dict_ids['sides'][s]]
            else: break 
  

    '''


class FrictionType:
       
    def __init__(self):

        self.fricTypeList = []
        self.number_of_friction_types = 0
        #specifies if the friction types will be written in an external file or .rcm file
        #default: location = 0, write the friction types to an external file
        self.location = 1
        self.ntype = 0             

    #------------------------------------------------------------------------------
    #add a friction type to the the friction type list
    #------------------------------------------------------------------------------

    def add_friction_type(self,type_id,a0 = 0.01,a1 = 0.1,r0 = 0.0,r1 = 0.0025,cd = 0.0,ah = 0.0,qinf = 0.0 ):
        
        i = 0
        while i < self.number_of_friction_types:
            if self.fricTypeList[i][0] == type_id:
                print "Warning: Friction Type = %s already exists.  Overwriting!" % type_id
                self.fricTypeList.pop(i)
                self.number_of_friction_types = self.number_of_friction_types - 1
            i += 1
                
        
        self.fricTypeList.append([ int(type_id), float(a0), float(a1), float(r0), float(r1),
                               float(cd), float(ah), float(qinf)])
        self.number_of_friction_types += 1
        
        if self.location == 0: self.ntype = -1 * self.number_of_friction_types
        elif self.location == 1: self.ntype = self.number_of_friction_types
   
    #------------------------------------------------------------------------------
    #add friction types from an existing friction file
    #type_id a0 a1 r0 r1 cd ah qinf 
    #------------------------------------------------------------------------------
    
    def add_friction_from_file(self, fricFilename = "friction.par"):
        
        self.number_of_friction_types = 0
        self.fricTypeList = []
        file_in = open(fricFilename,"r").readlines()        
        for f in file_in:
            self.number_of_friction_types += 1
            f = f.split()
            self.fricTypeList.append([ int(f[0]), float(f[1]), float(f[2]), float(f[3]), float(f[4]),
                                    float(f[5]), float(f[6]), float(f[7])])

        for f in self.fricTypeList:
                match = 0
                for row in self.fricTypeList:
                    if row[0] == f[0]:
                        match += 1
                        
                if match > 1:
                    print "Warning: Repeated friction types in input file. type_id = %s" % f[0]
        
        if self.location == 0: self.ntype = -1 * self.number_of_friction_types
        elif self.location == 1: self.ntype = self.number_of_friction_types       
    #------------------------------------------------------------------------------
    #Set the location if where the friction will be stored during a ricom run
    #location = 0: friction types in an external file
    #location = 1: friction types stored in .rcm file
    #------------------------------------------------------------------------------

    def set_location(self, location = 0):
        if location < 0 or location > 1:
            print "Invalid friction location.  Location = 0 for external file. location = 1 for to write in .rcm"
        else:
            self.location = location
            if self.location == 0: self.ntype = -1 * self.number_of_friction_types
            elif self.location == 1: self.ntype = self.number_of_friction_types

    #------------------------------------------------------------------------------
    #Write a friction (.par) text file for input into RiCOM
    #------------------------------------------------------------------------------
        
    def write_friction_file(self, fricFilename = "friction.par"):
        
        outfile = open(fricFilename,"w")

        for f in self.fricTypeList:
            outfile.write("%s   %s   %s   %s   %s   %s   %s   %s\n" % (f[0], f[1], f[2], f[3], f[4], f[5], 
                                                                       f[6], f[7]))  
        outfile.close()
                        
    def reset(self):
        self.number_of_friction_types = 0
        self.fricTypeList = []
        
        
##This CLASS contains a single okada fault rupture segment        
class OkadaFault:
    
    def __init__(self):
        
        #private variables
        self.patchList = []
        self.number_of_patches = 0
        
    #------------------------------------------------------------------------------
    #add a patch to the the patch list
    #------------------------------------------------------------------------------

    def add_patch(self, x, y, length, width, depth, delta, strike, slip, rake, lamda = 2e10, mu = 2e10):
        self.patchList.append([ float(x), float(y), float(length), float(width), float(depth),
                               float(delta), float(strike), float(lamda), float(mu), float(slip), float(rake)])
        self.number_of_patches += 1
   
    #------------------------------------------------------------------------------
    #add fault patches from an existing fault file
    #FORMAT:
    #number_of_patches
    #x y length width depth delta strike lamda mu slip rake 
    #------------------------------------------------------------------------------
    
    def add_patches_from_file(self, faultFilename = "fault.param"):
        
        file_in = open(faultFilename,"r").readlines()
        self.number_of_patches = int(file_in[0].split()[0])
        file_in.pop(0)
        
        i = 0
        while i < self.number_of_patches:
            p = file_in[i].split()
            self.patchList.append([ float(p[0]), float(p[1]), float(p[2]), float(p[3]), float(p[4]),
                                    float(p[5]), float(p[6]), float(p[7]), float(p[8]), float(p[9]), float(p[10])])
            i += 1

    #------------------------------------------------------------------------------
    #Write a fault text file for input into RiCOM
    #------------------------------------------------------------------------------
        
    def write_fault_file(self, faultFilename = "fault.param"):
        
        outfile = open(faultFilename,"w")
        
        outfile.write("%s\n" % self.number_of_patches)
        for p in self.patchList:
            outfile.write("%s %s %s %s %s %s %s %s %s %s %s\n" % (p[0], p[1], p[2], p[3], p[4], p[5], 
                                                                p[6], p[7], p[8], p[9], p[10]))
            
        outfile.close()
            
    
                           
    def reset(self):                                     
        self.patchList = []
        self.number_of_patches = 0


class Ricom:
    
    def __init__(	self, 
    				run_dir,
				):
        
        
        self.fault = OkadaFault()               #list of okada fault parameters - see class OkadaSegment
        self.run_filename = ""
        self.run_dir = run_dir
        #---------------------------------------------------------
        self.title = ""
        self.gridfilename = "GridName.nc"
        
        #Line 1      nprt,nsbc,nitn,isolve,nson
        self.nprt = 0
        self.nsbc = 0
        self.nitn = 2000
        self.isolve = 2
        self.nson = 0
        
        #Line 2   -   ivar(i),i=1,4), rvar(1),rvar(2)
        self.maxiter = 700
        self.maxNMiter = 2
        self.ifill = 4
        self.itopt = 2
        self.epsi = 1.0e-10
        self.epsiNM = 0.001
        
        #Line 3
        self.omega0 = -46.0
        self.elev = 0.0
        self.depmin = 0.1
        self.zminq = 0.0

        #Line 4
        self.iomega = 0.0
        self.tload = 1.0
        self.eqfac = 0.0
        self.nload = 0.0
        self.neqtide = 0.0
        self.fn = 1.0
        self.nu0 = 0.0
        
        #Line 5
        self.icoord = 1
        self.lat0 = -46.0
        self.long0 = 155.5
        self.latoff = 0.0
        self.longoff = 0.0
        
        
        #Line 6
        self.delt = 10
        self.tmax = 12.4206012
        self.tscale = 1.0
        self.gamma0 = 0.55
        
        #Line 7
        self.ntype = 1
        self.friction = []
        self.friction.append([1, 0.0,   0.1,   0.0,   0.0025,   0.0,   0.0,   0.0])          #list of friction types, DEFAULT Value for element code = 1
       
        #Line 8
        self.frictionFile = "friction.par"
        
        #Line 9
        self.npvx = 1
        self.izcoord = 0
        self.izgrid = 0
        
        #Line 10
        self.nbx = 0
        self.ibcfile = 0
        self.ncon = 8
        self.nbxfile = 0
        self.irampa = 0
        self.irampq = 0
        
        #Line 11
        self.iwind = 0

        #Line 12
        self.nsed = 0
        
        #Line 13
        self.nsol = 0
        
        #Line 14
        self.ifr = 4
        self.itn = 2
        self.iwn = 5
        self.ivfr = 4
        self.ivsf = 0
        self.ihfr = 0
        self.ifdrag = 0

        #Line 15
        self.irst = 214
        self.irstout = 1000

        #Line 16
        self.nopt = 1
        self.noptstart = 1000000
        self.nskip = 100

        #Line 17
        self.outputFileFull = "OutputFullDomain.dat"
        
        #Line 18
        self.ntsdata = 0
        self.ntsskip = 1
        self.jPprofile = 0
        
        #Line 19
        self.outputElementsETA = []     
         
        #Line 20
        self.ntsUdata = 0
        self.ntsUskip = 1
        self.jUprofile = 1
        
        #Line 21
        self.outputElementsVEL = []         
         
        #Line 22
        self.ntsCdata = 0
        self.ntsCskip = 1
        self.jCprofile = 1
        
        #Line 23
        self.outputElementsC = []       #??? what is this ???
        
        
   
    def read_rcm(self, rcmFilename = "tsunami.rcm"):
        """
        read_rcm()
        ------------------------------------------------------------------------------     
        Read an the passed .RCM file and set the RiCOM class variables
        INPUT: .RCM filename
        
        """
        file_in = open(rcmFilename,"r").readlines()



    def write_rcm(self, rcmFilename = "tsunami.rcm"):
        """
        write_rcm()
        ------------------------------------------------------------------------------     
        Write an a ricom control file (.RCM) based on the Ricom class variables
        INPUT: .RCM filename
        
        """
        
        outfile = open(rcmFilename,"w")
        outfile.write("%s\n%s\n" % (self.title,self.gridfilename))
        
        outfile.write("%s\n" % (self.run_filename))

        outfile.write("%s %s %s %s %s                       nprt,nsbc,nitn,isolve,nson\n" % 
                            (self.nprt,self.nsbc,self.nitn,self.isolve,self.nson))
        outfile.write("%s %s %s %s %s %s              maxiter, maxNMiter, ifill, itopt, epsi, epsiNM\n" % 
                            (self.maxiter,self.maxNMiter,self.ifill,self.itopt,self.epsi,self.epsiNM))
        
        outfile.write("%s %s %s %s                  omega0, elev, depmin, zminq\n" % (self.omega0, self.elev, self.depmin, self.zminq))
        outfile.write("%s %s %s %s %s %s %s      iomega,Tload,eqfac,nload,neqtide,fn,nu0\n" %
                            (self.iomega, self.tload, self.eqfac, self.nload, self.neqtide, self.fn, self.nu0))
        
        outfile.write("%s %s %s %s %s              icoord,lat0,long0,latoff,longoff\n" % 
                            (self.icoord, self.lat0,self.long0, self.latoff, self.longoff))
        
        outfile.write("%s %s %s %s             delt, tmax, tscale, gamma0\n" % (self.delt, self.tmax, self.tscale, self.gamma0))
        
        #--------------------Write FRICTION TYPES-------------------------
        if self.ntype != len(self.friction):
                self.ntype = len(self.friction)
                
        outfile.write("%s                                  ntype\n" % self.ntype)
        if self.ntype > 0:      #write friction types in .rcm file
            for f in self.friction:
                outfile.write("%s   %s   %s   %s   %s   %s   %s   %s\n" % (f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]))  
        
        
        
        #elif self.friction.ntype < 0:      #write friction types in external file
        #    self.friction.write_friction_file(self.run_dir + "/" + self.frictionFile)
        #    outfile.write(self.frictionFile + "\n")
        #------------------------------------------------------------------

    
        outfile.write("%s %s %s                              npvx, izcoord, izgrid\n" % (self.npvx, self.izcoord, self.izgrid))
        outfile.write("%s %s %s %s %s %s                        nbx, ibcfile, ncon, nbxfile, irampa, irampq\n" % 
                        (self.nbx, self.ibcfile, self.ncon, self.nbxfile, self.irampa, self.irampq))

        outfile.write("%s                                  iwind\n%s                                  nsed\n%s                                  nsol\n" %
                      (self.iwind, self.nsed, self.nsol))
        
        outfile.write("%s %s %s %s %s %s %s                      ifr, itn, iwn, ivfr, ivsf, ihfr, ifdrag\n" %
                      (self.ifr, self.itn, self.iwn, self.ivfr, self.ivsf, self.ihfr, self.ifdrag))
        
        #--------------------------OUTPUT OPTIONS-------------------------------

        outfile.write("%s %s                           irst, irstout\n" % (self.irst, self.irstout))   #set up intial conditions and secify restart
        outfile.write("%s %s %s                            nopt noptstart nskip\n" % (self.nopt, self.noptstart, self.nskip))
        if self.nopt > 0:
            outfile.write("%s\n" % self.outputFileFull)
        
        #output ETA at specified elements
        outfile.write("%s %s %s                              ntsdata, ntsskip, jPprofle\n" % (self.ntsdata, self.ntsskip, self.jPprofile))
        if self.ntsdata > 0:    #element numbers
            for ele in self.outputElementsETA:
                outfile.write("%s" % ele)
            outfile.write("\n") 

        #output ETA at specified elements
        outfile.write("%s %s %s                              ntsUdata, ntsUskip, jUPprofle\n" % (self.ntsUdata, self.ntsUskip, self.jUprofile))
        if self.ntsdata > 0:    #element numbers
            for ele in self.outputElementsVEL:
                outfile.write("%s" % ele)
            outfile.write("\n") 
            
            
        #output ETA at specified elements
        outfile.write("%s %s %s                              ntsCdata, ntsCskip, jCPprofle\n" % (self.ntsCdata, self.ntsCskip, self.jCprofile))
        if self.ntsdata > 0:    #element numbers
            for ele in self.outputElementsC:
                outfile.write("%s " % ele)
            outfile.write("\n")      
        
        
        outfile.close()      
        
        #-----------------------------------------------------------------------

    def write_fault_file(self, fault_filename):
        """
        
        
        """
        self.fault.write_fault_file(fault_filename)


    def add_fault_rupture_from_file(self, fault_filename):
        
        """
        ------------------------------------------------------------------------------
        Set the fault rupture intial conditions for the Ricom RUN and write the 
        fault.param input file
        
        irst = 214 makes RiCOM read an external fault rupture file = fault.param
        ------------------------------------------------------------------------------
        """
        self.fault.add_patches_from_file(fault_filename)        
        self.irst = 214

    
    def add_restart(self, restart_filename):
        
        """
        ------------------------------------------------------------------------------
		Use a restart file as the driving input to the model
		
		Function makes a copy of the restart file (restart_filename) and renames the file to restart.bin
		
		NOTE: Ricom will use whatever file is name restart.bin as the restart file
        ------------------------------------------------------------------------------
        """
        import shutil
        shutil.copyfile(restart_filename, self.run_dir + "/restart.bin")
        self.irst = 1
        

    def add_friction_from_file(self, fric_filename):
        """
        Add friction definitions to the rcm file
        """
      
        self.friction.add_friction_from_file(fric_filename)

        
        if(self.friction.location == 0):        #Friction definitions in external file
            self.ntype = self.friction.number_of_friction_types
        else:                                   #Friction definitions defined internally
            self.ntype = -1 * self.friction.number_of_friction_types



class Premod:
    
    
    """
    class Premod
    -----------------------------------------------------------
    Wrapping class to Premod (pre-processing program for RiCOM)
    Requires a working .dylib version of Premod (i.e. dynamic linked library)
    
    *building a .dylib with netcdf4 support
        ifort -dynamiclib -r8 -I/opt/local/include -L/opt/local/lib -lnetcdf -lnetcdff -lcurl -lhdf5 
        -lhdf5_hl /opt/local/lib/libnetcdf.a -o premod.dylib *.f90
    
    """
    
    def __init__(   self,
                    ngh_file="",
                    rtr_file="",
                    netcdf_file="",
                    xscale = 1.0,
                    yscale = 1.0,
                    zscale = 1.0,
                    zoffset = 0.0,
                    zlimit = 10,
                    izup = 1,
                    itest = 1,
                    ireorder = 0,
                    ifront = 0,
                    iord = 1,
                    ifmt = 1,
                    isideopt = 1,
                    ialfa = 0,
                    nbx = 0,
                    iprt = 1,
                    nopt = 2,
                    description = "Tsunami Premod File",
                    title = ""):  
          
        self.xscale = xscale
        self.yscale = yscale
        self.zscale = zscale
        self.zoffset = zoffset
        self.zlimit = zlimit
        self.izup = izup
        self.itest = itest
        self.ireorder = ireorder
        self.ifront = ifront
        self.iord = iord
        self.ifmt = ifmt
        self.isideopt = isideopt
        self.ialfa = ialfa
        self.nbx = nbx
        self.iprt = iprt
        self.nopt = nopt
        
        self.ngh_file = ngh_file
        self.rtr_file = rtr_file
        self.netcdf_file= netcdf_file
        self.description = description
        self.title = title
        
        self.premod_dylib = '/opt/local/bin/PhD/Ricom/premod.dylib'
        
    
    def run(self):
        
        from ctypes import byref, cdll
        import sys

        if(self.ngh_file == ""):
            print "Please specify an input .NGH file"
            sys.exit()
            
        if(self.rtr_file == ""):
            print "Please specify an input .RTR file"
            sys.exit()

        if(self.netcdf_file == ""):
            print "Please specify an output NetCDF file"  
            sys.exit()
 
        #write the Premod control file (.PMD)
        pmd = open("tsunami.pmd", "w")    
        pmd.write("%s\n" % (self.description))
        pmd.write("%s   %s   %s   %s   %s   %s   !XSCALE,YSCALE,ZSCALE,zoffset,zlimit,izup\n" % \
                  (self.xscale, self.yscale, self.zscale, self.zoffset, self.zlimit, self.izup))
        pmd.write("%s   %s   %s   %s   %s   %s   ! itest,ireorder,ifront,iord,ifmt,isideopt\n" % \
                  (self.itest, self.ireorder, self.ifront, self.iord, self.ifmt, self.isideopt))
        pmd.write("%s\n%s\n%s\n%s\n%s\n%s\n" % (self.ngh_file,self.rtr_file,self.ialfa,self.nbx,self.iprt,self.nopt))
        pmd.write("%s\n" % (self.netcdf_file))
        pmd.close()        
        
        premod = cdll.LoadLibrary(self.premod_dylib)    #call the premod function inside the Fortran .dylib
        premod.premod_() 
        
        del premod
           
        self._write_attributes_()                              
        

    def _write_attributes_(self):
        
        """
        Write the grid attributes to the Netcdf grid ouptutted from premod
        """
        #Open the Netcdf GRID file output from PreMOD
        try: dataset = Dataset(self.netcdf_file,'r+',format='NETCDF4')
        except Exception, e:
            print "ERROR: %s" % e
            sys.exit()

        dataset.title = self.title        
        dataset.description = self.description
        dataset.ngh_file = self.ngh_file
        dataset.rtr_file = self.rtr_file
        dataset.netcdf_file = self.netcdf_file
        dataset.epsg = 4326
        dataset.close()





#SOME USEFUL Functions



def LL2LocalRicom(pointsLL,lat0=0,long0=0,latoff=0,longoff=0):
    """
    Convert from Lat,long to local Ricom Coordinated


    """    
    import math
    bigr = 6378136.0
    dlong = longoff-long0
    dlat = latoff-lat0
    clat0 = math.cos(lat0*math.pi/180.0)
    numPts = len(pointsLL)
    n = 0
    points = []
    while n<numPts:
        x = (pointsLL[n][0]+dlong)*bigr*(math.pi/180.0)*clat0
        y = (pointsLL[n][1]+dlat)*bigr*(math.pi/180.0)
        points.append([x,y])
        n+=1
        
    return points


def ricomLocal2LL(points,lat0=0,long0=0,latoff=0,longoff=0):
    
    """
    Convert from local Ricom Coordinates to Lat,Long


    """ 
    import math
    bigr = 6378136.0
    dlong = longoff-long0
    dlat = latoff-lat0
    clat0 = math.cos(lat0*math.pi/180.0)
    numPts = len(points)
    n = 0
    pointsLL=[]
    while n<numPts:
        x = points[n][0]/(bigr*(math.pi/180.0)*clat0) - dlong
        y = points[n][1]/(bigr*(math.pi/180.0)) - dlat
        pointsLL.append([x,y])
        n+=1
    
    return pointsLL



def convert_points(pointsIN,epsgIN,epsgOUT):
       
    """
    Given a set of points ([x,y,z] or [x,y]) convert from espgIN to self.epsg
    
    """
    
    if(epsgIN != epsgOUT):
        
        coords_in = osr.SpatialReference()
        coords_in.ImportFromEPSG(epsgIN)
        coords_out = osr.SpatialReference() 
        coords_out.ImportFromEPSG(epsgOUT)    
        numPts = len(pointsIN)
        dimension = len(pointsIN[0])
        pointsOUT = []
        n=0
        while n<numPts:
            point = ogr.Geometry(type=ogr.wkbPoint)
            point.SetPoint(0, float(pointsIN[n][0]), float(pointsIN[n][1]))
            point.AssignSpatialReference(coords_in)
            point.TransformTo(coords_out)
            if dimension < 3:
                pointsOUT.append([float(point.GetX()),float(point.GetY())])
            else:
                pointsOUT.append([float(point.GetX()),float(point.GetY()),float(pointsIN[n][2])])
                
            n+=1
        
        return pointsOUT
    
    else:
        return pointsIN
    


def load_run_to_database(run_filename):
    """
    
    
    """
    from Grid_Class import GridClass

    try: dataset = Dataset(run_filename,'r',format='NETCDF4')
    except Exception, e:
        print "ERROR: %s" % e
        sys.exit()
    
    run_id = getattr(dataset,"run_id")          
    grid_dbname = getattr(dataset,"grid_dbname")          
    grid_filename = getattr(dataset,"grid_filename")          
    grid_db_epsg = getattr(dataset,"grid_db_epsg")          
    buildings_dbname = getattr(dataset,"buildings_dbname")          
    name = getattr(dataset,"name")          
    description = getattr(dataset,"description")          
    user = getattr(dataset,"user")          
        
    grid = GridClass(new_project = True, grid_filename=grid_filename,
                grid_dbname=grid_dbname,epsg=grid_db_epsg,buildings_dbname=buildings_dbname, user=user)
    
    
    grid.delete_run(run_id)

    
    grid.add_run(run_id, name, description,
                      run_filename,grid_filename)
    """--------------------------------------------------------------
    Reads the netcdf grid file and returns the side outputs 

    --------------------------------------------------------------"""

    buildings = dataset.groups['buildings']
    building_sides = buildings.groups['sides']

    eta_output_added = getattr(building_sides,'eta_output_added')
    uv_output_added = getattr(building_sides,'uv_output_added')
    eta = []
    uv = []
    sideIds = []
    time = []
    if(eta_output_added or uv_output_added ):
        time = building_sides.variables['time'][:].tolist()
        sideIds = building_sides.variables['id'][:].tolist()
        if eta_output_added: eta = building_sides.variables['eta'][:].tolist()
        if uv_output_added: uv = building_sides.variables['uv'][:].tolist()
    
    i = 0
    print "Adding side output to the database..."
        
    if sideIds != []:
        if (eta != [] and uv != []):
            for side_id in sideIds: 
                grid.add_output_at_side(run_id, side_id, eta[i], uv[i], time)
        elif (eta != [] and uv == []):
            for side_id in sideIds: 
                grid.add_output_at_side(run_id, side_id, eta[i], uv, time)
        elif (eta == [] and uv != []):
            for side_id in sideIds: 
                grid.add_output_at_side(run_id, side_id, eta, uv[i], time)
                
    grid.grid_pg.conn.commit()
   
    
    """--------------------------------------------------------------
    Reads the netcdf grid file and returns the element outputs 

    --------------------------------------------------------------"""

    buildings = dataset.groups['buildings']
    building_elements = buildings.groups['elements']

    eta_output_added = getattr(building_elements,'eta_output_added')
    uv_output_added = getattr(building_elements,'uv_output_added')
    eta = []
    uv = []
    elementIds = []
    time = []
            
    if(eta_output_added or uv_output_added ):
        time = building_elements.variables['time'][:].tolist()
        elementIds = building_elements.variables['id'][:].tolist()
        if eta_output_added: eta = building_elements.variables['eta'][:].tolist()
        if uv_output_added: uv = building_elements.variables['uv'][:].tolist()
    
    i = 0
    print "Adding element output to the database..."
    if elementIds != []:
        if (eta != [] and uv != []):
            for element_id in elementIds: 
                grid.add_output_at_element(run_id, element_id, eta[i], uv[i], time)
        elif (eta != [] and uv == []):
            for element_id in elementIds: 
                grid.add_output_at_element(run_id, element_id, eta[i], uv, time)
        elif (eta == [] and uv != []):
            for element_id in elementIds: 
                grid.add_output_at_element(run_id, element_id, eta, uv, time)
    
    grid.grid_pg.conn.commit()

    """--------------------------------------------------------------
    Reads the netcdf grid file and returns the nodal outputs at the requested building id

    --------------------------------------------------------------"""

    buildings = dataset.groups['buildings']
    building_nodes = buildings.groups['nodes']

    eta_output_added = getattr(building_nodes,'eta_output_added')
    uv_output_added = getattr(building_nodes,'uv_output_added')

    eta = []
    uv = []
    nodeIds = []
    time = []
    
    if(eta_output_added or uv_output_added ):
        time = building_nodes.variables['time'][:].tolist()
        nodeIds = building_nodes.variables['id'][:].tolist()
        if eta_output_added: eta = building_nodes.variables['eta'][:].tolist()
        if uv_output_added: uv = building_nodes.variables['uv'][:].tolist()

    print "Adding node output to the database..."     
    i = 0
    if nodeIds != []:
        if (eta != [] and uv != []):
            for node_id in nodeIds: 
                grid.add_output_at_node(run_id, node_id, eta[i], uv[i], time)
        elif (eta != [] and uv == []):
            for node_id in nodeIds: 
                grid.add_output_at_node(run_id, node_id, eta[i], uv, time)
        elif (eta == [] and uv != []):
            for node_id in nodeIds: 
                grid.add_output_at_node(run_id, node_id, eta, uv, time)


    grid.grid_pg.conn.commit()
