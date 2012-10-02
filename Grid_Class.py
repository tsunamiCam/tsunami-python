#Import the required modules


import psycopg2 #@UnresolvedImport
#import errorcodes 

from netCDF4 import Dataset #@UnresolvedImport
from numpy import *


import sys
from osgeo import ogr #@UnresolvedImport
from osgeo import osr #@UnresolvedImport

#module to allow commands to be sent to the command line (i.e. Terminal)
from subprocess import call


class GridClass:
    
    def __init__(   self, 
                    new_project,
                    grid_filename,
                    grid_dir,
                    grid_dbname,
                    epsg,
                    buildings_dbname,
                    gridEPSG = 4326,
                    user='tbone'):
        '''
        Initialization function for the GridClass object

        Create a new GridClass object
        ------------------------------
        Includes:   1) a GridNC object - interface to the netcdf grid file
                    2) a GridPG object - interface to the PostGIS database which stores the grid
        
        A copy of the footprints table in PostGIS database 'buildings_dbname' is added
        to PostGIS database 'grid_dbname'
        
        '''
        
        self.grid_dbname = grid_dbname
        self.epsg = epsg
        self.user = user
        self.grid_filename = grid_filename
        self.grid_dir = grid_dir
        self.buildings_dbname = buildings_dbname
        self.grid_pg = 0
        self.grid_nc = 0
        
        
        if (new_project == True):
            """
            1) Create new database 'grid_dbname'
            2) Add building footprints from 'buildings_dbname'
            3) Add sides location and indices to Grid file
            4) Add footprints from 'buildings_dbname' to Grid file 
            5) Add nodes, sides and elements to grid_dbname

            
            """
            
            
            # 1) SETUP grid PostGIS database
            print "Tsunami: Creating GridPG object (the grid PostGIS database)"
            self.grid_pg = GridPG(  grid_dbname = self.grid_dbname,
                                    epsg = self.epsg,
                                    user = self.user) 
            
            print "Tsunami: Adding buildings and Areas..."

            if(self.grid_pg.buildingsAdded() == False):
                status = self.grid_pg.add_buildings_from_db(buildings_dbname = self.buildings_dbname,epsgOUT=self.epsg,
                                                       user=self.user)

                if (status == -1 ):
                    print "ERROR: buildings could not be added to grid database. Exiting..."
                    sys.exit()
                else:
                    print "SUCCESS: Buildings added from \'%s\' table" % self.buildings_dbname

            else:
                print "WARNING: Buildings and areas already added to grid database"

            
            
            print "Tsunami: Creating GridNC object (the NetCDF interface to the grid file)"
            # 1) SETUP the grid Netcdf interface object              
            
            print self.grid_dir + "/" + self.grid_filename
            self.grid_nc = GridNC(grid_filename = self.grid_dir  + "/" +  self.grid_filename, epsg = self.epsg)
            
            
            print "Tsunami: Getting the topology of the grid file"

            if(self.grid_pg.domainAdded() == False):
                #GET all the nodes, elements etc. from the Netcdf grid file
                id, xyz, neighbours, ncode = self.grid_nc.get_grid_nodes()
                id, sides = self.grid_nc.get_grid_sides()
                id, elements,ecode = self.grid_nc.get_grid_elements()
                
                
                print "Tsunami: Adding Grid to PostGIS Database"
                
                #Add the grid to the PostGIS database
    
                self.grid_pg.add_domain(xyz=xyz, neighbours=neighbours, ntype=ncode,
                                        elements=elements, etype=ecode, sides=sides)
            
            
            else:
                print "WARNING: Grid already added to PostGIS database"
            
            
        
    def __del__ (self):
        """
        Class deconstructor
    
        """
        cam = 1 #TODO


    def get_building_output_locations(self,area_id,holes=True):
        """
        Get id's of the nodes, elements and sides that are inside the buildings polygons
        
        Input:
            area_id - id of the area polygon that you want to get the info inside
        
        
        """ 
        
        self.grid_pg.init_building_output_tables()

        dictionary = self.grid_pg.get_building_output_locations(area_id,holes)
        return dictionary

    def init_building_output_tables(self):
        """
        
        """
        self.grid_pg.init_building_output_tables()

    def get_closest_elements(self, xy, epsgIN):
        """
        Given a list of points [x,y], return a list of ids corresponding to the element that 
        is closest to each point 

        
        """ 
        elementIds, distanceList = self.grid_pg.get_closest_elements(xy,epsgIN)
        return elementIds, distanceList

    def get_closest_nodes(self, xy,epsgIN):
        """
        Given a list of points [x,y], return a list of ids corresponding to the node that 
        is closest to each point 
        
        """         
        nodeIds, distanceList,elevationList = self.grid_pg.get_closest_nodes(xy,epsgIN)
        return nodeIds, distanceList
    

    def get_grid_bounds(self):
        """
        
        Get the extent of the grid.
        
        """
        return self.grid_pg.get_grid_bounds()
        
        
    def get_area_bounds_by_id(self,area_id):
        """
        
        Given the Id of an area in the PostGIS areas table Get the bounds
        
        """
        bounds = self.grid_pg.get_area_bounds_by_id(area_id)        
        return bounds
        
        
    def add_run(self,run_id, name,description,run_filename, grid_filename):
        """
        Add a run to the database
        
        
        """

        self.grid_pg.add_run(run_id, name, description, run_filename, grid_filename)


    def delete_run(self,run_id):
        
        """
        Remove a run and all its associated data from the database
    
        """
        self.grid_pg.delete_run(run_id)
    
    
    def add_output_at_node(self,run_id, node_id, eta, uv, time):
        """
        Add output at a node
 
        """
        
        self.grid_pg.add_output_at_node(run_id, node_id, eta, uv, time)
        
    
    def add_output_at_element(self,run_id, element_id, eta,uv, time):
        """
        Add output at an element
        
        """
        self.grid_pg.add_output_at_element(run_id, element_id, eta,uv, time)    


    def add_output_at_side(self,run_id, side_id, eta,uv, time):
        """
        Add output at a side
        
        """
        self.grid_pg.add_output_at_side(run_id, side_id, eta,uv, time)    
                
        
    def is_run_defined(self,run_id):
        """
        Check to see if a run is defined in the Database
        
        INPUT: run_id = id of the run
        """
        
        return self.grid_pg.is_run_defined(run_id)  
    
    
    def get_output_at_building(self,run_id,building_id):
        
        self.grid_pg.get_output_at_building(run_id,building_id) 
    
    
    
    def get_elements_in_area(self,area_id = 0):
		"""
		Get all the elements inside a given area.  If area_id = 0, return all elements
		
		area_id - the id (in the postgis table) of the area
		
		
		"""
		return self.grid_pg.get_elements_in_area(area_id)
	
	
	
    def get_elements_inside_buildings(self,area_id = 0):
        """
        Get all the elements inside a given area and inside the building polygons.  If area_id = 0, return all elements
        
        area_id - the id (in the postgis table) of the area
        
        
        """
        return self.grid_pg.get_elements_inside_buildings(area_id)	
	
    
    def get_nodes_at_building_edges(self,area_id = 0):
        """
        Get all the elements inside a given area and inside the building polygons.  If area_id = 0, return all elements
        
        area_id - the id (in the postgis table) of the area
        
        
        """
        return self.grid_pg.get_nodes_at_building_edges(area_id)


    def adjust_ieadj(self):
        '''
        
        TEST - adjust the ieadj array in the grid for a building    
        
        '''
        edgeElements = self.grid_pg.get_elements_at_building_edges()
        ieadj = self.grid_nc.get_grid_ieadj()
        for e1 in edgeElements:
            adjElements = [ieadj[e1-1][1],ieadj[e1-1][2],ieadj[e1-1][3]]
            i = 0
            for e2 in adjElements:
                if e2 in edgeElements:
                    ieadj[e1 - 1][i+1] = 0
                
                i+=1
                    
            print ieadj[e1-1]
                    
        return ieadj
	
	
class GridNC:
    
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

    
    def __init__(self, grid_filename,epsg):
        '''


        Initializing function for GridNC object


        '''

        self.grid_filename = grid_filename
        self.epsg = epsg
        
        #Open the Netcdf GRID file output from PreMOD
        try: self.dataset = Dataset(grid_filename,'r+',format='NETCDF4')
        except Exception, e:
            print "ERROR: %s" % e            

        if(self.grid_filename == ""):
            print "ERROR: No Grid Filename provided"
            #return -1
        
        self.domain = self.dataset.groups['domain']

        try: xyz_utm = self.domain.createVariable(varname = 'xyz_utm',datatype = 'd', dimensions=('np','ncn',)) 
        except Exception, e:
            print "WARNING: %s" % e
        else:
            #store the grid points in coordinate system of the Case Study (i.e. rectanglular coords)
            print "Adding xyz pts in rectangular coords - see xyz_utm"
            xyz = self.domain.variables['xyz'][:]
            xyz = xyz.tolist()
            xyz = convert_points(xyz,4326,self.epsg)
            xyz_utm = xyz
                
        try: self.buildings = self.dataset.createGroup('buildings')
        except Exception, e:
            print "WARNING: %s" % e
            self.buildings = self.dataset.groups['buildings']
        
        try: self.output = self.dataset.createGroup('output')
        except Exception, e:
            print "WARNING: %s" % e            
            self.output = self.dataset.groups['output']       
        try: self.output_nodes = self.output.createGroup('output_nodes')
        except Exception, e:
            print "WARNING: %s" % e           
            self.output_nodes = self.output.groups['output_nodes']
        try: self.output_elements = self.output.createGroup('output_elements')
        except Exception, e:
            print "WARNING: %s" % e            
            self.output_elements = self.output.groups['output_elements']   
        try: self.output_sides = self.output.createGroup('output_sides')
        except Exception, e:
            print "WARNING: %s" % e            
            self.output_sides = self.output.groups['output_sides']
        
        #add the node numbers of the sides [v1, v2] to the netcdf file.  
        self._add_sides_()

        
    def __del__ (self):
        """
        Class deconstructor
    
        """
        self.dataset.close()

    def close(self):
        """
        Close the dataset
    
        """
        self.dataset.close()

    def _add_sides_(self):
        
        '''
        Find the node number for each of the sides in the grid and add it to the grid
        '''
        try: 
            sides_nc = self.domain.createVariable(varname = 'sides', datatype = 'i', dimensions = ('nsides','ntwo',))
            id, sides = self.get_grid_sides()
            sides_nc[:] = sides

        except Exception, e:
            print "WARNING: %s\n" % e
            print "WARNING: Sides variable already created."


    def get_grid_nodes(self):
        """
        Reads the netcdf grid file and returns the nodes - [x y z], node_id and node_type 
    
        """
        xyz = self.domain.variables['xyz'][:]
        xyz = xyz.tolist()
        code = self.domain.variables['nbc'][:]
        code = code.tolist()
        neighbours = self.domain.variables['nbrs'][:]
        neighbours = neighbours.tolist()
        
        i = 1
        id = []
        for pt in xyz:
            id.append(i)
            i += 1

        return id, xyz, neighbours, code    
    
    def get_grid_elements(self):
        """
        Reads the netcdf grid file and returns the elements - [v1 v2 v3], element_id and element_type
    
        """
        
        elements = self.domain.variables['nen'][:]
        elements = elements.tolist()
        code = self.domain.variables['iecode'][:]
        code = code.tolist()
        
        i = 1
        id = []
        for ele in elements:
            id.append(i)
            i += 1
        
        return id, elements, code
        

    def get_grid_ieadj(self):
        """
        Get the elements adjacency array from the Netcdf grid file
    
        """
        ieadj = self.domain.variables['ieadj'][:]
        ieadj.tolist()
        return ieadj
       
    def get_grid_sides(self):
        """
        Reads the netcdf grid file and returns the sides - [v1 v2], side_id
        
        """
       
        elements = self.domain.variables['nen'][:]
        iside = self.domain.variables['iside'][:]
        numside = self.domain.variables['numside'][:]
        
        
        #Determine the node numbers for each side
        id = []
        sides = []     
        nside = 1
        for row in iside:
            en0 = row[0]
            en1 = row[1]
            side = []
            i = 0
            if (en0 == 0): en = (en1-1)
            else: en = (en0-1)
            ele = elements[en]
            while i < 3:
                if numside[en,i] == nside:
                    if i == 2: side = [ele[i],ele[0]]
                    else: side = [ele[i],ele[i+1]] 
                    break    
                i += 1
            sides.append(side)
            id.append(nside)
            nside += 1
                    
        return id, sides
        
        
    def get_output_nodes(self, building_id = 0):
        """
        Reads the netcdf grid file and returns the nodal outputs at the requested building id
    
        """
        
        import matplotlib.pyplot as plt
        
        eta = self.output_nodes.variables['eta'][:]
        time = self.output_nodes.variables['time'][:]
        uv = self.output_nodes.variables['uv'][:]
        node_id = self.output_nodes.variables['id'][:]
        footprint_id = self.buildings.variables['footprint_id'][:]
        footprint_nodes = self.buildings.variables['footprint_nodes'][:]
        
        dict_id = {}        #dictionary maps the the node id to the location in the id array
        n = 0
        for id in node_id:
            dict_id[id] = n 
            n += 1
                    
        dict_nodes = {}     #dictionary maps the footprint id to the nodal output from the netcdf file
        dict_nodes['time'] = time.tolist()
        n=0
        for id in footprint_id:
            dict_nodes[id] = {'node_ids': [], 'eta': {}, 'uv': {}}
            nodesList = []
            for node in footprint_nodes[n]:
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
    
    def get_output_elements(self, building_id = 0):
    
        """
        Reads the netcdf grid file and returns the element outputs at the requested building id
    
        """
        eta = []
        uv = []
        
        return eta, uv
    
    def get_output_sides(self, building_id = 0):
    
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
        footprint_ids = self.buildings.variables['footprint_id'][:]
        node_id = self.output_nodes.variables['id'][:]
        time = self.output_nodes.variables['time'][:]
        footprint_nodes = self.buildings.variables['footprint_nodes'][:]
        xyz = self.domain.variables['xyz'][:]
        eta_all = self.output_nodes.variables['eta'][:]
        uv_all = self.output_nodes.variables['uv'][:]
        
        n = 0
        for id in node_id:
            dict_ids['nodes'][id] = n 
            n += 1
            
        n = 0
        for id in footprint_ids:
            dict_ids['footprints'][id] = n 
            n += 1
        
       
        i_footprint = 0
        for id in footprint_ids:            #iterate throught the footprints
            print "Footprint Id = %s\n" % str(id)
            nodes = footprint_nodes[i_footprint]
            i_time = 0
            maxData = []
            for step in time:
                fd_max = 0
                speed_max = 0
                for node in nodes:          #For current footprint: get the max node data for the time step
                    if node != 0:
                        #eta = self.output_nodes.variables['eta'][dict_ids['nodes'][node]][i_time]
                        eta = eta_all[dict_ids['nodes'][node]][i_time]
                        uv = uv_all[dict_ids['nodes'][node]][i_time]
                        #uv = self.output_nodes.variables['uv'][dict_ids['nodes'][node]][i_time]
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
        footprint_ids = self.buildings.variables['footprint_id'][:]
        node_id = self.output_nodes.variables['id'][:]
        element_id = self.output_elements.variables['id'][:]
        side_id = self.output_sides.variables['id'][:]

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
        for id in footprint_ids:
            dict_ids['footprints'][id] = n 
            n += 1
        
                
        nodes = self.buildings.variables['footprint_nodes'][dict_ids['footprints'][building_id]]
        elements = self.buildings.variables['footprint_elements'][dict_ids['footprints'][building_id]]
        sides = self.buildings.variables['footprint_sides'][dict_ids['footprints'][building_id]]
            
        for n in nodes:
            if n != 0:
                xyz = self.domain.variables['xyz'][n-1]
                eta = self.output_nodes.variables['eta'][dict_ids['nodes'][n]]
                i = 0
                while (i < len(eta)):    
                    if eta[i] - xyz[2] < 0.01: eta[i] = 0
                    else: eta[i] = eta[i] - xyz[2] 
                    i+=1
                    
                uv = self.output_nodes.variables['uv'][dict_ids['nodes'][n]]
            else: break            

        for e in elements:
            if e != 0:
                element_nodes = self.domain.variables['nen'][e-1]
                eta = self.output_elements.variables['eta'][dict_ids['elements'][e]]
                uv = self.output_elements.variables['uv'][dict_ids['elements'][e]]
            else: break 

        for s in sides:
            if s != 0:
                side_nodes = self.domain.variables['sides'][s-1]        #get the node numbers of the side
                xyz0 = self.domain.variables['xyz'][side_nodes[0]-1]    #get the node coordinates for node0 of the side
                xyz1 = self.domain.variables['xyz'][side_nodes[1]-1]    #get the node coordinates for node1 of the side
                eta = self.output_elements.variables['eta'][dict_ids['sides'][s]]
                uv = self.output_elements.variables['uv'][dict_ids['sides'][s]]
            else: break 
  
    
    def add_buildings(self,buildings):
        '''
        
        '''
        
        
        add = 1
    
    def add_output_from_dictionary(self, dictionary = {}):
        """
        Given a dictionary of building footprints and associated nodes,element and sides, add the values 
        to the netcdf grid file.
        
        The nodes, elements and sides associated with each footprint correspond to the there index in the RiCOM grid file
        
        Dictionary format:
        {id1: {'nodes': [n1, n2,...nn] }, {'elements': [e1,e2,...,en] },{'sides': [s1,s2,...,sn]}, id2: {}, id3 {}, ...., idn {} } 
        
        idn = the id of the building footprint that the node, elements and sides belong to
        
        """

        maxNodes = 0
        maxElements = 0
        maxSides = 0
        nodesAll = []
        elementsAll = []
        sidesAll = []
        id = []
        
        for row in dictionary.iteritems(): 
            id.append(row[0])              
            n = row[1]['nodes'] 
            e = row[1]['elements']
            s = row[1]['sides']
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
        try: self.output_nodes.createDimension('number_of_nodes',len(nodesAll))
        except Exception, e: print "WARNING: %s" % e
        try: self.output_elements.createDimension('number_of_elements',len(elementsAll))
        except Exception, e: print "WARNING: %s" % e
        try: self.output_sides.createDimension('number_of_sides',len(sidesAll))
        except Exception, e: print "WARNING: %s" % e
        
        
        #create variables
        try: footprint_id = self.buildings.createVariable(varname = 'footprint_id',datatype = 'i', dimensions=('number_of_buildings',)) 
        except Exception, e:
            footprint_id = self.buildings.variables['footprint_id']
            print "WARNING: %s" % e

        try: footprint_wkt = self.buildings.createVariable(varname = 'footprint_wkt',datatype = str, dimensions=('number_of_buildings',)) 
        except Exception, e:
            footprint_wkt = self.buildings.variables['footprint_wkt']            
            print "WARNING: %s" % e

        try: footprint_nodes = self.buildings.createVariable(varname = 'footprint_nodes',datatype = 'i', dimensions=('number_of_buildings','max_number_nodes',)) 
        except Exception, e:
            footprint_nodes = self.buildings.variables['footprint_nodes']            
            print "WARNING: %s" % e

        try: footprint_elements = self.buildings.createVariable(varname = 'footprint_elements',datatype = 'i', dimensions=('number_of_buildings','max_number_elements',)) 
        except Exception, e:
            footprint_elements = self.buildings.variables['footprint_elements']
            print "WARNING: %s" % e
        
        try: footprint_sides = self.buildings.createVariable(varname = 'footprint_sides',datatype = 'i', dimensions=('number_of_buildings','max_number_sides',)) 
        except Exception, e:
            footprint_sides = self.buildings.variables['footprint_sides']
            print "WARNING: %s" % e
        
        footprint_nodes[:] = nodes
        footprint_elements[:] = elements
        footprint_sides[:] = sides
        footprint_id[:] = array(id) 
        
        #Set the attributes
        self.output_nodes.start = 3200
        self.output_nodes.finish = 4800
        self.output_nodes.step = 1
        self.output_elements.start = 3200
        self.output_elements.finish = 4800
        self.output_elements.step = 1
        self.output_sides.start = 3200
        self.output_sides.finish = 4800
        self.output_sides.step = 1
        
        #assign the data
        output_ids = {'nodes': [], 'elements': [], 'sides': []}
        try: output_ids['nodes'] = self.output_nodes.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_nodes',))
        except Exception, e:
            output_ids['nodes'] = self.output_nodes.variables['id']
            print "WARNING: %s" % e
        try: output_ids['elements'] = self.output_elements.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_elements',))
        except Exception, e:
            output_ids['elements'] = self.output_elements.variables['id']
            print "WARNING: %s" % e
        try: output_ids['sides'] = self.output_sides.createVariable(varname = 'id',datatype = 'i', dimensions=('number_of_sides',))
        except Exception, e:
            output_ids['sides'] = self.output_sides.variables['id']
            print "WARNING: %s" % e


        output_ids['nodes'][:] = array(nodesAll)
        output_ids['elements'][:] = array(elementsAll)
        output_ids['sides'][:] =  array(sidesAll)
        
        print self.output_nodes.start



class GridPG:
    """
    class GridPG
    
    This class opens and manages an existing PostGIS database to store the my Tsunami Project data
    
    """
    def __init__(self, grid_dbname, epsg = 4326,user = "tbone"):
    
    
        '''
        
        INIT function for GridPG Class
        
        '''
        
        self.grid_dbname = grid_dbname            #the database name
        self.user = user    
        self.epsg = epsg
        self.conn = 0
        self.cur = 0    
        
        #PostGIS SQL script locations
        pg_dir = '/opt/local/share/postgresql90/contrib/postgis-1.5/'
        pg_1 = pg_dir + 'postgis.sql'
        pg_2 = pg_dir + 'spatial_ref_sys.sql'
        pg_3 = pg_dir + 'postgis_comments.sql'
        
        try:           #open connection to the requested database     
            conn = psycopg2.connect(database=self.grid_dbname, user=self.user);        
        except Exception, e:
            if e[0].find("could not connect to server") >= 0:
                print "Postgis database server not running. Please start postgres\n"
                sys.exit()   
        
            elif e[0].find("does not exist") >= 0:
                print "Creating new GRID PostGIS database: \"%s\"" % self.grid_dbname            
                p = call(['createdb', self.grid_dbname])
                print p
                p = call(['psql','-d',self.grid_dbname,'-f',pg_1,'-q'])
                print p
                p = call(['psql','-q','-d',self.grid_dbname,'-f',pg_2,'-q'])
                print p
                p = call(['psql','-q','-d',self.grid_dbname,'-f',pg_3,'-q'])  
                print p
                #Try connecting to the database again
                try:
                    conn = psycopg2.connect(database=self.grid_dbname, user=self.user);
                except:
                    "ERROR: cannot connect or create database: \"%s\"" % self.grid_dbname 
                    sys.exit()
                else: 
                    cur = conn.cursor()
                    self.cur = cur
                    self.conn = conn                    
                    self.cur.execute("CREATE TABLE nodes (id int4, xyz float[3], code int4);")
                    self.cur.execute("SELECT AddGeometryColumn('nodes', 'geom', %s, 'POINT', 2);" % self.epsg)    
                    self.cur.execute("CREATE TABLE elements (id int4, node_ids int4[3], code int4);")    
                    self.cur.execute("SELECT AddGeometryColumn('elements', 'geom', %s, 'POLYGON', 2);" % self.epsg)    
                    self.cur.execute("CREATE TABLE sides (id int4, node_ids int4[2]);")    
                    self.cur.execute("SELECT AddGeometryColumn('sides', 'geom', %s, 'LINESTRING', 2);" % self.epsg)  
                    self.buildings_table_created = False
                    self.cur.execute("CREATE TABLE buildings (pkey integer PRIMARY KEY,id int4);")     
                    self.cur.execute("SELECT AddGeometryColumn('buildings', 'geom', %s, 'POLYGON', 2);" % self.epsg)
                    #ADD extra PTVA columns to the 'buildings' table
                    self.cur.execute("ALTER TABLE buildings ADD COLUMN s float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN m float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN g float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN f float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN so float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN mo float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN pc float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN prot_br float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN prot_nb float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN prot_sw float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN prot_w float DEFAULT 0; \
                                    ALTER TABLE buildings ADD COLUMN comment varchar(128); \
                                    ALTER TABLE buildings ADD COLUMN nodes integer[]; \
                                    ALTER TABLE buildings ADD COLUMN elements integer[]; \
                                    ALTER TABLE buildings ADD COLUMN sides integer[];")            
                    
                    #CREATE RULE for when new records are added to the buildings table
                    #This rules makes the ID of new building record equal to PKEY
                    self.cur.execute("CREATE RULE force_id AS ON INSERT TO buildings \
                                    DO ALSO \
                                    UPDATE buildings SET id = NEW.pkey WHERE pkey = NEW.pkey;")  


                    self.cur.execute("CREATE TABLE ptva (building_id int4, s int4 DEFAULT 0, m int4 DEFAULT 0, g int4 DEFAULT 0, \
                                                    f int4 DEFAULT 0, so int4 DEFAULT 0, mo int4 DEFAULT 0, pc int4 DEFAULT 0, \
                                                    prot_br int4 DEFAULT 0, prot_nb int4 DEFAULT 0, prot_sw int4 DEFAULT 0, prot_w int4 DEFAULT 0);")
                    
                
                    self.cur.execute("CREATE TABLE lines (pkey integer PRIMARY KEY,id int4);")
                    self.cur.execute("SELECT AddGeometryColumn('lines', 'geom', %s, 'LINESTRING', 2);" % self.epsg)
                    self.cur.execute("ALTER TABLE lines ADD COLUMN name varchar(50);")
                    self.cur.execute("ALTER TABLE lines ADD COLUMN comment varchar(128);")
                    self.cur.execute("CREATE RULE force_id AS ON INSERT TO lines DO ALSO UPDATE lines SET id = NEW.pkey WHERE pkey = NEW.pkey;")
            
                    self.cur.execute("CREATE TABLE areas (pkey integer PRIMARY KEY,id int4);")
                    self.cur.execute("SELECT AddGeometryColumn('areas', 'geom', %s, 'POLYGON', 2);" % self.epsg)
                    self.cur.execute("ALTER TABLE areas ADD COLUMN name varchar(50);")
                    self.cur.execute("ALTER TABLE areas ADD COLUMN comment varchar(128);")
                    self.cur.execute("CREATE RULE force_id AS ON INSERT TO areas DO ALSO UPDATE areas SET id = NEW.pkey WHERE pkey = NEW.pkey;")  
                    
                    
                    #create table to store the Runs for this grid
                    self.cur.execute("CREATE TABLE runs (pkey integer PRIMARY KEY, id int4, name varchar(50), description varchar(128), \
                                         run_filename varchar(50), grid_filename varchar(50));")

                    #create table to store outputs of the runs
                    self.cur.execute("CREATE TABLE output_nodes (node_id int4, run_id int4, eta float[], uv float[][], time float[]);")
                    self.cur.execute("CREATE TABLE output_elements (element_id int4, run_id int4, eta float[], uv float[][], time float[]);")
                    self.cur.execute("CREATE TABLE output_sides (element_id int4, run_id int4, eta float[], uv float[][], time float[]);")
                    
                    
                    self.conn.commit()

                    self.buildings_table_created = True 
        else:
            cur = conn.cursor()
            self.cur = cur
            self.conn = conn
                                           
    
            
    def __del__ (self):
        """
        Class deconstructor
        
        Closes the connection to the PostGIS database
        """
        
        self.cur.close()       
        self.conn.close()
        
    
    def buildingsAdded(self):
        '''
        Returns true if building have been added to the buildings table already
        '''
        
        self.cur.execute("SELECT * FROM buildings LIMIT 1;" )
        if (len(self.cur.fetchall()) == 0):   
            return False     
        else:
            return True
    
    def domainAdded(self):
        '''
        Returns true if domain have been added to the database
        '''
        
        self.cur.execute("SELECT * FROM nodes LIMIT 1;" )
        a = len(self.cur.fetchall())
        self.cur.execute("SELECT * FROM elements LIMIT 1;" )
        b = len(self.cur.fetchall())
        self.cur.execute("SELECT * FROM sides LIMIT 1;" )
        c = len(self.cur.fetchall())


        if ((a == 0) or (b == 0) or (c==0)):   

#            self.cur.execute("DROP TABLE IF EXISTS nodes;")
#            self.cur.execute("DROP TABLE IF EXISTS elements;")
#            self.cur.execute("DROP TABLE IF EXISTS sides;")
#
#            self.cur.execute("CREATE TABLE nodes (id int4, xyz float[3], code int4);")
#            self.cur.execute("SELECT AddGeometryColumn('nodes', 'geom', %s, 'POINT', 2);" % self.epsg)    
#            self.cur.execute("CREATE TABLE elements (id int4, node_ids int4[3], code int4);")    
#            self.cur.execute("SELECT AddGeometryColumn('elements', 'geom', %s, 'POLYGON', 2);" % self.epsg)    
#            self.cur.execute("CREATE TABLE sides (id int4, node_ids int4[2]);")    
#            self.cur.execute("SELECT AddGeometryColumn('sides', 'geom', %s, 'LINESTRING', 2);" % self.epsg)  

            return False
        
         
        else:
            return True

    def get_grid_bounds(self):
        
        """
        
        Get the bounds of the grid
        
        """

        self.cur.execute("SELECT ST_Extent(geom) FROM elements;")
        bounds = self.cur.fetchall() 
        return bounds[0][0]
    
    def _drop_buildings_table_(self):
        '''
        Drop the buildings table
        '''
        self.cur.execute("DROP TABLE IF EXISTS buildings;")
    
    def add_building_results_table(self,new_tablename):
        """
        Make a copy of the base buildings table to which results of a model run are added to.
        """
        

        
        
        self.cur.execute("DROP TABLE IF EXISTS %s;" % new_tablename)
        self.cur.execute("CREATE TABLE %s (id int4);" % new_tablename)     
        self.cur.execute("SELECT AddGeometryColumn('%s', 'geom', %s, 'POLYGON', 2);" % (new_tablename,self.epsg))
    
        self.cur.execute("ALTER TABLE %s ADD COLUMN s float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN m float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN g float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN f float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN so float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN mo float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN pc float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN prot_br float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN prot_nb float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN prot_sw float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN prot_w float DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN damage int4 DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN comment varchar(128); \
                ALTER TABLE %s ADD COLUMN visibility_before int4 DEFAULT 0; \
                ALTER TABLE %s ADD COLUMN visibility_after int4 DEFAULT 0; \
                ALTER TABLE %s ADD wall_height float DEFAULT 0;" % (new_tablename,new_tablename,new_tablename,new_tablename,new_tablename,new_tablename,new_tablename,
                                                                    new_tablename,new_tablename,new_tablename,new_tablename,new_tablename,new_tablename,new_tablename,
                                                                    new_tablename,new_tablename))  
        

        self.cur.execute("ALTER TABLE %s ADD COLUMN etamax_nodes float;" % new_tablename)
        self.cur.execute("ALTER TABLE %s ADD COLUMN fdmax_nodes float;" % new_tablename)
        self.cur.execute("ALTER TABLE %s ADD COLUMN speedmax_nodes float;" % new_tablename)        
        self.cur.execute("ALTER TABLE %s ADD COLUMN exposure float;" % new_tablename)
        self.cur.execute("ALTER TABLE %s ADD COLUMN fragility float;" % new_tablename)
        self.cur.execute("ALTER TABLE %s ADD COLUMN z float;" % new_tablename)

        self.conn.commit()
    
        self.cur.execute("SELECT id,geom FROM buildings;") 
        buildings = self.cur.fetchall() 
        
         
        #get the PTVA attributes for the Yuriage Google buildings table
        try:
            connTEMP = psycopg2.connect(database='yuriage_google', user=self.user)
        
        except Exception, e:
            print e
            sys.exit()
              
        curTEMP = connTEMP.cursor()
        
        for b in buildings:
            curTEMP.execute("SELECT s,m,g,f,so,mo,pc,prot_br,prot_nb,prot_sw, prot_w,damage,\
                         visibility_before,visibility_after,wall_height FROM buildings WHERE id = %s;" % (b[0]))            
            ptva = curTEMP.fetchall() 

            p = ptva[0]

            print len(p)
            self.cur.execute("INSERT INTO %s (id,geom,s,m,g,f,so,mo,pc,prot_br,prot_nb,prot_sw, prot_w,damage,\
                         visibility_before,visibility_after,wall_height) \
                         VALUES ( \
                         %s,ST_Geometry('%s'), \
                         %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);" % (new_tablename,b[0],b[1],p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14]))
        
        self.conn.commit()

        
        curTEMP.close()
        connTEMP.close()   
            
    def add_buildings_from_db(self,buildings_dbname, epsgOUT, user='tbone'):
        """
        Copy the building table from 'buildings_dbname' to this grid database
        
        Returns: True if buildings are sucessfull added to the grid_db
                 False if buildings have already been added
        """    

       	print "BUILDING dbname = %s" %  buildings_dbname   
        if (self.buildingsAdded() == False):
            print "ADDING buildings table to database = \"%s\"" % self.grid_dbname
             #open connection to the buildings (buildings_dbname) database 
            try:              
                connTEMP = psycopg2.connect(database=buildings_dbname, user=user);
        
            except Exception, e:
                if e[0].find("could not connect to server") >= 0:
                    print "Postgis database server not running. Please start postgres\n"
                    return -1   
        
                elif e[0].find("does not exist") >= 0:
                    print "ERROR: Database \"%s\" doesn't exist. Exiting..." % buildings_dbname            
                    return -1

            curTEMP = connTEMP.cursor()
        
            #ST_AsText
            try:
                curTEMP.execute("SELECT pkey, id, ST_Transform(geom,%s) FROM buildings;" % (epsgOUT))
            except Exception, e:
                print "ERROR: buildings table not available in \"%s\"" % buildings_dbname 
                return -1
            else:            
                buildings = curTEMP.fetchall() 
         
                        
                sql = ""
                #ST_GeomFromText('%s',%s)
                for b in buildings:
                    sql = sql + "INSERT INTO buildings (pkey, id, geom) \
                                  VALUES (%s,%s,ST_Geometry('%s'));" % (b[0],b[1],b[2])
                
                if(sql != ""):
                    self.cur.execute(sql)
                    self.conn.commit()
            
            
            print "ADDING areas table to database = \"%s\"" % self.grid_dbname
            try:
                curTEMP.execute("SELECT pkey, id, ST_Transform(geom,%s), name, comment FROM areas;" % (epsgOUT))
            except Exception, e:
            	print e
                print "***ERROR: areas table not available in \"%s\"" % buildings_dbname
                return -1
            else:            
                areas = curTEMP.fetchall()            
                sql = ""
                #ST_GeomFromText('%s',%s)
                for b in areas:
                    sql = sql + "INSERT INTO areas (pkey, id, geom, name, comment) \
                                  VALUES (%s,%s,ST_Geometry('%s'),'%s','%s');" % (b[0],b[1],b[2],b[3],b[4])
                
                
                if(sql != ""):
                    self.cur.execute(sql)
                    self.conn.commit()
         
            print "ADDING lines table to database = \"%s\"" % self.grid_dbname
            
            try:
                curTEMP.execute("SELECT pkey, id, ST_Transform(geom,%s), name, comment FROM lines;" % (epsgOUT))
            except Exception, e:
                print e[0]
                print "ERROR: line table not available in \"%s\"" % buildings_dbname
                return -1
            else:            
                areas = curTEMP.fetchall()            
                sql = ""
                #ST_GeomFromText('%s',%s)
                for b in areas:
                    sql = sql + "INSERT INTO lines (pkey, id, geom, name, comment) \
                                  VALUES (%s,%s,ST_Geometry('%s'),'%s','%s');" % (b[0],b[1],b[2],b[3],b[4])
                
                if(sql != ""):
                    self.cur.execute(sql)
                    self.conn.commit()
                
            curTEMP.close()
            connTEMP.close()   

            return 1
        
        else:
            return -1        
        

    def add_domain(self, xyz, neighbours, ntype, elements, etype, sides, epsgIN = 4326):
        """
        Add a RiCOM domain (nodes, elements and sides) to the PostGIS table
        
        INPUT:  xyz - the x,y and z coordinates of the node points
                neighbours - the neighbouring nodes for each node point
                elements - the node points (node#) corresponding to each element [v1 v2 v3]
                sides - the node points (node #) corresponding to each side [v1 v2]
        
       
        """      
        if (self.domainAdded() == False):
            print "GridPG: Adding domain to PostGIS database."
            nNodes = len(xyz)
            nEle = len(elements)
            nSides = len(sides)
    
            nNeighs = len(neighbours[0])
            #some sanity checks
            if (nNodes != len(neighbours)) or (nNodes != len(ntype)):
                print "The length of the xyz and neighbours list does not match!"
                return 0
    
            if nEle != len(etype):
                print "The length of the elements and etype lists do not match!"
                return 0
            
            try: self.cur.execute("ALTER TABLE nodes ADD COLUMN neighbours float[%s]" % (nNeighs))
            except:
                print "Can't add neighbour column to the nodes table"
                return 0
            
            self.conn.commit()
    
            #convert the input points (xyz) to the Database spatial reference system (i.e. self.epsg)
            xyz = self._convert_points_(xyz, epsgIN=epsgIN,epsgOUT=self.epsg)
            
            #add the nodal data to the node table
            i = 0
            
            print "Adding nodes..."

            while i < nNodes:
                sql1 =  "INSERT INTO nodes (id, xyz, code, geom, neighbours)"
                sql2 = "VALUES (%s, ARRAY%s, %s, ST_GeomFromText('POINT(%s %s)',%s), ARRAY%s);" % (i+1, xyz[i], ntype[i], xyz[i][0], xyz[i][1],self.epsg, neighbours[i])
                sql = sql1 + sql2
                self.cur.execute(sql)
                i += 1
                

            #Create spatial index for the node table
            self.cur.execute("CREATE INDEX nodes_idx ON nodes USING GIST (geom);")

            '''
            print "Getting building elevations..."

            self.cur.execute("SELECT id, geom,ST_Perimeter(geom),ST_X(ST_Centroid(geom)),ST_Y(ST_Centroid(geom)) FROM buildings;")
            buildings = self.cur.fetchall()
            
            buildingsCentroid = []
            buildingIds = []
            for b in buildings:
                
                buildingsCentroid.append([float(b[3]),float(b[4])])
                buildingIds.append(b[0])
                
            #get the elevations of the nodes closest to each of the building footprint polygons
            ids,distances,elevations = self.get_closest_nodes(buildingsCentroid, self.epsg) 
                        
            self.cur.execute("ALTER TABLE buildings ADD COLUMN elevation float;")
            self.conn.commit()
            

            i = 0
            for id in buildingIds:
                self.cur.execute("UPDATE buildings SET elevation = %s WHERE id = %s" % (elevations[i],id))
                i+=1
    
            '''
            
            print "Adding elements..."

            #add the elements to the database
            i = 0
            while i < nEle:     
                v1 = xyz[elements[i][0]-1]
                v2 = xyz[elements[i][1]-1]
                v3 = xyz[elements[i][2]-1]
                sql1 =  "INSERT INTO elements (id, node_ids, code, geom) VALUES (%s, ARRAY%s, %s," % (i+1, elements[i],etype[i]) 
                sql2 = "ST_GeomFromText('POLYGON((%s %s,%s %s,%s %s,%s %s))',%s));" % (v1[0],v1[1],v2[0],v2[1],v3[0],v3[1],v1[0],v1[1],self.epsg)
                sql = sql1 + sql2
                self.cur.execute(sql)
                i += 1
            #Create spatial index for the node table
            self.cur.execute("CREATE INDEX elements_idx ON elements USING GIST (geom);")
            
            
            print "Adding sides..."

            #add the sides to the database
            i = 0
            while i < nSides:     
                v1 = xyz[sides[i][0]-1]
                v2 = xyz[sides[i][1]-1]
                sql1 =  "INSERT INTO sides (id, node_ids, geom) VALUES (%s, ARRAY%s," % (i+1, sides[i]) 
                sql2 = "ST_GeomFromText('LINESTRING(%s %s,%s %s)',%s));" % (v1[0],v1[1],v2[0],v2[1],self.epsg)
                sql = sql1 + sql2
                self.cur.execute(sql)
                i += 1
                4    


                
            #Create spatial index for the node table
            self.cur.execute("CREATE INDEX sides_idx ON sides USING GIST (geom);")
            self.conn.commit()
            return 1
        
        else:
            return -1
        
        
    def _convert_points_(self,pointsIN,epsgIN,epsgOUT):
           
        """
        Given a set of points ([x,y,z] or [x,y]) convert from espgIN to self.epsg
        PRIVATE function to the PGTsunami CLASS
        
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
 
 
    def _offset_polygon2_(self, polygon_AsText,offset=0.1):
        """
        Offset the passed polygon offset (in meters)
        PRIVATE function to the PGTsunami CLASS
        
        """
        #Assumes that the input polygon is closed (i.e. the endpoint is the same as the start)
        #Assumes that input polygon is Clockwise defines (is this the standard??)
        #polyIn is a list of x,y points
      
      
        import math
        id = []
        poly = polygon_AsText
        polygon = []

        i = poly.find('POLYGON((')
        if (i < 0):
            print "NCGrid: Invalid polygon"
            return -1
        
        j = poly.find('))')
        if (j < 0):
            print "NCGrid: Invalid polygon"
            return -1
        
        points = poly[i+9:j]
        points = points.replace(',',' ').split()
    
        i = 0
        while i < len(points):
            polygon.append([float(points[i]), float(points[i+1])])
            i += 2
       
        poly = polygon            
        numPts = len(poly)
        polyOffset = []
        n = 0
        while n < numPts:
            #calculate the angle between segment [ n, n+1] and segment [n+1, n+2]
            if n == 0:
                i0 = numPts-2
                i1 = n
                i2 = n+1          
            elif n == numPts-1:
                i0 = n-1
                i1 = n
                i2 = 1  
            else:
                i0 = n-1
                i1 = n
                i2 = n+1         
                
            #translate vector to the origin (the point at index = n)
            d1 = [(poly[i0][0] - poly[i1][0]),(poly[i0][1] - poly[i1][1])]
            d2 = [(poly[i2][0] - poly[i1][0]),(poly[i2][1] - poly[i1][1])]
            #make all vectors unit vectors
            m1 = math.sqrt(d1[0]*d1[0]+d1[1]*d1[1])
            m2 = math.sqrt(d2[0]*d2[0]+d2[1]*d2[1])
            u1 = [d1[0]/m1, d1[1]/m1]
            u2 = [d2[0]/m2, d2[1]/m2]
            #calculate the angle between the two vectors
            a = (math.atan2(d2[1], d2[0])-math.atan2(d1[1], d1[0]))
            if a > 0:
                a = (-2*math.pi+a)/2
            else:
                a = a/2
            
            #rotate the the first unit vector (point n-1) to half the angle between the two unit vectors
            R = [[math.cos(a), -math.sin(a)],[math.sin(a),math.cos(a)]]
            nx = u1[0]*R[0][0] + u1[1]*R[0][1]
            ny = u1[0]*R[1][0] + u1[1]*R[1][1]
    
            polyOffset.append([poly[n][0]+nx*offset,poly[n][1]+ny*offset])            
            n+=1    
        
        poly = polyOffset
        wkt = "POLYGON(("
        i = 0
        nPts = len(poly)
        while i < nPts - 1:
            ptStr = "%s %s," % (poly[i][0], poly[i][1])
            wkt = wkt + ptStr
            i += 1    
        ptStr = "%s %s))" % (poly[i][0], poly[i][1])
        wkt = wkt + ptStr
        
        return wkt
 
 
            
            
    def geom_to_kml(self, geomList=[],kml_file="PostGIS.kml", epsgIN=28355):
        """
        Writes a list of arbitrary PostGIS geometries to a KML file
        ---------------------------------------------------------------------------------
        INPUT:     A list of PostGIS geometry entities as strings (i.e. POINT, POLYGON, LINESTRING etc.).
        OUTPUT:    Output KML file
        ---------------------------------------------------------------------------------
        """
        if(geomList != []):    
            
            polyList = []
            pointsList = []
            lineStrList = []
            
            #parse a list of SQL geometries and extract the values (to be output as a KML)
            for geom in geomList:
                #extract the polygon verices from postgis polygon    
                if (geom.find('POLYGON') >= 0):
                    i = geom.find('POLYGON((')
                    j = geom.find('))')
                    pts = geom[i+9:j]
                    pts = pts.replace(',',' ').split()
                    i = 0
                    polygon = []
                    while i < len(pts):
                        polygon.append([float(pts[i]), float(pts[i+1])])
                        i += 2      
                    polyList.append(polygon)
                        
                elif (geom.find('POINT') >= 0):
                    i = geom.find('POINT(')
                    j = geom.find(')')
                    pt = geom[i+6:j].split()
                    #pt = pt.replace(',',' ').split()
                    pointsList.append([float(pt[0]),float(pt[1])])
                
                elif (geom.find('LINESTRING') >= 0):
                    i = geom.find('LINESTRING(')
                    j = geom.find(')')
                    pts = geom[i+11:j]
                    pts = pts.replace(',',' ').split()
                    lineStr = []
                    i = 0
                    while i < len(pts):
                        lineStr.append([float(pts[i]), float(pts[i+1])])
                        i += 2      
                    lineStrList.append(lineStr)                
                else:
                    print "Unhandled SQL geometry..."
        
            # set the spatial reference - of the polygon data
            coordsLL = osr.SpatialReference()
            coordsLL.ImportFromEPSG(4326)
            coordsIN = osr.SpatialReference()
            coordsIN.ImportFromEPSG(int(epsgIN))
        
            #write the polygon data to KML
            driver = ogr.GetDriverByName('KML')
            kmlData = driver.CreateDataSource(kml_file)
            layerPolygons =  kmlData.CreateLayer("polygons", coordsLL, ogr.wkbLineString)
            layerPoints =  kmlData.CreateLayer("points", coordsLL, ogr.wkbPoint)
            layerLines =  kmlData.CreateLayer("lines", coordsLL, ogr.wkbLineString)
        
            featurePolygons = ogr.Feature(layerPolygons.GetLayerDefn())
            featurePoints = ogr.Feature(layerPoints.GetLayerDefn())
            featureLines = ogr.Feature(layerLines.GetLayerDefn())
            
            #write the polygon data to the KML file
            for poly in polyList:
                line = ogr.Geometry(type=ogr.wkbLineString)
                for pt in poly:
                    line.AddPoint(pt[0],pt[1])
                      
                line.AssignSpatialReference(coordsIN)
                line.TransformTo(coordsLL)
                featurePolygons.SetGeometry(line)
                layerPolygons.CreateFeature(featurePolygons)
        
            #write the point data to the KML file
            for pt in pointsList:
                point = ogr.Geometry(type=ogr.wkbPoint) 
                point.SetPoint(0, pt[0], pt[1])
                point.AssignSpatialReference(coordsIN)
                point.TransformTo(coordsLL)
                featurePoints.SetGeometry(point)
                layerPoints.CreateFeature(featurePoints)
                
            #write the polygon data to the KML file
            for lstr in lineStrList:
                line = ogr.Geometry(type=ogr.wkbLineString)
                for pt in lstr:
                    line.AddPoint(pt[0],pt[1])
                      
                line.AssignSpatialReference(coordsIN)
                line.TransformTo(coordsLL)
                featureLines.SetGeometry(line)
                layerLines.CreateFeature(featureLines)
            
            kmlData.Destroy()
            featurePolygons.Destroy()
            featurePoints.Destroy()
            featureLines.Destroy()
            
        else:
            print "No point in PostGIS Geometry List!!"
        
    
    def get_area_by_id(self,area_id):
        '''
        Given an area id (i.e. polygon defined in the areas table)
        
        Returns: the corresponding area geometry (polygon) as an ST_Geometry object (i.e. PostGIS)
        
        '''
        try: self.cur.execute("SELECT ST_AsText(geom) FROM areas WHERE id = %s;" % str(area_id))
        except Exception, e:
            print e[0]
            print "ERROR: areas table not available"
            return ""
        else:
            area = self.cur.fetchall()
            if (len(area) != 0):
                return area[0][0] 
            else:
                print "ERROR: Requested area id not available"
                return ""

    def get_area_bounds_by_id(self,area_id):
        '''
        Given an area id (i.e. polygon defined in the areas table)
        
        Returns: The bounds of the area geometry as an ST_Geometry
        
        '''
        try: self.cur.execute("SELECT ST_Extent(geom) FROM areas WHERE id = %s;" % str(area_id))
        except Exception, e:
            print e[0]
            print "ERROR: areas table not available"
            return -1
        else:
            bounds = self.cur.fetchall()
            if (len(bounds) != 0):
                return bounds[0][0]         
            else:
                print "ERROR: Requested area id not available"
                return -1
        
      
    
    def get_building_footprints(self):
		"""
		Returns a list of the building footprint polygons in WKT format
		"""     

		self.cur.execute("SELECT (ST_AsText(geom)) FROM buildings;")
		footprintsList = []
		
		footprints = self.cur.fetchall() 
		
		for row in footprints:
			footprintsList.append(row[0])
		
		return footprintsList


    def get_building_ids(self):
        """
        Returns a list of the building ids defined in the domain
        """     

        self.cur.execute("SELECT id FROM buildings;")
        buildings_ids = []
        
        buildings = self.cur.fetchall() 
        
        for b in buildings:
            buildings_ids.append(b[0])
        
        return buildings_ids


    def get_closest_elements(self, xy, epsgIN):
        """
        Given a list of points [x,y], return a list of ids corresponding to the element that 
        is closest to each point 
        
        Return the list of points with the corresponding elements ids and distance (in meters to the element centroid)
        return [(x,y,id,distance),...,]
        
        """ 
        #To RETURN:  the ids of closest node to each xy point
        elementIds = []
        distanceList = [] 
        for pt in xy:
            wkb = 'POINT(%s %s)' % (pt[0],pt[1])

            self.cur.execute("SELECT (ST_Distance(ST_Centroid(geom),ST_Transform(ST_GeomFromText('%s',%s),%s)), id) \
                              FROM elements \
                              WHERE ST_Contains(geom,ST_Transform(ST_GeomFromText('%s',%s),%s))" % 
                              (wkb,epsgIN,self.epsg,wkb,epsgIN,self.epsg))
            
                          
           # ST_AsText(ST_Transform(ST_Centroid(geom),%s))
            element = self.cur.fetchall()
            if(len(element) == 0):
                elementIds.append(0)
                distanceList.append(0)
                
            else:
                closestEl = element[0][0]
                i = closestEl.find(',')
                j = closestEl.find(')')
                id = closestEl[i+1:j]
                i = closestEl.find('(')
                j = closestEl.find(',')
                distance = closestEl[i+1:j]
                elementIds.append(int(id))
                distanceList.append(float(distance))
                
        return elementIds,distanceList
                

    def get_closest_nodes(self, xy,epsgIN):
        """
        Given a list of points [x,y], return a list of ids corresponding to the node that 
        is closest to each point 
        
        Return the list of points with the corresponding node ids and distance (in meters)
        return [(x,y,id,distance),...,]
        
        """
        #To RETURN:  the ids of closest node to each xy point
        nodeIds = [] 
        distanceList = []
        elevationList = [] 
        xy = self._convert_points_(pointsIN=xy,epsgIN=epsgIN,epsgOUT=self.epsg)

        
        i = 0
        for pt in xy:

            wkb = 'POINT(%s %s)' % (pt[0],pt[1])
            radius = 10
            x = float(pt[0])
            y = float(pt[1])
            
            ptsFound = False
            while(ptsFound == False):
                boundingBOX = 'POLYGON((%s %s, %s %s, %s %s, %s %s))' % (x-radius,y-radius, x+radius,y-radius, x+radius, y+radius, x-radius, y-radius)            
                self.cur.execute("SELECT \
                                  ST_Distance( geom,ST_GeomFromText('%s',%s) ), id, ST_X(geom),ST_Y(geom), xyz[3] \
                                  FROM nodes \
                                  WHERE ST_Contains(ST_GeomFromText('%s',%s),geom) \
                                  AND ST_DWithin(geom,ST_GeomFromText('%s',%s),%s) \
                                  ORDER BY 1" % 
                                  (wkb,self.epsg,boundingBOX,self.epsg,wkb,self.epsg, radius))
                                
                pts = self.cur.fetchall()                 
                if(len(pts) == 0):
                   radius = radius*2        #increase the search radius and try again
                   if(radius > 100000):     #Stop searching if there are no points nearby
                       break
                else:
                    ptsFound = True

            if(ptsFound == True):       #pts found
                x = pts[0][2]
                y = pts[0][3]
                z = pts[0][4]
                distance = pts[0][0]
                id = pts[0][1]
                                    
                nodeIds.append(int(id))
                distanceList.append(float(distance))
                elevationList.append(float(z))
                i+=1
                
            else:
                nodeIds.append(0)   #no point has been found within an acceptable distance
                distanceList.append(0)
                print "get_closest_nodes() : NO POINTS WITHIN required distance!!"
                print radius
                print pt
                sys.exit()
                

        return nodeIds,distanceList,elevationList


    def add_run(self,run_id, name,description,run_filename, grid_filename):
        
        """
        Add a Ricom run record to the run table
        
        
        """
        try: 
            self.cur.execute("INSERT INTO runs (pkey, id, name, description,run_filename,grid_filename) \
                            VALUES (%s, %s, '%s', '%s', '%s', '%s');" % (run_id, run_id, name, description,run_filename,grid_filename))
        except Exception, e:
            print "PostGIS ERROR: %sPlease change the RiCOM run identifiers." % e
            sys.exit() 
            return -1
        else:
            self.conn.commit()
            return 1
    
    def delete_run(self,run_id):
        
        """
        Remove a run and all its associated data from the database
        
        
        """
        #remove the run from the runs table
        print "GRIDPG: Deleting run_id = % from database." % run_id
        try: self.cur.execute("DELETE FROM runs WHERE id = %s;" % run_id)
        except Exception, e:
            print "PostGIS ERROR: %s" % e
            return -1
        else:
            #remove all the data from output data tables
            self.cur.execute("DELETE FROM output_nodes WHERE run_id = %s" % run_id)
            self.cur.execute("DELETE FROM output_elements WHERE run_id = %s" % run_id)
            self.cur.execute("DELETE FROM output_sides WHERE run_id = %s" % run_id)
            self.conn.commit()
            print "RUN with run_id = %s deleted from the PostGIS database" % run_id
            return 1
   
   
    def is_run_defined(self,run_id):
        
        """
        Check to see if a run is defined in the Database
        
        INPUT: run_id = id of the run
        """
        
        #check if the given Run_id is defined in the database 
        self.cur.execute("SELECT * FROM runs WHERE id = %s;" % (run_id))
        run = self.cur.fetchall() 
        if(len(run) != 1):
            print "Run with ID = %s not available in database." % run_id
            return False
        else:
            return True
      

    
  
    def get_elements_inside_buildings(self,area_id = 0):
        """
        Get all the elements corresponding to the buildings that are with the given area polygon.
        
        area_id - the id (in the postgis table) of the area that the buildings of interest are inside

        

        """
        
        self.cur.execute("DROP TABLE IF EXISTS drag_elements;")
        self.cur.execute("CREATE TABLE drag_elements (id int4);")    
        self.cur.execute("SELECT AddGeometryColumn('drag_elements', 'geom', %s, 'POLYGON', 2);" % self.epsg) 
        self.conn.commit()       

        elements = []
        elements_building = []
        buildings_dict = {}
        
        if area_id == 0:
            #search the entire domain
                
            print "Get_elements_inside_buildings: Invalid area_id = 0"


                            
            
        else:
            a = self.get_area_by_id(area_id)        #get the area as a text object
            if (a != ""):                           #check if a VALID area geometry has been found in areas table       
        
                #SELECT all the buildings that are inside the area polygon
                self.cur.execute("SELECT b.id, ST_AsText(b.geom) FROM buildings AS b \
                                    WHERE ST_Contains(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" % (a,self.epsg))
                buildings = self.cur.fetchall() 
 
                for b in buildings:
                    
                    elements_building = []

                    id = int(b[0])
                    p = self._offset_polygon2_(b[1],0.1)
                
                
                    #SELECT all the elements intersecting  the building footprint polygon (i.e. b offset by -0.1)
#                    self.cur.execute("SELECT e.id, e.code,geom FROM elements AS e \
#                                        WHERE ST_Contains(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;" \
#                                        % (p, self.epsg))

                     

                    self.cur.execute("SELECT e.id, e.code,e.geom,e.node_ids FROM elements AS e, buildings AS b \
                                        WHERE ST_DWithin(b.geom,e.geom,0.05) \
                                        AND ST_Contains(b.geom,ST_Centroid(e.geom)) \
                                        AND b.id = %s \
                                        ORDER BY e.id;" % (id)) 
                    
                    r_elements = self.cur.fetchall()
                                    
                    for el in r_elements:
                        elements.append([el[0],el[1]])
                        elements_building.append(el[0])
                        self.cur.execute("INSERT INTO drag_elements (id,geom) VALUES (%s,ST_Geometry('%s'));" % (el[0],el[2]))
                    
                    buildings_dict[id] = {'drag_elements': elements_building}
                        
                        
            else:
                print "ERROR: Invalid area id selected"
                return elements
        
        self.conn.commit()       

        return elements,buildings_dict


    def get_nodes_at_building_edge(self,building_id):
        """
        
        
        Given a building_id returns all the nodes that are on the edge of the building
        
        """



 

        self.cur.execute("SELECT ST_AsText(geom), ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)) FROM buildings WHERE id=%s;" % building_id)
        b = self.cur.fetchall()[0]
        
        polygon = b[0]

        radius = 1000
        x = float(b[1])
        y = float(b[2])
        
        boundingBOX = 'POLYGON((%s %s, %s %s, %s %s, %s %s))' % (x-radius,y-radius, x+radius,y-radius, x+radius, y+radius, x-radius, y-radius)    


        i = polygon.find('POLYGON((')
        j = polygon.find('))')
        pts = polygon[i+9:j]
        pts = pts.replace(',',' ').split()
        i = 0
        line_string = 'LINESTRING('                    
        while i < len(pts)-2:
            line_string = line_string + "%s %s," % (pts[i], pts[i+1])
            i +=2
        
        line_string = line_string + "%s %s)" % (pts[i], pts[i+1]) 
        
                  
        self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                                WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) AND \
                                ST_Distance(ST_GeomFromText('%s',%s),n.geom) < 0.1;" % (boundingBOX,self.epsg,line_string,self.epsg))

        r_nodes = self.cur.fetchall()

        nodes = []
        
        for n in r_nodes:
            nodes.append(n[0])
        
        return nodes

    def get_nodes_at_building_edges(self,area_id = 0):
        """
        Get all the elements corresponding to the buildings that are with the given area polygon.
        
        area_id - the id (in the postgis table) of the area that the buildings of interest are inside

        
        """
        
        self.cur.execute("DROP TABLE IF EXISTS buildings_no_edges;")
        self.cur.execute("CREATE TABLE buildings_no_edges (building_id int4);")    
        self.cur.execute("SELECT AddGeometryColumn('buildings_no_edges', 'geom', %s, 'POLYGON', 2);" % self.epsg) 
        
        outfile = open('NodesAtBuildings_code2_v2.ngh',"w")
        nodes_at_buildings = []
        nodes_all = []
        
        if area_id == 0:
            #search the entire domain
                
            #SELECT all the buildings that are inside the domain 
            self.cur.execute("SELECT id, ST_AsText(geom) FROM buildings;")
            buildings = self.cur.fetchall() 

            for b in buildings:
                id = int(b[0])
                p1 = self._offset_polygon2_(b[1],0.1)
                p2 = self._offset_polygon2_(b[1],-0.1)
            
                #SELECT all the elements intersecting  the building footprint polygon (i.e. b offset by -0.1)
                self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                                    WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) AND NOT ST_Contains(ST_GeomFromText('%s',%s),n.geom);" \
                                    % (p1, self.epsg,p2, self.epsg))
                r_nodes = self.cur.fetchall()
            
                for node in r_nodes:
                    nodes_at_buildings.append([node[0],node[1]])


                            
            
        else:
            a = self.get_area_by_id(area_id)        #get the area as a text object
            if (a != ""):                           #check if a VALID area geometry has been found in areas table       
        
                #SELECT all the buildings that are inside the area polygon
                self.cur.execute("SELECT b.id, ST_AsText(b.geom), geom FROM buildings AS b \
                                    WHERE ST_Contains(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" % (a,self.epsg))
                buildings = self.cur.fetchall() 
                print "Number of Buildings = %s" % len(buildings)
                
#                self.cur.execute("SELECT s.id, s.node_ids FROM sides AS s, buildings AS b \
#                    WHERE ST_Distance(b.geom,s.geom) < 0.2;")
#                
#                sides = self.cur.fetchall() 
#                
#                print "Number of Sides = %s" % len(sides)


                #b = buildings[0]
                for b in buildings:
                    #convert building polygon to a linestring
                    #-----------------------------------------------
                    #WHY? - PostGIS functions like ST_Distance - return the distance to the 
                    #Polygon as an area not as a line
                    i = b[1].find('POLYGON((')
                    j = b[1].find('))')
                    pts = b[1][i+9:j]
                    pts = pts.replace(',',' ').split()
                    i = 0
                    line_string = 'LINESTRING('                    
                    while i < len(pts)-2:
                        line_string = line_string + "%s %s," % (pts[i], pts[i+1])
                        i +=2
                    
                    line_string = line_string + "%s %s)" % (pts[i], pts[i+1])   
                    
                    #print b[1]
                    #print line_string                
                    
                    id = int(b[0])
                    #p1 = self._offset_polygon2_(b[1],0.5)
                    p2 = self._offset_polygon2_(b[1],-0.2)
                
                    #SELECT all the elements intersecting  the building footprint polygon (i.e. b offset by -0.1)
                    #self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                    #                    WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) AND NOT ST_Contains(ST_GeomFromText('%s',%s),n.geom);" \
                    #                    % (p1, self.epsg,p2, self.epsg))
                    '''
    
                    self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                                        WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) \
                                        AND ST_Distance(ST_GeomFromText('%s',%s),n.geom) < 0.2 ;" % (p1, self.epsg,p1, self.epsg))
                                        
                    '''
#                    
#                    self.cur.execute("SELECT n.id, n.code FROM nodes AS n, buildings AS b \
#                                        WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) AND \
#                                        ST_Distance(b.geom,n.geom) < 0.2 \
#                                        AND NOT ST_Contains(ST_GeomFromText('%s',%s),n.geom) \
#                                        AND b.id = %s;" % (a,self.epsg,p2,self.epsg, id))
    
                    self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                                            WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) AND \
                                            ST_Distance(ST_GeomFromText('%s',%s),n.geom) < 0.1;" % (a,self.epsg,line_string,self.epsg))

                    r_nodes = self.cur.fetchall()
                    
                    if (len(r_nodes) < 4):
                        
                        self.cur.execute("INSERT INTO buildings_no_edges (building_id,geom) VALUES (%s,ST_Geometry('%s'));" % (id,b[2]))

                        
                        print "ID = %s, L = %s *****************    " % (id, len(r_nodes))

                        
    
                    for n in r_nodes:
                        nodes_at_buildings.append([n[0],n[1]])
                            
            else:
                print "ERROR: Invalid area id selected"
                return nodes_at_buildings
            
        
        print "# of Nodes at Building edges = %s" % len(nodes_at_buildings)

        
        #Get all the nodes in the domain
        self.cur.execute("SELECT id, xyz, code, neighbours FROM nodes")
        r_nodes = self.cur.fetchall()
        
        for n in r_nodes:
            nodes_all.append([n[0],n[1][0],n[1][1],n[2],n[1][2],n[3]])
            
             
        for n in nodes_at_buildings:
            nodes_all[n[0]-1][3] = 2
        
        
        print "Domain Size of Nodes = %s" % len(nodes_all)
        
        print "Writing the new NGH file..."
        
        outfile.write("#NGH\n       0.00000       0.00000       0.00000       0.00000         0\n      %s\n            %s\n" % (len(nodes_all), len(nodes_all[0][5]) ))
        for n in nodes_all:
            if n[3] == 10 or n[3] == 100:
                n[3] = 1
            line = ""
            line = "%s %s %s %s %s    " % (str(n[0]).ljust(10),str(n[1]).ljust(15),str(n[2]).ljust(15),str(n[3]).ljust(4),str(n[4]).ljust(10))
            
            for nbr in n[5]:
                line = line + "%s" % (str(int(nbr)).ljust(8))
            
            line = line + "\n"
            
            outfile.write(line)
        
        outfile.close()        
#        for n in r_nodes:
#            nodes_all.append([n[0],n[1]])      
#            
        self.conn.commit()       

        return nodes_all,nodes_at_buildings



    def get_side_lengths(self):
        """
        Returns a list of side lengths of each of the sides in the domain
        """

        self.cur.execute("SELECT ST_Length(geom) FROM sides AS s;")
        r_sides = self.cur.fetchall()
        side_lengths = []
        for s in r_sides:
            side_lengths.append(s[0])
            
        return side_lengths
            
    def get_element_areas(self):
        """
        Returns the areas of each element in the domain
        """

        self.cur.execute("SELECT ST_Area(geom) FROM elements;")
        r_elements = self.cur.fetchall()
        element_areas = []
        for e in r_elements:
            element_areas.append(e[0])
            
        return element_areas 
        
        

    def get_elements_at_building_edges(self):
        '''
        
        
        '''
        
        
        self.cur.execute("SELECT b.id, ST_AsText(b.geom) FROM buildings AS b \
                            WHERE b.id = 473;")
        
        
        buildings = self.cur.fetchall()
        b = buildings[0]

              
        #convert building polygon to a linestring
        #-----------------------------------------------
        #WHY? - PostGIS functions like ST_Distance - return the distance to the 
        #Polygon as an area not as a line
        i = b[1].find('POLYGON((')
        j = b[1].find('))')
        pts = b[1][i+9:j]
        pts = pts.replace(',',' ').split()
        i = 0
        line_string = 'LINESTRING('                    
        while i < len(pts)-2:
            line_string = line_string + "%s %s," % (pts[i], pts[i+1])
            i +=2
        
        line_string = line_string + "%s %s)" % (pts[i], pts[i+1])  


        self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                            WHERE ST_Distance(ST_GeomFromText('%s',%s),n.geom) < 0.1;" % (line_string,self.epsg))
        r_nodes = self.cur.fetchall()
        
        nodes = []
        for n in r_nodes:
            nodes.append(n[0])

         
        self.cur.execute("SELECT e.id, e.node_ids,e.code FROM elements AS e \
                            WHERE ST_Distance(ST_GeomFromText('%s',%s),e.geom) < 0.1;" % (line_string,self.epsg))

        r_elements = self.cur.fetchall()
        elements = []
        for e in r_elements:
            elements.append([e[0], e[1]])
            
            
        elementEdge = []
        for e in elements:
            i = 0
            for n in e[1]:
                if n in nodes:
                    i+=1
            
            if i == 2:
                elementEdge.append(e[0])
            

        print "All close elements to 473: %s" % (len(elements))
        for e in elements:
            print "%s" % e[0]
            
        print "Edge elements to Building 473: %s \n %s" % (len(elementEdge),elementEdge)
        return elementEdge

  

    def get_elements_in_area(self,area_id = 0):
		"""
		Get all the elements inside a given area.
		
		area_id - the id (in the postgis table) of the area
		
		return: An array of the ids and codes corresponding to elements inside the area
		
		NOTE: if area_id = 0 return all elements within the area
		
		"""
		elements = []        
		if area_id == 0:				#select the whole domain
			self.cur.execute("SELECT id, code FROM elements")
			r_elements = self.cur.fetchall()
			for el in r_elements:
				elements.append([el[0],el[1]])
				
			
		else:
			a = self.get_area_by_id(area_id)        #get the area as a text object
			if (a != ""):                           #check if a VALID area geometry has been found in areas table       
			
				'''
				#get all the elements in the domain and their corresponding element codes
				self.cur.execute("SELECT id, code FROM elements")
				r_elements = self.cur.fetchall()
				for el in r_elements:
					elements.append([el[0],el[1]])
				'''
	
				#SELECT all the elements inside the given area 
				self.cur.execute("SELECT e.id, code FROM elements AS e WHERE ST_Contains(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;"  % (a,self.epsg) )
				r_elements = self.cur.fetchall()
				for el in r_elements:
					elements.append([el[0],el[1]])
					
			else:
				print "ERROR: Invalid area id selected"
				return elements
			
		return elements
    
    def init_building_output_tables(self):
        """
        
        
        
        """
        
        self.cur.execute("DROP TABLE IF EXISTS building_sides;")
        self.cur.execute("CREATE TABLE building_sides (id int4);")    
        self.cur.execute("SELECT AddGeometryColumn('building_sides', 'geom', %s, 'LINESTRING', 2);" % self.epsg) 

        self.cur.execute("DROP TABLE IF EXISTS building_elements;")
        self.cur.execute("CREATE TABLE building_elements (id int4);")    
        self.cur.execute("SELECT AddGeometryColumn('building_elements', 'geom', %s, 'POLYGON', 2);" % self.epsg) 
 
 
        self.cur.execute("DROP TABLE IF EXISTS building_nodes;")
        self.cur.execute("CREATE TABLE building_nodes (id int4);")    
        self.cur.execute("SELECT AddGeometryColumn('building_nodes', 'geom', %s, 'POINT', 2);" % self.epsg) 
 
        self.conn.commit()            
        

    def get_building_output_locations(self,area_id,type):
        """
        For each building polygon get id's of the nodes, elements and sides that are inside it
        
        area_id: id of the area where building output is requested

        type:        BUILDINGS_AS_HOLES
                    BUILDINGS_AS_POINTS
                    BUILDINGS_GRIDDED
        
        Returns:    a dictionary of sides,elements and nodes corresponding to each footprint
                    Depending on the type of grid, different criteria for extracting the output locations is used
        



        
        """                        
  
        validArea = False
        if area_id == 0:
            self.cur.execute("SELECT b.id, ST_AsText(b.geom), ST_Perimeter(b.geom),ST_X(ST_Centroid(b.geom)),ST_Y(ST_Centroid(b.geom)) FROM buildings AS b ORDER BY b.id;")
            buildings = self.cur.fetchall() 
            validArea = True

        else:   
            a = self.get_area_by_id(area_id)        #get the area as a text object
            if (a != ""):                           #VALID area geometry has been found in areas table         
                self.cur.execute("SELECT b.id, ST_AsText(b.geom), ST_Perimeter(b.geom),ST_X(ST_Centroid(b.geom)),ST_Y(ST_Centroid(b.geom)) FROM buildings AS b WHERE ST_Contains(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" % (a,self.epsg))
                buildings = self.cur.fetchall()
                validArea = True


        if validArea:
            
            nodesDict = {}
            eleDict = {}
            sidesDict = {}
            dictionary = {}

            print "Getting building output locations" 


            if type == "BUILDINGS_AS_POINTS":
                for b in buildings:                
                    id = int(b[0])
                    building_nodes,distances,elevations = self.get_closest_nodes([[float(b[3]),float(b[4])]], self.epsg) 
                    dictionary[id] = {'nodes': building_nodes, 'elements': [], 'sides': [], 'perimeter': float(b[2]), 'type':type}     #building is Gridded in mesh

            
            elif type == "BUILDINGS_AS_HOLES" or type == "BUILDINGS_GRIDDED": 
                for b in buildings:    
                    #-----------------------------------------------
                    #Convert building polygon to a LINESTRING
                    #-----------------------------------------------
                    #WHY? - PostGIS functions like ST_Distance - return the distance to the 
                    #Polygon as an area not as a line
                    #-----------------------------------------------
                    i = b[1].find('POLYGON((')
                    j = b[1].find('))')
                    pts = b[1][i+9:j]
                    pts = pts.replace(',',' ').split()
                    i = 0
                    line_string = 'LINESTRING('                    
                    while i < len(pts)-2:
                        line_string = line_string + "%s %s," % (pts[i], pts[i+1])
                        i +=2
                    
                    line_string = line_string + "%s %s)" % (pts[i], pts[i+1])    
       
                    id = int(b[0])
                    building_nodes = []
                    building_nodes_neighbours_dict = {}
                    
                    self.cur.execute("SELECT n.id, n.code, n.neighbours,n.geom,n.xyz FROM nodes AS n, buildings AS b \
                                                WHERE ST_DWithin(b.geom,n.geom,0.05) \
                                                AND b.id = %s;" % (id))
    
                    r_nodes = self.cur.fetchall()
                    for n in r_nodes:
                        building_nodes.append(n[0])                
                        building_nodes_neighbours_dict[n[0]] = n[2]
                        self.cur.execute("INSERT INTO building_nodes (id,geom) VALUES (%s,ST_Geometry('%s'));" % (n[0],n[3]))
                                                
                    building_sides = []
                    building_sides_db = []
                    
                    
                    
                    if type == "BUILDINGS_AS_HOLES":
                        #GET the sides that make up the building edge - NOT interested in the interior sides
    
                        self.cur.execute("SELECT s.id, s.node_ids, geom FROM sides AS s \
                                                WHERE ST_DWithin(ST_GeomFromText('%s',%s),s.geom,0.05);" % (line_string,self.epsg))
        
                        r_sides = self.cur.fetchall()
                        
                        #check to see if the side is actually part of the edge
                        for s in r_sides:
                            node_ids = s[1]
                            #only choose nodes which are part of the nodes list above (i.e. nodes on the buildings edge)
                            if node_ids[0] in building_nodes and node_ids[1] in building_nodes:  
                                n1_set = set(building_nodes_neighbours_dict[node_ids[0]])                        
                                n2_set = set(building_nodes_neighbours_dict[node_ids[1]])
                                
                                #check if the side spans across the corner of the building (i.e. not actually a building side)
                                matches = n1_set.intersection(n2_set)       
                                matches.discard(0)
                                if len(matches) == 1:
                                    building_sides.append(s[0])
                                    building_sides_db.append([s[0],s[2]])
        
                        for side in building_sides_db:                    
                            self.cur.execute("INSERT INTO building_sides (id,geom) VALUES (%s,ST_Geometry('%s'));" % (side[0],side[1]))
                   
                    else:
                        self.cur.execute("SELECT s.id, s.node_ids, geom FROM sides AS s \
                                            WHERE ST_DWithin(ST_GeomFromText('%s',%s),s.geom,0.05) \
                                            AND ST_DWithin(ST_GeomFromText('%s',%s),ST_Centroid(s.geom),0.05);" % (line_string,self.epsg,line_string,self.epsg))
                        r_sides = self.cur.fetchall()
                        for s in r_sides:
                            building_sides.append(s[0])
                            building_sides_db.append([s[0],s[2]])
                        
                    for side in building_sides_db:                    
                        self.cur.execute("INSERT INTO building_sides (id,geom) VALUES (%s,ST_Geometry('%s'));" % (side[0],side[1]))          
                    
                    building_elements = []
                    #IF Buildings are HOLES select the elements that are along the edges of the building
                    if type == "BUILDINGS_AS_HOLES":  
                        self.cur.execute("SELECT e.id, e.geom FROM elements AS e, buildings AS b \
                                WHERE ST_DWithin(b.geom,e.geom,0.1) \
                                AND b.id = %s;" % (id))
                    
                    #IF Buildings are Gridded select the elements that are inside te building
                    else:
                        self.cur.execute("SELECT e.id,e.geom FROM elements AS e, buildings AS b \
                                            WHERE ST_DWithin(b.geom,e.geom,0.1) \
                                            AND ST_Contains(b.geom,ST_Centroid(e.geom)) \
                                            AND b.id = %s \
                                            ORDER BY e.id;" % (id))             
    
    
                    r_elements = self.cur.fetchall()
                    
                    for e in r_elements:
                        building_elements.append(e[0])
                        self.cur.execute("INSERT INTO building_elements (id,geom) VALUES (%s,ST_Geometry('%s'));" % (e[0],e[1]))
                        
                    if( len(building_nodes) > 0 or len(building_sides) > 0 or len(building_elements) > 0):
                        
                        dictionary[id] = {'nodes': [], 'elements': [], 'sides': [], 'perimeter': float(b[2]), 'type':type}
                        if(building_nodes != []):
                            dictionary[id]['nodes'] = building_nodes           #add to dictionary                   
                        if(building_sides != []):                    
                            dictionary[id]['sides'] = building_sides           #add to dictionary 
                        if(building_elements != []):                    
                            dictionary[id]['elements'] = building_elements           #add to dictionary                     

            else:
                print "get_building_output_locations(): Invalid building type given."
                sys.exit()
                
            self.conn.commit()
            return dictionary
        
        else:
            print "ERROR: No footprints in requested area"
            return -1;




"""
SOME USEFUL FUNCTIONS

"""

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





"INSERT INTO nodes (id, xyz, code, geom, neighbours)"

'''  NOTES - Create Lines and AREAS tables -  TODO!!!

CREATE TABLE lines (pkey integer PRIMARY KEY,id int4);
SELECT AddGeometryColumn('lines', 'geom', 4326, 'LINESTRING', 2);
ALTER TABLE lines ADD COLUMN name varchar(50);
ALTER TABLE lines ADD COLUMN comment varchar(128);
CREATE TABLE areas (pkey integer PRIMARY KEY,id int4);
SELECT AddGeometryColumn('areas', 'geom', 4326, 'POLYGON', 2);
ALTER TABLE areas ADD COLUMN name varchar(50);
ALTER TABLE areas ADD COLUMN comment varchar(128);

CREATE RULE force_id AS ON INSERT TO lines DO ALSO UPDATE lines SET id = NEW.pkey WHERE pkey = NEW.pkey;
CREATE RULE force_id AS ON INSERT TO areas DO ALSO UPDATE areas SET id = NEW.pkey WHERE pkey = NEW.pkey;

'''


#TODO
#    Add premod GRID
#    Determine output required at nodes,elements and sides:
#        Use the footprint polygons to determine (done by spatial searches in the postgis database):
#            The sides in the grid that make up a footprint
#            The elements in the grid that make up a footprint
#            The nodes in the grid that make up a footprint
#        Write the output locations to the PostGIS database
#        Write the output locations to the Netcdf grid file

#        