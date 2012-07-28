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
            
            print self.grid_dir + self.grid_filename
            self.grid_nc = GridNC(grid_filename = self.grid_dir + self.grid_filename, epsg = self.epsg)
            
            
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


    def get_building_output_locations(self,area_id):
        """
        Get id's of the nodes, elements and sides that are inside the buildings polygons
        
        Input:
            area_id - id of the area polygon that you want to get the info inside
        
        
        """ 
        dictionary = self.grid_pg.get_building_output_locations(area_id)
        return dictionary
        
    def get_element_output_locations(self, xy, epsgIN):
        """
        Given a list of points [x,y], return a list of ids corresponding to the element that 
        is closest to each point 

        
        """ 
        elementIds, distanceList = self.grid_pg.get_element_output_locations(xy,epsgIN)
        return elementIds, distanceList

    def get_node_output_locations(self, xy,epsgIN):
        """
        Given a list of points [x,y], return a list of ids corresponding to the node that 
        is closest to each point 
        
        """         
        nodeIds, distanceList = self.grid_pg.get_node_output_locations(xy,epsgIN)
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
	
	
	
	
    def get_elements_in_buildings(self,area_id = 0):
		"""
		Get all the elements inside a given area and inside the building polygons.  If area_id = 0, return all elements
		
		area_id - the id (in the postgis table) of the area
		
		
		"""
		return self.grid_pg.get_elements_in_buildings(area_id)
	
	
    def get_elements_in_buildings2(self,area_id = 0):
        """
        Get all the elements inside a given area and inside the building polygons.  If area_id = 0, return all elements
        
        area_id - the id (in the postgis table) of the area
        
        
        """
        return self.grid_pg.get_elements_in_buildings2(area_id)	
	
    
    def get_nodes_at_building_edges(self,area_id = 0):
        """
        Get all the elements inside a given area and inside the building polygons.  If area_id = 0, return all elements
        
        area_id - the id (in the postgis table) of the area
        
        
        """
        return self.grid_pg.get_nodes_at_building_edges(area_id)    		
	
	
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
        
        
    def add_max_values(self, footprints_max):
        """
        Adds the max values (Flood Depth and flow Speed) at each footprint to the postgis database (buildings table)
        
        IN:  footprints_max - dictionary of fd_max and speed_max values for each footprint
        
        """     
        
        
        #ADD new max value columns to the buildings table
#        try: self.cur.execute("ALTER TABLE buildings ADD COLUMN fd_max float DEFAULT 0;")
#        except:
#            print "add_max_values: Can't add fd_max column to the buildings table"
#            return 0
#        
#        try: self.cur.execute("ALTER TABLE buildings ADD COLUMN speed_max float DEFAULT 0;")
#        except:
#            print "add_max_values: Can't add speed_max column to the buildings table"
#            return 0
        
        for footprint_id, data in footprints_max.iteritems():
            self.cur.execute("UPDATE buildings SET fd_max = %s, speed_max = %s WHERE id = %s;" % (data['fd_max'],data['speed_max'],footprint_id))
            
        self.conn.commit()
        
    def _add_table_buildings_(self):	
		try: self.cur.execute("CREATE TABLE buildings (pkey integer PRIMARY KEY,id int4);")     
		except: print "Table \"buildings\" already exists."
		else: 
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
								ALTER TABLE buildings ADD COLUMN comment varchar(128);")        	
				
			#CREATE RULE for when new records are added to the buildings table
			#This rules makes the ID of new building record equal to PKEY
			self.cur.execute("CREATE RULE force_id AS ON INSERT TO buildings \
								DO ALSO \
								UPDATE buildings SET id = NEW.pkey WHERE pkey = NEW.pkey;")
			self.conn.commit()	
		
    def add_footprints(self, filename, layername = "", boundaryID = 1):
        """
        This function imports the building footprints for the project
        
        IN: filename = file name containing footprints POLYGONS for the study area
            Input file is generally a KML or ESRI Shapefile containing building
            footprints as POLYGONS or MULTIPOLYGONS
            layername = name of layer containing footprints
            If no layername is given function assumes that the first layer (index=0)
            contains the footprints data
            boundaryID - the id of the boundary polygon inside which the footprints are to be imported (i.e. the study area)
            
        OUT: function returns # of geometries imported if successful or 0 if unsuccessful
        
        """
        n = 0
        self.cur.execute("DROP TABLE IF EXISTS buildings;")
        self._add_table_buildings_()					#Create new buildings table
        #NOTE the PRIMARY KEY is pkey is required for Geoserver to write to the table via WFS_T
        #self.cur.execute("CREATE TABLE buildings (pkey integer PRIMARY KEY, id int4);")
        #self.cur.execute("SELECT AddGeometryColumn('buildings', 'geom', %s, 'POLYGON', 2);" % self.epsg)      
        self.cur.execute("SELECT * FROM buildings LIMIT 1;" )
        if len(self.cur.fetchall()) == 0:                   #If table is empty then fill up with footprints
            datasource_in = ogr.Open(filename)    
            if datasource_in is None:
                print "Could not open footprints file.\n"
                return n        
            #print datasource_in.GetLayerCount()
            #layer_in = datasource_in.GetLayerByName("footprints")
            
            if datasource_in.GetLayer(0) < 1:
                print "No data in the input file. \n"
                return n
            
            layer_in = datasource_in.GetLayer(0)
            layer_in.ResetReading()
            
            srs_in = layer_in.GetSpatialRef()       # the footprints spatial reference
            srs = osr.SpatialReference()            # the spatial reference of the PostGIS database (i.e. self.dbname)
            srs.ImportFromEPSG(self.epsg)
                 
            transform = 0
            if srs_in != srs:
                transform = 1
    
                       
            m = 0
            i = 0
            sql = ""
            sqlPTVA = ""
    
            for feature in layer_in:
                geom = feature.GetGeometryRef()
                #geom.wkbFlatten()               #flatten geometry if it is defined with the 25D flag
                if geom is not None and (geom.GetGeometryType() ==  ogr.wkbPolygon or geom.GetGeometryType() ==  ogr.wkbMultiPolygon):
                    
                    if (geom.GetGeometryType() ==  ogr.wkbPolygon):
                        valid = geom.IsValid()
                        if valid:
                            if transform:
                                geom.TransformTo(srs) 
                            
                            wkt = geom.ExportToWkt()
                            self.cur.execute("SELECT b.id FROM boundary as b WHERE (ST_Contains(b.geom,ST_GeomFromText('%s',%s))) AND b.id=%s" % (wkt, self.epsg, boundaryID))        
                            #if the current footprint is contained within the study boundary polygon then add to SQL command list
                            if len(self.cur.fetchall()) != 0:
                                sql = sql + "INSERT INTO buildings (pkey,id,geom) VALUES (%s,%s,ST_GeomFromText('%s',%s));\n" % (n+1,n+1, wkt, self.epsg)                             
                                sqlPTVA = sqlPTVA + "INSERT INTO ptva (building_id) VALUES (%s);\n" % (n+1)
                                n += 1
                        else:
							i += 1
                    elif (geom.GetGeometryType() ==  ogr.wkbMultiPolygon):
                        print "Warning: footprints layer contains MULTIPOLYGONS."
                        m += 1
    
                feature.Destroy()
            
            print "Number of building POLYGON footprints imported = %s" % n
            print "Number of building MULTIPOLYGON footprints NOT imported = %s" % m
            print "Number of invalid geometries = %s" % i
            
            if (sql != "") and (sqlPTVA != ""):
                self.cur.execute(sql)
                self.cur.execute(sqlPTVA)
                print "Creating spatial index on buildings...."
                self.cur.execute("CREATE INDEX buildings_idx ON buildings USING GIST (geom);")

                
                #ADD extra PTVA columns to the 'buildings' table
                '''
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
                                ALTER TABLE buildings ADD COLUMN comment varchar(128);")
                '''  


				#CREATE RULE for when new records are added to the buildings table
				#This rules makes the ID of new building record equal to PKEY
				#self.cur.execute("CREATE RULE force_id AS ON INSERT TO buildings DO ALSO UPDATE buildings SET id = NEW.pkey WHERE pkey = NEW.pkey;")
                self.conn.commit()

    
            datasource_in.Destroy()
        

        
            return n
        
        else: 
            print "Footprints table already full."
            return n

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
 
 
            
    def _offset_polygon_(self, offset=0.1, table = 'buildings',column = 'geom'):
        """
        Offset the building footprint polygons by passed amount offset (in meters)
        PRIVATE function to the PGTsunami CLASS
        
        """
        #Assumes that the input polygon is closed (i.e. the endpoint is the same as the start)
        #Assumes that input polygon is Clockwise defines (is this the standard??)
        #polyIn is a list of x,y points
      
      
        import math
        self.cur.execute("SELECT id, ST_AsText(%s) FROM %s;" %(column, table))
        polySQL = self.cur.fetchall()
        polyList = []
        id = []
        
        
        for poly in polySQL:    #extract the polygon vertices from the POSTGIS polygon    
            id.append(poly[0])
            poly = poly[1]
            polygon = []
            i = poly.find('POLYGON((')
            j = poly.find('))')
            points = poly[i+9:j]
            points = points.replace(',',' ').split()
    
            i = 0
            while i < len(points):
                polygon.append([float(points[i]), float(points[i+1])])
                i += 2
           
            polyList.append(polygon)
            
        polyOffsetList = []        
        for poly in polyList:
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
            
            polyOffsetList.append(polyOffset)
            
        #convert to wkt formated polygon list
        polyOffsetListWKT = []
        for poly in polyOffsetList:
            wkt = "POLYGON(("
            i = 0
            nPts = len(poly)
            while i < nPts - 1:
                ptStr = "%s %s," % (poly[i][0], poly[i][1])
                wkt = wkt + ptStr
                i += 1    
            ptStr = "%s %s))" % (poly[i][0], poly[i][1])
            wkt = wkt + ptStr
            polyOffsetListWKT.append(wkt)
        
        #add offset building footprints to a new table
        self.cur.execute("DROP TABLE IF EXISTS offset1;")
        self.cur.execute("CREATE TABLE offset1 (id int4);")
        self.cur.execute("SELECT AddGeometryColumn('offset1', 'geom', %s, 'POLYGON', 2);" % self.epsg)
        
        n = 0
        for poly in polyOffsetListWKT:
            self.cur.execute("INSERT INTO offset1 (id, geom) VALUES (%s, ST_GeomFromText('%s',%s));" % (polySQL[n][0], poly, self.epsg))            
            n += 1
            
        self.conn.commit()

        
        return id, polyOffsetListWKT
            
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
        
    def add_boundary(self,ply_in="../output/SQL/KBBuildingBoundingPoly2.ply",epsgIN = 28355):
        """
        """  
        ply_file = open(ply_in, "r").readlines() 
        boundingPoly = []
        for row in ply_file:
            row =  row.split()
            boundingPoly.append(row)
            
            print "row = %s" % row
            
     	print "self.epsg = %s" % self.epsg   
        
        boundingPoly = self._convert_points_(boundingPoly, epsgIN,epsgOUT=self.epsg) 
           
        self.cur.execute("DROP TABLE IF EXISTS boundary;")
        self.cur.execute("CREATE TABLE boundary (id int4, description varchar(64));")
        self.cur.execute("SELECT AddGeometryColumn('boundary', 'geom', %s, 'POLYGON', 2);" % self.epsg)
        
        sql = "ST_GeomFromText('POLYGON(("
        numPts = len(boundingPoly)
        i = 0
        while i < (numPts - 1):
            ptText = "%s %s, " % (str(boundingPoly[i][0]),str(boundingPoly[i][1]))
            sql = sql + ptText
            i += 1
            
        ptText = "%s %s))',%s));" % (str(boundingPoly[i][0]),str(boundingPoly[i][1]),self.epsg)
        sql = sql + ptText
    
        ptText = "INSERT INTO boundary (ID,description,geom) VALUES (1,'Ishinomaki Study Area',"
        sql = ptText + sql
        self.cur.execute(sql)  
        self.conn.commit()

        self.cur.execute("SELECT ST_AsText(geom) FROM boundary;")  
        boundary = self.cur.fetchall() 
        self.geom_to_kml(boundary[0],kml_file="boundary.kml", epsgIN=self.epsg)

        return
  
    
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

    def get_element_output_locations(self, xy, epsgIN):
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
                

    def get_node_output_locations(self, xy,epsgIN):
        """
        Given a list of points [x,y], return a list of ids corresponding to the node that 
        is closest to each point 
        
        Return the list of points with the corresponding node ids and distance (in meters)
        return [(x,y,id,distance),...,]
        
        """
        #To RETURN:  the ids of closest node to each xy point
        nodeIds = [] 
        distanceList = [] 
        for pt in xy:

            wkb = 'POINT(%s %s)' % (pt[0],pt[1])
            radius = 10
            
            ptsFound = False
            while(ptsFound == False):
                self.cur.execute("SELECT ( \
                                  ST_Distance(geom,ST_Transform(ST_GeomFromText('%s',%s),%s)), id) \
                                  FROM nodes \
                                  WHERE ST_DWithin(geom,ST_Transform(ST_GeomFromText('%s',%s),%s),%s) \
                                  ORDER BY 1" % 
                                  (wkb,epsgIN,self.epsg,wkb,epsgIN,self.epsg, radius))
                                
                pts = self.cur.fetchall() 
                if(len(pts) == 0):
                   radius = radius*2        #increase the search radius and try again
                   if(radius > 100000):     #Stop searching if there are no points nearby
                       break
                else:
                    ptsFound = True

            if(ptsFound == True):       #pts found
                closestPt = pts[0][0]
                i = closestPt.find(',')
                j = closestPt.find(')')
                id = closestPt[i+1:j]
                i = closestPt.find('(')
                j = closestPt.find(',')
                distance = closestPt[i+1:j]

                nodeIds.append(int(id))
                distanceList.append(float(distance))
            else:
                nodeIds.append(0)   #no point has been found within an acceptable distance
                distanceList.append(0)

        return nodeIds,distanceList


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
      

    def get_output_at_building(self,run_id,building_id):
        """
        
        
        
        """
        nodes_dictionary = {}
        if(self.is_run_defined(run_id)):
            #Get node outputs
            self.cur.execute("SELECT nodes FROM buildings WHERE id = %s;" % (building_id))
            nodes_dictionary = {'node_ids': [], 'time': [], 'eta': {}, 'uv': {}}
            nodes = self.cur.fetchall() 
            nodes_dictionary['node_ids'] = nodes[0][0]

            for id in nodes[0][0]:
                self.cur.execute("SELECT eta, uv, time FROM output_nodes WHERE run_id = %s AND node_id = %s;" % (run_id, id))
                output = self.cur.fetchall() 
                if(len(output) != 0):
                    nodes_dictionary['eta'][id] = output[0][0]
                    nodes_dictionary['uv'][id] = output[0][1]

            if(len(output) != 0): nodes_dictionary['time'] = output[0][2]

            #Get element outputs
            self.cur.execute("SELECT elements FROM buildings WHERE id = %s;" % (building_id))
            elements_dictionary = {'element_ids': [], 'time': [], 'eta': {}, 'uv': {}}
            elements = self.cur.fetchall() 
            elements_dictionary['element_ids'] = elements[0][0]

            for id in elements[0][0]:
                self.cur.execute("SELECT eta, uv, time FROM output_elements WHERE run_id = %s AND element_id = %s;" % (run_id, id))
                output = self.cur.fetchall()
                if(len(output) != 0):
                    elements_dictionary['eta'][id] = output[0][0]
                    elements_dictionary['uv'][id] = output[0][1]
            
            if(len(output) != 0): elements_dictionary['time'] = output[0][2]

            #Get side outputs
            self.cur.execute("SELECT sides FROM buildings WHERE id = %s;" % (building_id))
            sides_dictionary = {'side_ids': [], 'time': [], 'eta': {}, 'uv': {}}
            sides = self.cur.fetchall() 
            sides_dictionary['side_ids'] = sides[0][0]

            for id in sides[0][0]:
                self.cur.execute("SELECT eta, uv, time FROM output_sides WHERE run_id = %s AND side_id = %s;" % (run_id, id))
                output = self.cur.fetchall()
                if(len(output) != 0):
                    sides_dictionary['eta'][id] = output[0][0]
                    sides_dictionary['uv'][id] = output[0][1]

            if(len(output) != 0): sides_dictionary['time'] = output[0][2]


        return nodes_dictionary,elements_dictionary,sides_dictionary

    
    def add_output_at_node(self, run_id, node_id, eta, uv, time):
        """
        Add output at an node
        
        """        
        
        #check if the given Node_id is defined in the grid
        self.cur.execute("SELECT id FROM nodes WHERE id = %s;" % (node_id))
        node = self.cur.fetchall()    
        if(len(node) != 1):
            print "Node with ID = %s not available in database." % node_id
            return -1
        
        self.cur.execute("SELECT node_id FROM output_nodes WHERE run_id = %s AND node_id = %s LIMIT 1;" % (run_id, node_id))
        r = self.cur.fetchall()    
        if(len(r) == 1):
            print "Data from run_id = %s already in Database.  Please delete run and try again;" % run_id
            sys.exit()


        if eta == []: eta = [float(0.0)]
        
        if uv == []: uv = [[float(0.0),float(0.0)]]    
        
        self.cur.execute("INSERT INTO output_nodes (run_id, node_id, eta,uv, time) \
                                 VALUES (%s, %s, ARRAY%s, ARRAY%s,ARRAY%s);" 
                                 % (run_id, node_id, eta,uv, time))             
        
        #self.conn.commit()




        
    
    def add_output_at_element(self,run_id, element_id, eta,uv, time):
        """
        Add output at an element
        
        
        """
        #check if the given Node_id is defined in the grid
        self.cur.execute("SELECT id FROM elements WHERE id = %s;" % (element_id))
        element = self.cur.fetchall()    
        if(len(element) != 1):
            print "Element with ID = %s not available in database." % element_id
            return -1
        
        self.cur.execute("SELECT element_id FROM output_elements WHERE run_id = %s AND element_id = %s LIMIT 1;" % (run_id, element_id))
        r = self.cur.fetchall()    
        if(len(r) == 1):
            print "Data from run_id = %s already in Database.  Please delete run and try again;" % run_id
            sys.exit()
    
        
        if eta == []:
            eta = [float(0.0)]

        if uv == []:
            uv = [[float(0.0),float(0.0)]]    

        self.cur.execute("INSERT INTO output_elements (run_id, element_id, eta, uv, time) \
                         VALUES (%s, %s, ARRAY%s, ARRAY%s,ARRAY%s);" 
                         % (run_id, element_id, eta, uv,time))             
        
#        self.conn.commit()




 
    
    def add_output_at_side(self,run_id, side_id, eta,uv, time):
        """
        Add output at a side
        
        
        """
        #check if the given side_id is defined in the grid
        self.cur.execute("SELECT id FROM sides WHERE id = %s;" % (side_id))
        side = self.cur.fetchall()    
        if(len(side) != 1):
            print "Node with ID = %s not available in database." % side_id
            return -1
        
        self.cur.execute("SELECT side_id FROM output_sides WHERE run_id = %s AND side_id = %s LIMIT 1;" % (run_id, side_id))
        r = self.cur.fetchall()    
        if(len(r) == 1):
            print "Data from run_id = %s already in Database.  Please delete run and try again;" % run_id
            sys.exit()

        
        if eta == []:
            eta = [float(0.0)]

        if uv == []:
            uv = [[float(0.0),float(0.0)]]    
        
        self.cur.execute("INSERT INTO output_sides (run_id, side_id, eta, uv, time) \
                                 VALUES (%s, %s, ARRAY%s, ARRAY%s, ARRAY%s);" 
                                 % (run_id, side_id, eta,uv, time))             
        
#        self.conn.commit()        
        
        
    def add_builiding_topology(self):
        """
        INSERT into the "buildings" table all the nodes, elements and sides in the grid that correspond to each footprints
        
        """
        
        self.cur.execute("SELECT b.id, ST_AsText(b.geom) FROM buildings AS b;")
        footprints = self.cur.fetchall() 
        
        
        n = 0
        for f in footprints:
                   
            id = int(f[0])
            p = self._offset_polygon2_(f[1])
            #iterate through the building footprint polygon list
                
            #Get all the nodes contained in the building footprint polygon (offset)
            self.cur.execute("SELECT n.id FROM nodes AS n WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) ORDER BY n.id;" % (p, self.epsg ))
            
            r_nodes = self.cur.fetchall()
            
            #Get all the elements intersecting  the building footprint polygon (offset)
            self.cur.execute("SELECT e.id FROM elements AS e \
                    WHERE ST_Intersects(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;" \
                    % (p, self.epsg))
            r_elements = self.cur.fetchall()
            
            #Get all the sides contained in the building footprint polygon (offset)
            self.cur.execute("SELECT s.id FROM sides AS s\
                    WHERE ST_Contains(ST_GeomFromText('%s',%s),s.geom) ORDER BY s.id;" \
                    % (p, self.epsg))
            
            r_sides = self.cur.fetchall()
            
            nodes = []
            elements = []
            sides = []

            
            for row in r_nodes: nodes.append(row[0])
            for row in r_elements: elements.append(row[0])
            for row in r_sides: sides.append(row[0])
            
            self.cur.execute("UPDATE buildings SET nodes = ARRAY%s WHERE id = %s;" % (nodes,id))
            self.cur.execute("UPDATE buildings SET elements = ARRAY%s WHERE id = %s;" % (elements,id))
            self.cur.execute("UPDATE buildings SET sides = ARRAY%s WHERE id = %s;" % (elements,id))
            
            self.conn.commit()

            
            n += 1
        return




    def get_elements_in_buildings(self,area_id = 0):
        """
        Get all the elements corresponding to the buildings that are with the given area polygon.
        
        area_id - the id (in the postgis table) of the area that the buildings of interest are inside

        
        """
        
        elements = []
        
        if area_id == 0:
        	#search the entire domain
			#SELECT all the elements intersecting  the building footprint polygon (i.e.b offset by -0.1)
			self.cur.execute("SELECT e.id, e.code FROM elements AS e, buildings AS b WHERE ST_Contains(b.geom,e.geom) ORDER BY e.id;" )
			r_elements = self.cur.fetchall()
			
			for el in r_elements:
				elements.append([el[0],el[1]])
				        	
        	
        else:
			a = self.get_area_by_id(area_id)        #get the area as a text object
			if (a != ""):                           #check if a VALID area geometry has been found in areas table       
	
				#SELECT all the elements intersecting  the building footprint polygon (i.e.b offset by -0.1)
				self.cur.execute("SELECT e.id, e.code FROM elements AS e, buildings AS b WHERE ST_Contains(b.geom,e.geom) AND ST_Contains(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;"  % (a,self.epsg) )
				r_elements = self.cur.fetchall()
				
				for el in r_elements:
					elements.append([el[0],el[1]])

				#SELECT all the elements intersecting  the building footprint polygon (i.e.b offset by -0.1)
				self.cur.execute("SELECT e.id, e.code FROM elements AS e, buildings AS b WHERE ST_Contains(b.geom,e.geom) AND ST_Intersects(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;"  % (a,self.epsg) )
				r_elements = self.cur.fetchall()
				
				for el in r_elements:
					elements.append([el[0],el[1]])

				'''
				self.cur.execute("SELECT e.id, e.code FROM elements AS e, buildings AS b WHERE ST_Intersects(b.geom,e.geom) ORDER BY e.id;")
				r_elements = self.cur.fetchall()
				'''
					
			else:
				print "ERROR: Invalid area id selected"
				return elements
        	
        return elements


    def get_elements_in_buildings2(self,area_id = 0):
        """
        Get all the elements corresponding to the buildings that are with the given area polygon.
        
        area_id - the id (in the postgis table) of the area that the buildings of interest are inside

        
        """
        
        elements = []
        
        if area_id == 0:
            #search the entire domain
                
            #SELECT all the buildings that are inside the domain 
            self.cur.execute("SELECT id, ST_AsText(geom) FROM buildings;")
            buildings = self.cur.fetchall() 

            for b in buildings:
                id = int(b[0])
                p = self._offset_polygon2_(b[1],0.1)
            
                #SELECT all the elements intersecting  the building footprint polygon (i.e. b offset by -0.1)
                self.cur.execute("SELECT e.id, e.code FROM elements AS e \
                                    WHERE ST_Contains(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;" \
                                    % (p, self.epsg))
                r_elements = self.cur.fetchall()
            
                for el in r_elements:
                    elements.append([el[0],el[1]])


                            
            
        else:
            a = self.get_area_by_id(area_id)        #get the area as a text object
            if (a != ""):                           #check if a VALID area geometry has been found in areas table       
        
                #SELECT all the buildings that are inside the area polygon
                self.cur.execute("SELECT b.id, ST_AsText(b.geom) FROM buildings AS b \
                                    WHERE ST_Contains(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" % (a,self.epsg))
                buildings = self.cur.fetchall() 
 
                for b in buildings:
                    id = int(b[0])
                    p = self._offset_polygon2_(b[1],0.1)
                
                    #SELECT all the elements intersecting  the building footprint polygon (i.e. b offset by -0.1)
                    self.cur.execute("SELECT e.id, e.code FROM elements AS e \
                                        WHERE ST_Contains(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;" \
                                        % (p, self.epsg))
                    r_elements = self.cur.fetchall()
                
                    for el in r_elements:
                        elements.append([el[0],el[1]])
                        
            else:
                print "ERROR: Invalid area id selected"
                return elements
            
        return elements



    def get_nodes_at_building_edges(self,area_id = 0):
        """
        Get all the elements corresponding to the buildings that are with the given area polygon.
        
        area_id - the id (in the postgis table) of the area that the buildings of interest are inside

        
        """
        
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
                self.cur.execute("SELECT b.id, ST_AsText(b.geom) FROM buildings AS b \
                                    WHERE ST_Contains(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" % (a,self.epsg))
                buildings = self.cur.fetchall() 
                print "Number of Buildings = %s" % len(buildings)
                
#                self.cur.execute("SELECT s.id, s.node_ids FROM sides AS s, buildings AS b \
#                    WHERE ST_Distance(b.geom,s.geom) < 0.2;")
#                
#                sides = self.cur.fetchall() 
#                
#                print "Number of Sides = %s" % len(sides)


                
                for b in buildings:
                    id = int(b[0])
                    p1 = self._offset_polygon2_(b[1],0.3)
                    p2 = self._offset_polygon2_(b[1],-0.3)
                
                    #SELECT all the elements intersecting  the building footprint polygon (i.e. b offset by -0.1)
                    self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                                        WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) AND NOT ST_Contains(ST_GeomFromText('%s',%s),n.geom);" \
                                        % (p1, self.epsg,p2, self.epsg))
                    '''

                    self.cur.execute("SELECT n.id, n.code FROM nodes AS n \
                                        WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) \
                                        OR ST_Distance(ST_GeomFromText('%s',%s),n.geom) < 0.2;" % (p1, self.epsg,p1, self.epsg))
                                        
                    '''
                    
#                    self.cur.execute("SELECT n.id, n.code FROM nodes AS n, buildings AS b \
#                                        WHERE ST_Distance(b.geom,n.geom) < 0.2 AND b.id = %s;" % (id))

                    r_nodes = self.cur.fetchall()
                    
                    
                    #if len(r_nodes) == 0:
                    
                    if (len(r_nodes) < 4):
                        print "ID = %s, L = %s *****************    " % (id, len(r_nodes))
                    else:
                        print "ID = %s, L = %s" % (id, len(r_nodes))
                        

                    for n in r_nodes:
                        nodes_at_buildings.append([n[0],n[1]])
                        
            else:
                print "ERROR: Invalid area id selected"
                return nodes_at_buildings
            
            
        #Get all the nodes in the domain
        self.cur.execute("SELECT id, code FROM nodes")
        r_nodes = self.cur.fetchall()
        for n in r_nodes:
            nodes_all.append([n[0],n[1]])      
            
        return nodes_all,nodes_at_buildings

  

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
    
    def get_building_output_locations(self,area_id):
        """
        For each building polygon get id's of the nodes, elements and sides that are inside it
        
        boundary_WKB: bounding polygon of the area where output is requested
        
        Returns: a dictionary of sides,elements and nodes corresponding to each footprint
        
        """                     


    
        
        a = self.get_area_by_id(area_id)        #get the area as a text object
        if (a != ""):                           #VALID area geometry has been found in areas table       
        
            self.cur.execute("SELECT b.id, ST_AsText(b.geom) FROM buildings AS b WHERE ST_Contains(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" % (a,self.epsg))
            footprints = self.cur.fetchall() 

            
            nodesDict = {}
            eleDict = {}
            sidesDict = {}
            dictionary = {}
            
            n = 0

            
            for f in footprints:
                
                
                id = int(f[0])
                p = self._offset_polygon2_(f[1])

                #iterate through the building footprint polygon list
                            
                #Get all the nodes contained in the building footprint polygon (offset)
                self.cur.execute("SELECT n.id, ST_AsText(n.geom) FROM nodes AS n \
                                WHERE ST_Contains(ST_GeomFromText('%s',%s),n.geom) ORDER BY n.id;"
                                % (p, self.epsg ))
                
                r_nodes = self.cur.fetchall()
    
                #Get all the elements intersecting  the building footprint polygon (offset)
                self.cur.execute("SELECT e.id, ST_AsText(e.geom) FROM elements AS e \
                                WHERE ST_Intersects(ST_GeomFromText('%s',%s),e.geom) ORDER BY e.id;" \
                                % (p, self.epsg))
                r_elements = self.cur.fetchall()

                #Get all the sides contained in the building footprint polygon (offset)
                self.cur.execute("SELECT s.id, ST_AsText(s.geom) FROM sides AS s\
                                WHERE ST_Contains(ST_GeomFromText('%s',%s),s.geom) ORDER BY s.id;" \
                                % (p, self.epsg))

                r_sides = self.cur.fetchall()
            

                if( len(r_nodes) > 0 or len(r_elements) > 0 or len(r_sides) > 0 ):
                    
                    dictionary[id] = {'nodes': [], 'elements': [], 'sides': []}

                    nodes = []
                    for row in r_nodes: nodes.append(row[0])
                    if r_nodes != []: dictionary[id]['nodes'] = nodes           #add to dictionary                   

                    elements = []
                    for row in r_elements: elements.append(row[0])
                    if r_elements != []: dictionary[id]['elements'] = elements     #add to dictionary
                   
                    sides = []
                    for row in r_sides: sides.append(row[0])
                    if r_sides != []: dictionary[id]['sides'] = sides           #add to dictionary   
                
                n += 1
                #self.conn.commit()
            
            return dictionary
            
        else:
            print "ERROR: No footprints in requested area"
            return -1;


    def add_shapefile(self, filename, tablename = "vectorlayer", boundaryID = 1):
        """
        This function imports a shapefile into a the postgis database
        
        IN: filename = shapefile filename - assume shapefile consists of LINESTRING only
            tablename = name of the table to be created in the database
            boundaryID - the id of the boundary polygon inside which the shapefile features are to be imported (i.e. the study area)
            
        OUT: function returns # of geometries imported if successful or 0 if unsuccessful
        
        """
        
        
        
        n = 0
        self.cur.execute("DROP TABLE IF EXISTS roads;")
        self.cur.execute("CREATE TABLE %s (id int4);" % tablename)
        self.cur.execute("SELECT AddGeometryColumn('%s', 'geom', %s, 'LINESTRING', 2);" % (tablename, self.epsg))      
        self.cur.execute("SELECT * FROM %s LIMIT 1;" % tablename )
        if len(self.cur.fetchall()) == 0:                   #If table is empty then fill up with footprints
            datasource_in = ogr.Open(filename)    
            if datasource_in is None:
                print "Could not open shapefile file.\n"
                return n        
            #print datasource_in.GetLayerCount()
            #layer_in = datasource_in.GetLayerByName("footprints")
            
            if datasource_in.GetLayer(0) < 1:
                print "No data in the input file. \n"
                return n
            
            layer_in = datasource_in.GetLayer(0)
            layer_in.ResetReading()
            
            srs_in = layer_in.GetSpatialRef()       # the footprints spatial reference
            srs = osr.SpatialReference()            # the spatial reference of the PostGIS database (i.e. self.dbname)
            srs.ImportFromEPSG(self.epsg)
                 
            transform = 0
            if srs_in != srs:
                transform = 1
    
                       
            m = 0
            i = 0
            sql = ""
            sqlPTVA = ""
    
            for feature in layer_in:
                geom = feature.GetGeometryRef()
                #geom.wkbFlatten()               #flatten geometry if it is defined with the 25D flag
                if geom is not None and (geom.GetGeometryType() ==  ogr.wkbLineString):
					valid = geom.IsValid()
					if valid:
						if transform:
							geom.TransformTo(srs) 
						
						wkt = geom.ExportToWkt()
						self.cur.execute("SELECT b.id FROM boundary as b WHERE (ST_Contains(b.geom,ST_GeomFromText('%s',%s))) AND b.id=%s" % (wkt, self.epsg, boundaryID))        
						#if the current linestring is contained within the study boundary polygon then add to SQL command list
						if len(self.cur.fetchall()) != 0:
							sql = sql + "I %s (id,geom) VALUES (%s,ST_GeomFromText('%s',%s));\n" % (tablename,n+1, wkt, self.epsg)                             
							n += 1
					else:
						i += 1

    
                feature.Destroy()
            
            print "Number of building POLYGON footprints imported = %s" % n
            print "Number of building MULTIPOLYGON footprints NOT imported = %s" % m
            print "Number of invalid geometries = %s" % i
            
            if (sql != ""):
                self.cur.execute(sql)
                print "Creating spatial index on %s...." % tablename
                self.cur.execute("CREATE INDEX %s_idx ON %s USING GIST (geom);" % (tablename, tablename))

                self.conn.commit()

    
            datasource_in.Destroy()

        
            return n
        
        else: 
            print "Table: %s table already full." % (tablename)
            return n


    def grid_buildings(self, filename = 'triangle.poly', bounding_polygon_name = '', epsgOUT=4326,  holes = True):
     
        '''
        Create a grid with the buildings in this database using Triangle.c
        
        INPUT: bounding_polygon_name - The name of the bounding polygon (i.e. 'name' 
                                    in table 'areas')
           holes: True (default) if buildings are to be holes, False if buildings are
                  to be gridded
            
           epsgOUT:  the srs of the Triangle.c mesh
            
        
        This function extracts the buildings that lie within the input bounding polygon
        (boundary_polygon_name) and sets up triangle.c control files.
        '''
    
        if (bounding_polygon_name == ''):
            print 'ERROR: Invalid polygon name'
            return -1
        
        
        outfile = open(filename, "w")    

        #self.cur.execute("SELECT (id,ST_AsText(geom)) FROM areas WHERE name='%s';" % (boundary_polygon_name))
        self.cur.execute("SELECT geom FROM areas WHERE name='%s';" % (bounding_polygon_name))        
        boundingPolyST = self.cur.fetchall() 
        boundingPolyST = boundingPolyST[0]        
        
        self.cur.execute("SELECT ST_AsText(geom) FROM areas WHERE name='%s';" % (bounding_polygon_name))        
        boundingPolyTXT = self.cur.fetchall() 
        boundingPolyTXT = boundingPolyTXT[0][0]
        
        #get the building footprints and a point inside for all buildings that
        #are inside the bounding polygon
   
        #NOTE: ERROR in ST_PointOnSurface() function due to self intersection polygon in the buildings table
        #see: http://postgis.refractions.net/pipermail/postgis-users/2010-August/027487.html
        #fixed problem by deleting the invalid geometry.
        
    
        self.cur.execute("SELECT b.id, ST_AsText(b.geom), ST_AsText(ST_PointOnSurface(b.geom)) FROM buildings AS b \
                        WHERE (ST_Contains('%s',b.geom)) ORDER BY b.id;" % (boundingPolyST))
        buildings = self.cur.fetchall()
         
        buildingsList = []
        buildingsPtsInside = []
        
        #Convert the building polygons from TXT to python list
        for b in buildings:
            
            i = b[2].find('POINT(')
            j = b[2].find(')')
            pointInside = b[2][i+6:j]
            pointInside = pointInside.replace(',',' ').split()
            buildingsPtsInside.append(pointInside)            
            i = b[1].find('POLYGON((')
            j = b[1].find('))')
            pts = b[1][i+9:j]
            pts = pts.replace(',',' ').split()
            i = 0
            polygon = []
            while i < len(pts):
                polygon.append([float(pts[i]), float(pts[i+1])])
                i += 2      
            polygon.pop()                                                           #remove repeated last item in building polygons
            polygon = self._convert_points_(polygon,self.epsg,epsgOUT)              #convert to target SRS
            buildingsList.append(polygon)

        buildingsPtsInside = self._convert_points_(buildingsPtsInside,self.epsg,epsgOUT) 
        #convert the bounding polygon from TXT to python list
        i = boundingPolyTXT.find('POLYGON((')
        j = boundingPolyTXT.find('))')
        pts = boundingPolyTXT[i+9:j]
        pts = pts.replace(',',' ').split()       

        i=0
        boundingPoly = []
        while i < len(pts):
            boundingPoly.append([float(pts[i]), float(pts[i+1])])
            i += 2  
        
        boundingPoly.pop()                                                          #remove repeated last item in polygon
        boundingPoly = self._convert_points_(boundingPoly,self.epsg,epsgOUT)        #convert to target SRS

      
        numVerticies = len(boundingPoly)
        numHoles = len(buildingsList)

        for b in buildingsList:
            numVerticies = numVerticies + len(b)
        
        outfile.write("%s 2 0 1\n" % (str(numVerticies)))
        
        
        segments = []
        #write the vertices in the building and bounding polygon lists
        i = 1
        for pt in boundingPoly:
            #ptYshift = float(pt[1]) - 5000000
            outfile.write("%s   %s    %s   1\n" % (str(i), pt[0], pt[1]))
            if i != len(boundingPoly):
                segments.append([i,i+1,1])
                
            else:
                index = i - len(boundingPoly) + 1
                segments.append([i,index, 1])
            
            i += 1
    
        polyIndex = 2
        for poly in buildingsList:
            lengthPoly = len(poly)
            j = 1
            for pt in poly:
                #ptYshift = float(pt[1]) - 5000000
                outfile.write("%s   %s    %s    %s\n" % (str(i), pt[0], pt[1], str(polyIndex)))
               
                if j != lengthPoly:
                    segments.append([i,i+1,polyIndex])
                    
                else:
                    index = i - lengthPoly + 1
                    segments.append([i,index,polyIndex])
                
                j += 1
                i += 1
            
            polyIndex += 1            
        
        #write the segements
        outfile.write("%s 1\n" % (str(len(segments))))
        i = 1
        for seg in segments:
            outfile.write("%s   %s   %s   %s\n" % (str(i), str(seg[0]), str(seg[1]), str(seg[2])))
            i += 1
    
    
        #write the holes
        if(holes == True):
            outfile.write("%s\n" % (str(len(buildingsPtsInside))))
            i = 1
            for pt in buildingsPtsInside:
                #ptYshift = float(pt[1]) - 5000000
                outfile.write("%s   %s   %s\n" % (str(i), str(pt[0]), str(pt[1])))
                i += 1
        else:
            outfile.write("0\n")
    
    
        
        outfile.close()
        
        print len(buildingsPtsInside)
        print len(buildingsList)
        print buildingsPtsInside[0]
        print buildingsList[0]


    def copy_record(self,tablename = 'areas',id = 1):
        '''
        Creates a copy of the geom with supplied 'id' in table 'tablename'
        
        INSERTS the geometry back into the supplied table 'tablename'
        
        
        '''
        
        #get the highest pkey value and set the new record to one larger
        self.cur.execute("SELECT pkey FROM %s ORDER BY pkey;" % (tablename))
        ids = self.cur.fetchall() 
        pkeyNew =  ids[len(ids)-1][0] + 1
        
        self.cur.execute("SELECT ST_AsText(geom) FROM %s WHERE id=%s;" % (tablename, id))        
        geom = self.cur.fetchall() 
        geom = geom[0][0]
        self.cur.execute("INSERT INTO %s (pkey,id,geom,name) VALUES (%s, %s, ST_GeomFromText('%s',%s), 'copy');" % (tablename,pkeyNew,pkeyNew, geom, self.epsg))     
        self.conn.commit()
        



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