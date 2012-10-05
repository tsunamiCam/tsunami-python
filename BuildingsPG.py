#Import the required modules


import psycopg2 #@UnresolvedImport
#import errorcodes 


import sys
from osgeo import ogr #@UnresolvedImport
from osgeo import osr #@UnresolvedImport

from subprocess import call
from numpy import *


class BuildingsPG:
    """
    class BuildingsPG
    
    This class opens and manages an existing PostGIS database of buildings
    
    """

    
    def __init__(self, 
                 dbname, 
                 user = "tbone",
                 epsg = 4326):

        self.dbname = dbname            #the database name
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
            self.conn = psycopg2.connect(database=self.dbname, user=self.user);
        except Exception, e:
            if e[0].find("could not connect to server") >= 0:
                print "Postgis database server not running. Please start postgres\n"
                sys.exit()   
        
            elif e[0].find("does not exist") >= 0:
                print "Creating new BUILDINGS PostGIS database: \"%s\"" % self.dbname            
                p = call(['createdb', self.dbname])
                print p
                p = call(['psql','-d',self.dbname,'-f',pg_1])
                print p
                p = call(['psql','-d',self.dbname,'-f',pg_2])
                print p
                p = call(['psql','-d',self.dbname,'-f',pg_3])  
                print p
                #Try connecting to the database again
                try:
                    self.conn = psycopg2.connect(database=self.dbname, user=self.user);
                except:
                    "ERROR: cannot connect or create database: \"%s\"" % self.dbname 
                    sys.exit()
                else: 
                    self.cur = self.conn.cursor()
        
                    #------------------CREATE TABLES--------------------        
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
                    				ALTER TABLE buildings ADD COLUMN comment varchar(128);")        	
                    
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
                    self.conn.commit()
        else:
            self.cur = self.conn.cursor()

        
        #cur.close()
        #conn.close()
    
    def __del__ (self):
        """
        Class deconstructor
        
        Closes the connection to the PostGIS database
        """
        
        self.cur.close()
        self.conn.close()

        
    def calculate_ptva(self):
        """
        Using the data from the buildings table - calculate the PTVA RVI score for each building
        
        
        NOTE: the ptva model is not properly represented at present.  Trying to get good numbers for presentation
        
        SV is not scaled properly
        RVI score doesn't include WV (water vulnerability)
        bv scaling changed for effect
        """
        
        #Create the ptva column in the buildings table
        try: self.cur.execute("ALTER TABLE buildings ADD COLUMN ptva float DEFAULT 0;")
        except:
            print "calculate_ptva: Can't add ptva column to the buildings table.  Already exists.  This is just a warning. Trying to write..."
            self.conn.commit()
                
        
        #BV constants
        bv_w1 = 0.236
        bv_w2 = 0.189
        bv_w3 = 0.149
        bv_w4 = 0.142
        bv_w5 = 0.121
        bv_w6 = 0.109
        bv_w7 = 0.054
        
        
        #Prot constants
        prot_w1 = 0.332
        prot_w2 = 0.243
        prot_w3 = 0.243
        prot_w4 = 0.183
        
        
        self.cur.execute("SELECT id, s, m, g, f, so, mo, pc, prot_br, prot_nb, prot_sw, prot_w, fd_max FROM buildings ORDER BY id;" )
        attributes = self.cur.fetchall()
        for row in attributes:

            bv = bv_w1*row[1] + bv_w2*row[2] + bv_w3*row[3] + bv_w4*row[4] + bv_w5*row[5] + bv_w6*row[6] + bv_w7*row[7]
            prot = prot_w1*row[8] + prot_w2*row[9] + prot_w3*row[10] + prot_w4*row[11]
            ex = row[12]
            
            
            print "id = %s, ex = %s, bv = %s, prot = %s" % (str(row[0]), str(ex), str(bv), str(prot))

            #Scale prot to be between 1 and 5
            if prot >= 0 and prot < 0.2 : prot = 1
            elif prot >= 0.2 and prot < 0.4 : prot = 2
            elif prot >= 0.4 and prot < 0.6 : prot = 3
            elif prot >= 0.6 and prot < 0.8 : prot = 4
            elif prot >= 0.8 and prot < 1 : prot = 5

            #Scale bv to be between 1 and 5        
            #NOTE: THESE AREN't THE PTVA VALUES - changed for effect
            if bv >= -1 and bv < -0.6 : bv = 1
            elif bv >= -0.6 and bv < 0 : bv = 2
            elif bv >= 0 and bv < 0.2 : bv = 3
            elif bv >= 0.2 and bv < 0.6 : bv = 4
            elif bv >= 0.6 and bv < 1 : bv = 5

            #Scale ex to be between 1 and 5
            if ex > 0 and ex < 1 : ex = 1
            elif ex >= 1 and ex < 2 : ex = 2
            elif ex >= 2 and ex < 3 : ex = 3
            elif ex >= 3 and ex < 4 : ex = 4
            elif ex >= 4: ex = 5  
            else: ex = 0      
            
            
            print "ex = %s, bv = %s, prot = %s" % (str(ex), str(bv), str(prot))
            #  print "ex = %s, bv = %s, prot = %s, sv = %s, WV = %s, RVI = %s" % (str(ex), str(bv), str(prot), str(SV), str(WV), str(RVI))

            
            SV = bv*ex*prot
            
            
            #comment out to make numbers look better
            #Scale ex to be between 1 and 5
#            if SV >= 1 and SV < 25 : SV = 1
#            elif SV >= 25 and SV < 50 : SV = 2
#            elif SV >= 50 and SV < 75 : SV = 3
#            elif SV >= 75 and SV < 100 : SV = 4
#            elif SV >= 100 and SV < 125 : SV = 5
            
            if row[1] == 1: numLevels = 1
            elif row[1] == 0.5: numLevels = 2
            elif row[1] == 0: numLevels = 3
            elif row[1] == -0.5: numLevels = 4
            elif row[1] == -1: numLevels = 5

            if row[12] == 0: numInundatedLevels = 0
            elif (row[12] > 0) and (row[12] < 2): numInundatedLevels = 1
            elif (row[12] >= 2) and (row[12] < 4): numInundatedLevels = 2
            elif (row[12] >= 4) and (row[12] < 6): numInundatedLevels = 3
            elif (row[12] >= 6) and (row[12] < 8): numInundatedLevels = 4
            
            WV = numInundatedLevels/numLevels
            
            #Scale WV to be between 1 and 5
            if WV >= 0 and WV < 0.2 : WV = 1
            elif WV >= 0.2 and WV < 0.4 : WV = 2
            elif WV >= 0.4 and WV < 0.6 : WV = 3
            elif WV >= 0.6 and WV < 0.8 : WV = 4
            elif WV >= 0.8 and WV < 1 : WV = 5
            
#            RVI = (2/3)*SV + (1/3)*WV
            RVI = SV            #Not true just experimenting

        
            self.cur.execute("UPDATE buildings SET ptva = %s WHERE id = %s;" % (RVI,row[0]))
        
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
                                ALTER TABLE buildings ADD COLUMN wall_height float DEFAULT 0; \
                                ALTER TABLE buildings ADD COLUMN water_dept float DEFAULT 0; \
                                ALTER TABLE buildings ADD COLUMN inund_lev float DEFAULT 0; \
                                ALTER TABLE buildings ADD COLUMN elevation float DEFAULT 0; \
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
                print n
            
            return pointsOUT
        
        else:
            return pointsIN


    def export_line_as_SHP(self, filename_out,line_id, epsgOUT):
        """
        Export a line as a SHP file
        
        """
        
        self.cur.execute("SELECT ST_AsText(geom) FROM lines WHERE id = %s;" % (line_id))
        lineST = self.cur.fetchall() 
        geom = lineST[0][0]
        
        lineStr = [] 
        if(geom.find('LINESTRING') >= 0):
            i = geom.find('LINESTRING(')
            j = geom.find(')')
            pts = geom[i+11:j]
            pts = pts.replace(',',' ').split()
            i = 0
            while i < len(pts):
                lineStr.append([float(pts[i]), float(pts[i+1])])
                i += 2      
 
        srs_in = osr.SpatialReference()
        srs_in.ImportFromEPSG(self.epsg)
        srs_out = osr.SpatialReference()            
        srs_out.ImportFromEPSG(epsgOUT) 
        driver = ogr.GetDriverByName('ESRI Shapefile')
        outfile = driver.CreateDataSource(filename_out)
        
        layerLinestring = outfile.CreateLayer("Linestring", srs_out, ogr.wkbLineString)    
        featureLinestring = ogr.Feature(layerLinestring.GetLayerDefn())
        line = ogr.Geometry(type=ogr.wkbLineString)        

        for p in lineStr:
            line.AddPoint(p[0],p[1])

        line.AssignSpatialReference(srs_in)
        line.TransformTo(srs_out)
        featureLinestring.SetGeometry(line)
        layerLinestring.CreateFeature(featureLinestring) 
        outfile = None

        featureLinestring.Destroy()
                
                
    def export_area_as_SHP(self, filename_out,area_id, epsgOUT):
        """
        Export an area as a SHP file
        Converts the area from a linestring to a polygon
        
        """
	
		
        self.cur.execute("SELECT ST_AsText(geom) FROM areas WHERE id = %s;" % (area_id))
        areaST = self.cur.fetchall() 
        geom = areaST[0][0]
		
        area = []
        if (geom.find('POLYGON') >= 0):
            i = geom.find('POLYGON((')
            j = geom.find('))')
            pts = geom[i+9:j]
            pts = pts.replace(',',' ').split()
            i = 0
            while i < len(pts):
                area.append([float(pts[i]), float(pts[i+1])])
                i += 2

		
		srs_in = osr.SpatialReference()
		srs_in.ImportFromEPSG(self.epsg)
		srs_out = osr.SpatialReference()			
		srs_out.ImportFromEPSG(epsgOUT)
		
		#write the polygon data to KML
        driver = ogr.GetDriverByName('ESRI Shapefile')
        outfile = driver.CreateDataSource(filename_out)
        layerPolygons = outfile.CreateLayer("Polygons", srs_out, ogr.wkbPolygon)	
        featurePolygons = ogr.Feature(layerPolygons.GetLayerDefn())
        '''
        layerLinestring = outfile.CreateLayer("Linestring", srs_out, ogr.wkbLineString)    
        featureLinestring = ogr.Feature(layerLinestring.GetLayerDefn())
        line = ogr.Geometry(type=ogr.wkbLineString)
        '''
        #print outfile.GetLayerCount()
        
        #NOTE: To create a polygon you first need to create a ring
        #		then add it to the polygon
        ptCount = len(area)
        firstPt = area[0]
        lastPt = area[ptCount-1]
        polygon = ogr.Geometry(type=ogr.wkbPolygon)			#create polygon
        ring = ogr.Geometry(type=ogr.wkbLinearRing)			#create ring
		
        #Only convert the feature to polygon if is a closed linestring
        if firstPt[0] == lastPt[0] and firstPt[1] == lastPt[1]:
            print "Valid LINESTRING"	  
            i = 0
            while i < ptCount-1:
                p = area[i]
                ring.AddPoint(p[0],p[1])
                i += 1
            ring.CloseRings()
            polygon.AddGeometry(ring)						#add ring to polygon	
            polygon.AssignSpatialReference(srs_in)
            polygon.TransformTo(srs_out)
            featurePolygons.SetGeometry(polygon)
            layerPolygons.CreateFeature(featurePolygons) 
        
        '''
        for p in area:
            line.AddPoint(p[0],p[1])

        line.AssignSpatialReference(srs_in)
        line.TransformTo(srs_out)
        featureLinestring.SetGeometry(line)
        layerLinestring.CreateFeature(featureLinestring) 
        
        featureLinestring.Destroy()
        
        '''
        
        outfile = None
        featurePolygons.Destroy()

    def line_to_dxf(self, dxf_filename,line_id, epsgOUT):
        
        self.cur.execute("SELECT ST_AsText(geom) FROM lines WHERE id = %s;" % (line_id))
        lineST = self.cur.fetchall() 
        geom = lineST[0][0]
        
        lineStr = [] 
        if(geom.find('LINESTRING') >= 0):
            i = geom.find('LINESTRING(')
            j = geom.find(')')
            pts = geom[i+11:j]
            pts = pts.replace(',',' ').split()
            i = 0
            while i < len(pts):
                lineStr.append([float(pts[i]), float(pts[i+1])])
                i += 2      


        # set the spatial reference - of the polygon data
        coordsPG = osr.SpatialReference()
        coordsPG.ImportFromEPSG(self.epsg)
        coordsOUT = osr.SpatialReference()
        coordsOUT.ImportFromEPSG(int(epsgOUT))
    
        #write the polygon data to KML
        driver = ogr.GetDriverByName('DXF')
        dxfData = driver.CreateDataSource(dxf_filename)
        layerPolygons =  dxfData.CreateLayer("polygons", coordsOUT, ogr.wkbLineString)    
        featurePolygons = ogr.Feature(layerPolygons.GetLayerDefn())
        
        #write the polygon data to the KML file
        line = ogr.Geometry(type=ogr.wkbLineString)
        for pt in lineStr:
            line.AddPoint(pt[0],pt[1])
              
        line.AssignSpatialReference(coordsPG)
        line.TransformTo(coordsOUT)
        featurePolygons.SetGeometry(line)
        layerPolygons.CreateFeature(featurePolygons)
        dxfData.Destroy()
        featurePolygons.Destroy()


    def area_to_dxf(self, dxf_filename,area_id, epsgOUT):
		
		self.cur.execute("SELECT ST_AsText(geom) FROM areas WHERE id = %s;" % (area_id))
		areaST = self.cur.fetchall() 
		geom = areaST[0][0]
		
		polygon = []
		if (geom.find('POLYGON') >= 0):
			i = geom.find('POLYGON((')
			j = geom.find('))')
			pts = geom[i+9:j]
			pts = pts.replace(',',' ').split()
			i = 0
			while i < len(pts):
				polygon.append([float(pts[i]), float(pts[i+1])])
				i += 2
		
		print polygon

		# set the spatial reference - of the polygon data
		coordsPG = osr.SpatialReference()
		coordsPG.ImportFromEPSG(self.epsg)
		coordsOUT = osr.SpatialReference()
		coordsOUT.ImportFromEPSG(int(epsgOUT))
	
		#write the polygon data to KML
		driver = ogr.GetDriverByName('DXF')
		dxfData = driver.CreateDataSource(dxf_filename)
		layerPolygons =  dxfData.CreateLayer("polygons", coordsOUT, ogr.wkbLineString)	
		featurePolygons = ogr.Feature(layerPolygons.GetLayerDefn())
		
		#write the polygon data to the KML file
		line = ogr.Geometry(type=ogr.wkbLineString)
		for pt in polygon:
			line.AddPoint(pt[0],pt[1])
			  
		line.AssignSpatialReference(coordsPG)
		line.TransformTo(coordsOUT)
		featurePolygons.SetGeometry(line)
		layerPolygons.CreateFeature(featurePolygons)
		dxfData.Destroy()
		featurePolygons.Destroy()
		
			
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
        
    def add_area(self,ply_in,epsgIN,id):
        """
        """  
        ply_file = open(ply_in, "r").readlines() 
        boundingPoly = []
        for row in ply_file:
            row =  row.split()
            boundingPoly.append(row)
            
            #print "row = %s" % row        
        
        boundingPoly = self._convert_points_(boundingPoly, epsgIN,epsgOUT=self.epsg) 
        '''    
        self.cur.execute("DROP TABLE IF EXISTS boundary;")
        self.cur.execute("CREATE TABLE boundary (id int4, description varchar(64));")
        self.cur.execute("SELECT AddGeometryColumn('boundary', 'geom', %s, 'POLYGON', 2);" % self.epsg)
        '''
        geom = "ST_GeomFromText('POLYGON(("
        numPts = len(boundingPoly)
        i = 0
        while i < (numPts - 1):
            ptText = "%s %s, " % (str(boundingPoly[i][0]),str(boundingPoly[i][1]))
            geom = geom + ptText
            i += 1
            
        ptText = "%s %s))',%s)" % (str(boundingPoly[i][0]),str(boundingPoly[i][1]),self.epsg)
        geom = geom + ptText
    	print geom
        sql = "INSERT INTO areas (pkey,geom) VALUES (%s, %s);" % (id,geom)
        self.cur.execute(sql)  
        self.conn.commit()

        #self.cur.execute("SELECT ST_AsText(geom) FROM areas;")  
        #boundary = self.cur.fetchall() 
        #self.geom_to_kml(boundary[0],kml_file="boundary.kml", epsgIN=self.epsg)

        #return
  
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
							sql = sql + "INSERT INTO %s (id,geom) VALUES (%s,ST_GeomFromText('%s',%s));\n" % (tablename,n+1, wkt, self.epsg)                             
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



    def replace_buildings_from_kml(self,kml_filename, area_id):
  
  
  
  
        datasource_in = ogr.Open(kml_filename)    
        if datasource_in is None:
            print "Could not open KML file.\n"
            return n        
        
        if datasource_in.GetLayer(0) < 1:
            print "No data in the input file. \n"
            return n
        
        layer_in = datasource_in.GetLayer(0)
        layer_in.ResetReading()
        
        srs_in = layer_in.GetSpatialRef()       # the footprints spatial reference
        srs = osr.SpatialReference()            # the spatial reference of the PostGIS database (i.e. self.dbname)
        srs.ImportFromEPSG(self.epsg)


        #self.cur.execute("SELECT (id,ST_AsText(geom)) FROM areas WHERE name='%s';" % (boundary_polygon_name))
        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (area_id))        
        area = self.cur.fetchall() 
        area = area[0]        
        
        print area
        
    
        self.cur.execute("DELETE FROM buildings \
                        WHERE (ST_Contains('%s',geom));" % (area))

        self.cur.execute("DELETE FROM buildings \
                        WHERE id > 4999;")


        self.conn.commit()
             
        transform = 0
        if srs_in != srs:
            transform = 1
        
        idStart = 5000
                   
        m = 0
        i = 0
        sql = ""
        sqlPTVA = ""
        j = 0
        p=0
        for feature in layer_in:
            geom = feature.GetGeometryRef()
            p+=1
            #geom.wkbFlatten()               #flatten geometry if it is defined with the 25D flag
            if transform:
                geom.TransformTo(srs)             
            
            
            wkt = geom.ExportToWkt()

            #wkt = wkt.wktFlatten()

            self.cur.execute("INSERT INTO buildings (pkey,geom,visibility_after,damage,prot_nb,prot_sw,prot_w,so,mo) \
                             VALUES (%s, ST_Force_2D(ST_GeomFromText('%s','4326')), \
                             3,7,1,1,1,0.5,0.5);" % (idStart,wkt))     
            
            idStart += 1

            feature.Destroy()


            
            '''
            if geom is not None and (geom.GetGeometryType() ==  ogr.wkbPolygon or geom.GetGeometryType() ==  ogr.wkbMultiPolygon):
                
                if (geom.GetGeometryType() ==  ogr.wkbPolygon):
                    valid = geom.IsValid()
                    if valid:
                        if transform:
                            geom.TransformTo(srs) 
                        j+=1
                    else:
                        i += 1
                elif (geom.GetGeometryType() ==  ogr.wkbMultiPolygon):
                    print "Warning: footprints layer contains MULTIPOLYGONS."
                    m += 1
    
            feature.Destroy()
            '''
        
        #if (geom.GetGeometryType() == ogr.wkbPolygon):
        #    print "YES" 


        self.conn.commit()

        
        print wkt
        datasource_in.Destroy()

    
        print i,j,m,p
        
        
        



    def grid_buildings(self, filename, bounding_area_id, epsgOUT,  holes = True):
     
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
        
        outfile = open(filename, "w")    

        #self.cur.execute("SELECT (id,ST_AsText(geom)) FROM areas WHERE name='%s';" % (boundary_polygon_name))
        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (bounding_area_id))        
        boundingPolyST = self.cur.fetchall() 
        boundingPolyST = boundingPolyST[0]        
        
        self.cur.execute("SELECT ST_AsText(geom) FROM areas WHERE id='%s';" % (bounding_area_id))        
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
         
        intersectingBuildings = 0
    	for poly in buildings: #check if there are any intersecting buildings   
			self.cur.execute("SELECT b.id FROM buildings AS b \
							WHERE ST_Intersects(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" \
							% (poly[1], self.epsg))
			intersect = self.cur.fetchall()
			if len(intersect) > 1:
				print "Building = %s intersects Building = %s" % (poly[0],intersect[1])
				intersectingBuildings = 1

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
        
        #write the segments
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



    def grid_buildings2(self, filename, bounding_area_id,outer_area_id,inner_area_id, epsgOUT,  holes = True):
     
        '''
        Create a grid with the buildings which area inside outer area
        
        INPUT: filename - The name of the bounding polygon (i.e. 'name' 
                                    in table 'areas')
           holes: True (default) if buildings are to be holes, False if buildings are
                  to be gridded
            
           epsgOUT:  the srs of the Triangle.c mesh
            
        
        This function extracts the buildings that lie within the input bounding polygon
        (boundary_polygon_name) and sets up triangle.c control files.
        '''
        
        outfile = open(filename, "w")    

        
        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (bounding_area_id))        
        boundingArea = self.cur.fetchall()[0][0]
        
        self.cur.execute("SELECT ST_AsText(geom) FROM areas WHERE id='%s';" % (bounding_area_id))        
        boundingPolyTXT = self.cur.fetchall() 
        boundingPolyTXT = boundingPolyTXT[0][0]

        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (outer_area_id))        
        outerArea = self.cur.fetchall()[0][0]
        
        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (inner_area_id))        
        innerArea = self.cur.fetchall()[0][0]
        
        
        self.cur.execute("SELECT b.id, ST_AsText(b.geom), ST_AsText(ST_PointOnSurface(b.geom)) FROM buildings AS b \
                        WHERE \
                        (ST_Contains('%s',b.geom) AND NOT ST_Contains('%s',b.geom)) \
                        OR \
                        ST_Contains('%s',b.geom) \
                        ORDER BY b.id;" % (boundingArea, outerArea, innerArea))
        buildings = self.cur.fetchall()


       #get the building footprints and a point inside for all buildings that
        #are inside the bounding polygon
   
        #NOTE: ERROR in ST_PointOnSurface() function due to self intersection polygon in the buildings table
        #see: http://postgis.refractions.net/pipermail/postgis-users/2010-August/027487.html
        #fixed problem by deleting the invalid geometry.

         
        intersectingBuildings = 0
        for poly in buildings: #check if there are any intersecting buildings   
            self.cur.execute("SELECT b.id FROM buildings AS b \
                            WHERE ST_Intersects(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" \
                            % (poly[1], self.epsg))
            intersect = self.cur.fetchall()
            if len(intersect) > 1:
                print "Building = %s intersects Building = %s" % (poly[0],intersect[1])
                intersectingBuildings = 1

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
        
        #write the segments
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



    def grid_buildings3(self, filename, bounding_area_id,outer_area_id,inner_area_id, epsgOUT):
     
        '''
        Create a grid with a mixture of buildings as holes and buildings that are gridded
        
        BUILDINGS as HOLES: All buildings which are inside both bounding_area_id & inner_area_id 
        GRIDDED BUILDINGS: ALL buildings which are inside outer_area_id AND outside inner_area
        
        INPUT: filename - name of the triangle input file to be written
        
           epsgOUT:  the srs of the Triangle.c mesh
            
        
        This function extracts the buildings that lie within the input bounding polygon
        (boundary_polygon_name) and sets up triangle.c control files.
        '''
        
        outfile = open(filename, "w")    

        
        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (bounding_area_id))        
        boundingArea = self.cur.fetchall()[0][0]
        
        self.cur.execute("SELECT ST_AsText(geom) FROM areas WHERE id='%s';" % (bounding_area_id))        
        boundingPolyTXT = self.cur.fetchall() 
        boundingPolyTXT = boundingPolyTXT[0][0]

        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (outer_area_id))        
        outerArea = self.cur.fetchall()[0][0]
        
        self.cur.execute("SELECT geom FROM areas WHERE id='%s';" % (inner_area_id))        
        innerArea = self.cur.fetchall()[0][0]
        
        
        self.cur.execute("SELECT b.id, ST_AsText(b.geom), ST_AsText(ST_PointOnSurface(b.geom)) FROM buildings AS b \
                        WHERE \
                        (ST_Contains('%s',b.geom) AND NOT ST_Contains('%s',b.geom)) \
                        OR \
                        ST_Contains('%s',b.geom) \
                        ORDER BY b.id;" % (boundingArea, outerArea, innerArea))
        buildingsHoles = self.cur.fetchall()

        self.cur.execute("SELECT b.id, ST_AsText(b.geom), ST_AsText(ST_PointOnSurface(b.geom)) FROM buildings AS b \
                        WHERE \
                        (ST_Contains('%s',b.geom) AND NOT ST_Contains('%s',b.geom)) \
                        ORDER BY b.id;" % (outerArea, innerArea))
        buildingsGridded = self.cur.fetchall()


       #get the building footprints and a point inside for all buildings that
        #are inside the bounding polygon
   
        #NOTE: ERROR in ST_PointOnSurface() function due to self intersection polygon in the buildings table
        #see: http://postgis.refractions.net/pipermail/postgis-users/2010-August/027487.html
        #fixed problem by deleting the invalid geometry.

         
        intersectingBuildings = 0
        for poly in buildingsGridded: #check if there are any intersecting buildings   
            self.cur.execute("SELECT b.id FROM buildings AS b \
                            WHERE ST_Intersects(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" \
                            % (poly[1], self.epsg))
            intersect = self.cur.fetchall()
            if len(intersect) > 1:
                print "Building = %s intersects Building = %s" % (poly[0],intersect[1])
                intersectingBuildings = 1


        intersectingBuildings = 0
        for poly in buildingsHoles: #check if there are any intersecting buildings   
            self.cur.execute("SELECT b.id FROM buildings AS b \
                            WHERE ST_Intersects(ST_GeomFromText('%s',%s),b.geom) ORDER BY b.id;" \
                            % (poly[1], self.epsg))
            intersect = self.cur.fetchall()
            if len(intersect) > 1:
                print "Building = %s intersects Building = %s" % (poly[0],intersect[1])
                intersectingBuildings = 1


        buildingsListHoles = []
        buildingsPtsInsideHoles = []
        buildingsListGridded = []
        buildingsPtsInsideGridded = []


        #------------------------------------------------------------------------------
        #Convert the building polygons that will be HOLES from TXT to python list
        #------------------------------------------------------------------------------
        for b in buildingsHoles:
            
            i = b[2].find('POINT(')
            j = b[2].find(')')
            pointInside = b[2][i+6:j]
            pointInside = pointInside.replace(',',' ').split()
            buildingsPtsInsideHoles.append(pointInside)            
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
            buildingsListHoles.append(polygon)

        buildingsPtsInsideHoles = self._convert_points_(buildingsPtsInsideHoles,self.epsg,epsgOUT) 

        #------------------------------------------------------------------------------
        #Convert the building polygons that will be HOLES from TXT to python list
        #------------------------------------------------------------------------------
        for b in buildingsGridded:
            
            i = b[2].find('POINT(')
            j = b[2].find(')')
            pointInside = b[2][i+6:j]
            pointInside = pointInside.replace(',',' ').split()
            buildingsPtsInsideGridded.append(pointInside)            
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
            buildingsListGridded.append(polygon)

        buildingsPtsInsideGridded = self._convert_points_(buildingsPtsInsideGridded,self.epsg,epsgOUT) 
        #------------------------------------------------------------------------------




        #convert the bounding polygon from TXT to python list
        #------------------------------------------------------------------------------

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
        numHoles = len(buildingsListHoles+buildingsListGridded)

        for b in (buildingsListHoles+buildingsListGridded):
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
        for poly in (buildingsListHoles+buildingsListGridded):
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
        
        #write the segments
        outfile.write("%s 1\n" % (str(len(segments))))
        i = 1
        for seg in segments:
            outfile.write("%s   %s   %s   %s\n" % (str(i), str(seg[0]), str(seg[1]), str(seg[2])))
            i += 1
    
    
        #write the holes
        outfile.write("%s\n" % (str(len(buildingsPtsInsideHoles))))
        i = 1
        for pt in buildingsPtsInsideHoles:
            outfile.write("%s   %s   %s\n" % (str(i), str(pt[0]), str(pt[1])))
            i += 1
        outfile.write("0\n")
    
    
        
        outfile.close()
    
    
    

    def init_buildings_block_table(self):
        self.cur.execute("DROP TABLE IF EXISTS buildings_block;")
        self.cur.execute("CREATE TABLE buildings_block (id int4);")    
        self.cur.execute("SELECT AddGeometryColumn('buildings_block', 'geom', %s, 'POLYGON', 2);" % self.epsg) 
    
    def add_buildings_block(self, xCorner,yCorner, angle,blockX,blockY,bldgX,bldgY, nx, ny):
        '''
        
        
        
        [x0,y0] - building block orgin
        angle - angle of the building block (i.e. angle of the grouping of buildings)
        blockLength - length of the building block
        blockWidth - width of the building block
        L - Length of each building
        W - Width of each building
        nx
        ny
        '''
#        self.cur.execute("DROP TABLE IF EXISTS buildings_block;")
#        self.cur.execute("CREATE TABLE buildings_block (id int4);")    
#        self.cur.execute("SELECT AddGeometryColumn('buildings_block', 'geom', %s, 'POLYGON', 2);" % self.epsg) 

        #calculate the space between buildings
        spaceX = (blockX - nx*bldgX)/(nx-1)
        spaceY = (blockY - ny*bldgY)/(ny-1)
        
        print spaceX
        print spaceY
        
        #the rotatation matrix
        angle = angle*math.pi/180        
        R = []
        R.append([math.cos(angle), -math.sin(angle)]) 
        R.append([math.sin(angle), math.cos(angle)]) 

        #create the base building footprint
        bldg = []
        bldg.append([0,0])
        bldg.append([bldgX,0])
        bldg.append([bldgX,bldgY])
        bldg.append([0,bldgY])
        bldg.append([0,0])
        
        #rotate the base buildign footprint
        bldgR = []
        i = 0
        for pt in bldg:
            x = pt[0]*R[0][0] + pt[1]*R[0][1]
            y = pt[0]*R[1][0] + pt[1]*R[1][1]
            bldgR.append([x,y])
            i+=1


        ix=0
        iy=0
        while iy < ny:
            ix =0
            while ix < nx:
                #origin of the building (block coordinates)
                x0 = ix*bldgX+spaceX*ix
                y0 = iy*bldgY+spaceY*iy
                #rotate the orgin
                xr = x0*R[0][0] + y0*R[0][1]
                yr = x0*R[1][0] + y0*R[1][1]
                
                x0 = xr
                y0 = yr

                i=0
                bldgOUT = []
                for pt in bldgR:
                    x = pt[0] + x0 + xCorner
                    y = pt[1] + y0 + yCorner
                    bldgOUT.append([x,y])

                geom = "POLYGON(("
                i = 0
                while i < len(bldgOUT):
                    pt = bldgOUT[i]
                    if i == 4:
                        ptText = "%s %s))" % (pt[0],pt[1])
                    else:
                        ptText = "%s %s, " % (pt[0],pt[1])
                    geom = geom + ptText
                    i+=1
                
                ix+=1
                
                self.cur.execute("INSERT INTO buildings_block (geom) VALUES (ST_GeomFromText('%s',%s));" % (geom,self.epsg))         
 
            iy+=1

        self.conn.commit()
          





         
        
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
        



    def get_record(self,source_db,tablename = 'areas',id = 1):
        '''
        Creates a copy of the geom with supplied 'id' in table 'tablename'
        
        INSERTS the geometry back into the supplied table 'tablename'
        
        
        '''
        
        connTEMP = psycopg2.connect(database=source_db, user=self.user);
        curTEMP = connTEMP.cursor()
        curTEMP.execute("SELECT ST_AsText(geom) FROM %s WHERE id=%s;" % (tablename, id))        
        geom = curTEMP.fetchall()[0][0]
        self.cur.execute("INSERT INTO %s (pkey,id,geom) VALUES (%s, %s, ST_GeomFromText('%s',%s));" % (tablename,id,id, geom, self.epsg))     
        self.conn.commit()
        


    
    
    def transform_table(self, tablename, epsgOUT):
        
        tablename_new = tablename + "_" + str(epsgOUT)
        t = tablename_new
        self.cur.execute("DROP TABLE IF EXISTS %s;" % (tablename_new))




        try: 
            self.cur.execute("CREATE TABLE %s (pkey integer PRIMARY KEY,id int4);" % (tablename_new))     
        except: 
            print "Table \"%s\" already exists." % tablename_new
        else: 
            self.cur.execute("SELECT AddGeometryColumn('%s', 'geom', %s, 'POLYGON', 2);" % (tablename_new,epsgOUT))
            #ADD extra PTVA columns to the 'buildings' table
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
                            ALTER TABLE %s ADD COLUMN comment varchar(128);" % (t,t,t,t,t,t,t,t,t,t,t,t) )            
            
            #CREATE RULE for when new records are added to the buildings table
            #This rules makes the ID of new building record equal to PKEY
            self.cur.execute("CREATE RULE force_id AS ON INSERT TO %s \
                            DO ALSO \
                            UPDATE %s SET id = NEW.pkey WHERE pkey = NEW.pkey;" % (t,t))

        #self.cur.execute("CREATE RULE force_srs AS ON INSERT TO %s \
        #                  DO INSTEAD \
        #                  INSERT INTO %s (geom) VALUES ST_Transform(NEW.geom,%s)" % (t,t,epsgOUT))       
        
        
        self.conn.commit() 
        
        self.cur.execute("SELECT pkey, ST_AsText(ST_Transform(geom,%s)) FROM %s;" % (epsgOUT,tablename))
        #self.cur.execute("SELECT pkey, ST_AsText(geom) FROM %s;" % (tablename))
        
        geom_transformed = self.cur.fetchall() 
        
        sql = ""
        for g in geom_transformed:
            sql = sql + "INSERT INTO %s (pkey, geom) VALUES (%s, ST_GeomFromText('%s',%s));" % (tablename_new, g[0], g[1],epsgOUT)
        
        self.cur.execute(sql)
        self.conn.commit() 
        


    def add_footprints_manly(self, filename):
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
                #geom.WkbFlatten()               #flatten geometry if it is defined with the 25D flag
                #if geom is not None and (geom.GetGeometryType() ==  ogr.wkbPolygon or geom.GetGeometryType() ==  ogr.wkbMultiPolygon):
                    
#                    if (geom.GetGeometryType() ==  ogr.wkbPolygon):

                geom.FlattenTo2D()
                valid = geom.IsValid()

                if valid:
                    if transform:
                        geom.TransformTo(srs) 
                    
                    wkt = geom.ExportToWkt()
                    
                    if (geom.GetGeometryType() ==  ogr.wkbPolygon):

                        #if the current footprint is contained within the study boundary polygon then add to SQL command list
                        #sql = sql + "INSERT INTO buildings (pkey,id,geom) VALUES (%s,%s,ST_GeomFromText('%s',%s));\n" % (n+1,n+1, wkt, self.epsg)                             

                        n+=1 
                        s = feature.GetFieldAsString('s')
                        mat = feature.GetFieldAsString('m')
                        g = feature.GetFieldAsString('g')
                        f = feature.GetFieldAsString('f')
                        so = feature.GetFieldAsString('so')
                        mo = feature.GetFieldAsString('mo')
                        pc = feature.GetFieldAsString('pc')
                        prot_sw = feature.GetFieldAsString('prot_sw')
                        prot_w = feature.GetFieldAsString('prot_w')
                        prot_nb = feature.GetFieldAsString('prot_nb')
                        prot_br = feature.GetFieldAsString('prot_br')
                        wall_heigh = feature.GetFieldAsString('wall_heigh')
                        water_dept = feature.GetFieldAsString('water_dept')
                        inund_lev = feature.GetFieldAsString('inund_lev')
                        elevation = feature.GetFieldAsString('elevation')
                        sql = sql + "INSERT INTO buildings (pkey,id,geom,s,m,g,f,so,mo,pc,prot_sw,prot_w,prot_nb,prot_br,wall_height,water_dept,inund_lev,elevation) \
                        VALUES (%s,%s,ST_GeomFromText('%s',%s),%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);\n" % (n+1,n+1, wkt, self.epsg,s,mat,g,f,so,mo,pc,prot_sw,prot_w,prot_nb,prot_br,wall_heigh,water_dept,inund_lev,elevation)

                    elif (geom.GetGeometryType() ==  ogr.wkbMultiPolygon):
                        print "Warning: footprints layer contains MULTIPOLYGONS."
                        m+=1
                else:
					i += 1
 #                   elif (geom.GetGeometryType() ==  ogr.wkbMultiPolygon):
 #                       print "Warning: footprints layer contains MULTIPOLYGONS."
 #                       m += 1
    
                feature.Destroy()
            
            print "Number of building POLYGON footprints imported = %s" % n
            print "Number of building MULTIPOLYGON footprints NOT imported = %s" % m
            print "Number of invalid geometries = %s" % i
                        
            if (sql != ""):
                self.cur.execute(sql)
                print "Creating spatial index on buildings...."
                self.cur.execute("CREATE INDEX buildings_idx ON buildings USING GIST (geom);")
                self.conn.commit()

    
            datasource_in.Destroy()
        

        
            return n
        
        else: 
            print "Footprints table already full."
            return n


      
        

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