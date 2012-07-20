def pts_in_polygon(polygon, pointsAll, pointsInside=[], pointsOutside=[]):
# given a set of points and a polygon, return the points that are within the polygon
# see algorithm description at: http://www.ics.uci.edu/~eppstein/161/960307.html#intest

    npoly = len(polygon)
    nInside = 0
    nOutside = 0
    #iterate through the points
    
    ptsNum = 0
    ptsTotal = len(pointsAll)
    for row in pointsAll:
        ptsNum += 1
        i = 0
        x = round(row[0],1)      # pt(x, y) = checking if this point is inside the polygon 
        y = round(row[1],1)
        crossings = 0
        onBoundary = 0
        if(ptsNum % 10000 == 0):
            print "%s of %s" % (str(ptsNum), str(ptsTotal))
        #iterate through the polygon points
        while i < (npoly):
            x1 = round(polygon[i][0],1)  #pt(x1,y1)  = point 1 in the polygon segment
            y1 = round(polygon[i][1],1)  
            
            if (i == (npoly - 1)):
                x2 = round(polygon[0][0],1)  #pt(x2,y2) = point 2 in the polygon segment
                y2 = round(polygon[0][1],1)
            else:
                x2 = round(polygon[i+1][0],1)  #pt(x2,y2) = point 2 in the polygon segment
                y2 = round(polygon[i+1][1],1)

            #check to see if pt(x,y) intersects with polygon segment
            if ( ((x > x1) and (x < x2)) or ((x < x1) and (x > x2)) ): 
                #find the intersection point (i.e. where point ray intersects on the y axis)
                t = float((x - x2) / (x1 - x2))
                yc = round((t*y1 + (1-t)*y2),1)     #round to the nearest centimeter, assuming input data is in meters)
                                
                if (y == yc):
                    onBoundary = 1
                    break
                elif (y > yc):
                    crossings += 1
                                #need to check for the case where the pt(x,y) grazes a vertex
            elif (x == x1) and (y >= y1):     #if the pt(x,y) goes through the first point in the segment
                
                #determine previous poly segment index
                if (i == 0): i_prev = npoly - 1
                else: i_prev = i - 1
               
                if (y == y1):       # check if point lies on vertex of segment
                    onBoundary = 1
                    break
                else:
                    #perturb Left and check crossings of poly segment P(i+1) and P(i-1)
                    x0 = polygon[i_prev][0]
                    y0 = polygon[i_prev][1]

                    if ((x == x0) and (y < y0)):       #pt(x,y) lies on a vertical line
                        onBoundary = 1
                        break
                    elif ( ((x0 < x1) and (x2 < x1)) or ((x0 > x1) and x2 > x1)):
                        z = 0
                        #grazes a point, no crossing
                    elif ((x2 > x1) or (x0 > x1)):
                        #xp = (x2 - x1)/2 + x1
                        crossings += 1
                        
            i += 1      #iterate the polygon vertex index
        
    
        if ((onBoundary != 1) and ((crossings % 2) == 1)):
            #if: the point is not on the boundary and the number of crossings is odd
            #then: point is inside the polygon - append to the pointsInside array
            pointsInside.append(row)
            nInside += 1


        else:
            #point is outside the polygon - append to the pointsOutside array
            pointsOutside.append(row)
            nOutside += 1                                
            
def asciiGrid_to_array(file_in,bounds,points=[], reduction=10):
# Input: file_in - an ASCII grid file (i.e. ESRI)
#        points -  an empty array which will hold the points from the ASC grid file
#        reduction - the amount to reduce the input grid (i.e reduction=10 only writes every 10 pts)

    ascfile = open(file_in, "r").readlines() 
    ncols = int(ascfile[0].split()[1])
    nrows = int(ascfile[1].split()[1])
    xllcenter = float(ascfile[2].split()[1])
    yllcenter = float(ascfile[3].split()[1])
    cell_size = float(ascfile[4].split()[1])
    nodata_value = float(ascfile[5].split()[1])

    i = 0
    n = 0
    j = nrows
    
    while n < nrows:
        i = 0
        row = ascfile[n+6].split()
        while i < ncols:
            #start in the top left corner (i = 0, j = nrows)
            x = float(xllcenter + i*cell_size)
            y = float(yllcenter + (j-1)*cell_size)
            z = round(float(row[i]),4)
            #only append points that are within the bounding box
            if (((x > bounds[0]) and (x < bounds[1])) and ((y > bounds[2]) and (y < bounds[3]))):
                if ((i % int(reduction) == 0) and (j % int(reduction) == 0)):
                    points.append([x,y,z])
            i = i + 1
        j = j - 1
        n += 1

def NODfile_to_polygonlist(file_in, polygonlist=[]):

    #TODO - check if the .nod file is valid (see header below)
    #    #NOD
    #    0.0 0.0 1.0 1.0 0
    #     1640 
    #     2             0 
    #     333 

    polygon = []
    plyfile = open(file_in, "r").readlines() 
    npoly = int(plyfile[3].split()[0])    
    nrow = 4
    
    n = 0   #polygon index
    while (n < npoly):
        numpts = int(plyfile[nrow].split()[0])  
        i = 0
        nrow += 1
        polygon = []
        while (i < numpts):
            line = plyfile[nrow].split()
            vertex = [float(line[0]),float(line[1]),float(line[2])] 
            polygon.append(vertex)
            i += 1
            nrow += 1  
        polygonlist.append(polygon)
        n += 1
    
    


def NODfile_write(file_out, pointslist, polygonlist):  

# Input: point - points to write to NOD file, file_in
#        polygonlist - List of polygons to write to NOD file
#                      Assume the first polygon contains all the other polygons
#        pointslist - a set of points arrays

#TODO:  Check polygonlist to make sure that the first polygon contains all further ones

    
    npolys = len(polygonlist)
    outfile = open(file_out, "w")    
    #line = "#NOD \n0.00000       0.00000       1.00000       1.00000    0\n"
    #outfile.write(line)
    line = "VARIABLES = \"X\", \"Y\", \"Z\"\n"
    outfile.write(line)
    numPoints = 0
    for points in pointslist:
        numPoints += len(points)
    
    numPolyPts = 0
    for polygon in polygonlist:
        numPolyPts += len(polygon)
    
    #line = "%s \n%s       0\n" % (str(numPolyPts+numPoints), str(npolys))
    #outfile.write(line)
    
    #write the polygons (i.e. boundaries and islands
#    for polygon in polygonlist:
##        line = "%s\n" % (str(len(polygon)))
##        outfile.write(line)
#        for row in polygon:
#            line = "%s %s %s\n" % (str(row[0]).ljust(14), str(row[1]).ljust(14), str(row[2]))
#            outfile.write(line)
    
    #write the points
#    line = "%s\n" % (str(numPoints))
#    outfile.write(line)            
    
    for points in pointslist:
        for row in points:
            line = "%s %s %s\n" % (str(row[0]).ljust(14), str(row[1]).ljust(14), str(row[2]))
            outfile.write(line)

    outfile.close()     
    
    
    
def get_polygon_bounds(polygon,bounds):  
    #find the bounds (xmin,xmax,ymin,ymax) of a polygon
    
    
    xmin = polygon[0][0]
    xmax = polygon[0][0]
    ymin = polygon[0][1]
    ymax = polygon[0][1]
    
    for vertex in polygon:
        if (vertex[0] < xmin): xmin = vertex[0]
        if (vertex[0] > xmax): xmax = vertex[0]
        if (vertex[1] < ymin): ymin = vertex[1]
        if (vertex[1] > ymax): ymax = vertex[1]
        
    bounds.append(xmin)
    bounds.append(xmax)
    bounds.append(ymin)
    bounds.append(ymax)
    


def createFaultPatch(fileout, epsg_in = 4326, epsg_out = 32758):
    
    
    from osgeo import ogr
    from osgeo import osr
    import math
    

    originLL = []                   # LAT/LONG coordinates to the lower left corner of the fault patch
    originUTM = []                  # UTM coordinates to the lower left corner of the fault patch
    originLL = [float(165.59), float(-46.31)]
    
    epsgLL = 4326
    epsgUTM = 32758
    faultWidth = 85.0
    faultLength = 125.0
    segmentWidthRupture = 5000.0
    segmentLength = 5000.0
    
    segmentWidthSurface = (54000.0/80000.0)*5000.0
    
    strike = -math.pi/7.4
    strikeDegrees = (180/math.pi)*strike
    
    outfile = open(fileout, "w")    
    
    # set the spatial reference - of the input data
    coordsLL = osr.SpatialReference()
    coordsLL.ImportFromEPSG(epsgLL)
    # set the spatial reference - of the output data    
    coordsUTM = osr.SpatialReference() 
    coordsUTM.ImportFromEPSG(epsgUTM)
    point = ogr.Geometry(type=ogr.wkbPoint) 
       
       
    # convert LL fault-slip origin point to UTM
    point.SetPoint(0, originLL[0], originLL[1], float(0.0))
    point.AssignSpatialReference(coordsLL)
    point.TransformTo(coordsUTM)
    originUTM = [float(point.GetX()),float(point.GetY())]
    
    
    # Create new shapefile
    driver = ogr.GetDriverByName('KML')
    shapeData = driver.CreateDataSource('./test/test2.kml')
    #create layer
    layer =  shapeData.CreateLayer("FaultPatch", coordsLL, ogr.wkbPoint)
    
    numSegX = 17
    numSegY = 25
        
    R = [[math.cos(strike), -math.sin(strike)],[math.sin(strike),math.cos(strike)]]
    
    currentSegX = 0.0
    dip = 0.0
    dipDegrees = 0.0
    
    line = "X         Y         Length    Width     Depth     Dip       Strike    C1        C2        Rake\n"
    outfile.write(line)
    
    i = 0
    while i < numSegX:
        
        print "CurrentSegX = %s" % (str(currentSegX))   
        
        dip = math.atan(0.000006*currentSegX - 0.3231)
        dipDegrees = (180/math.pi)*dip
        depth = 0.000003*currentSegX*currentSegX + 0.3231*currentSegX + 5678.9

        j = 0        
        while j < numSegY:
            x0 = currentSegX
            y0 = j*segmentLength
            xUTM = x0*R[0][0] + y0*R[0][1] + originUTM[0]
            yUTM = x0*R[1][0] + y0*R[1][1] + originUTM[1]
            point = ogr.Geometry(type=ogr.wkbPoint) 
            point.SetPoint(0, xUTM, yUTM)
            point.AssignSpatialReference(coordsUTM)
            point.TransformTo(coordsLL)
            xLL = float(point.GetX())
            yLL = float(point.GetY()) 
            #faultPatchUTM[i][j] = [xUTM,yUTM]
            #faultPatchLL[i][j] = [xLL,yLL]

            feature = ogr.Feature(layer.GetLayerDefn())
            feature.SetGeometry(point)
            layer.CreateFeature(feature)
#            feature.Destroy()   


            line = "%s  %s  %s  %s  %s  %s  %s  2.e10     2.e10   %s  90.0\n" % (str(xLL).ljust(8), \
                                         str(yLL).ljust(8), str(segmentLength).ljust(8), \
                                         str(segmentWidthRupture).ljust(8), str(depth).ljust(8), \
                                         str(dipDegrees).ljust(8), str(strikeDegrees).ljust(8), \
                                         str(10.0).ljust(8) )
            outfile.write(line)
            j += 1        
       
        startSeg = segmentWidthRupture*i
        endSeg = segmentWidthRupture*(i+1)
        d1 = 0.000003*startSeg*startSeg + 0.3231*startSeg + 5678.9
        d2 = 0.000003*endSeg*endSeg + 0.3231*endSeg + 5678.9
        segmentWidthSurface = math.sqrt(segmentWidthRupture*segmentWidthRupture - (d2 - d1)*(d2 - d1))
        currentSegX = currentSegX + segmentWidthSurface   
        print "    segmentWidthSurface = %s\n" % (str(segmentWidthSurface))
     
        i += 1
        
    shapeData.Destroy()
    feature.Destroy()   
        
    outfile.close()
    





def convert_coords(nghfile_in, nghfile_out, epsg_in = 28355, epsg_out = 4326):
    
    
    from osgeo import ogr
    from osgeo import osr

    

    file_in = open(nghfile_in, "r").readlines() 
    outfile = open(nghfile_out, "w")    
    numNodes = int(file_in[2].split()[0])
    numNeighs = int(file_in[3].split()[0])

    print "Number of Nodes = %s" % numNodes

    print "Number of Neighs = %s" % numNeighs

    # set the spatial reference - of the input data
    coords_in = osr.SpatialReference()
    coords_in.ImportFromEPSG(epsg_in)
    # set the spatial reference - of the output data    
    coords_out = osr.SpatialReference() 
    coords_out.ImportFromEPSG(epsg_out)
    point = ogr.Geometry(type=ogr.wkbPoint25D)
    
    #write header of new NGH file
    outfile.write("#NGH\n")
    outfile.write(file_in[1])
    outfile.write(file_in[2])
    outfile.write(file_in[3])

    #convert each point from zone 55UTM to Latlong
    n = 0
    while n < numNodes:
        nghLine = file_in[n+4].split()
        x = float(nghLine[1])
        y = float(nghLine[2])
        z = float(nghLine[4])
        
        #Create OGR point object and transform the coordinates
        point.SetPoint(0, x, y, int(z))
        point.AssignSpatialReference(coords_in)
        point.TransformTo(coords_out)
        x = float(point.GetX())
        y = float(point.GetY())

        line = "%s %s %s %s %s   " % ((nghLine[0]).ljust(9), str(x).ljust(14), str(y).ljust(14), str(nghLine[3]), str(z).ljust(10))
        i = 0
        while i < numNeighs:
            line = line + str(nghLine[i+5]).ljust(8)
            i+=1
            
        line = line + "\n"
        outfile.write(line)
        n+=1
        
    outfile.close()    




def write_PLYfile(polygonlist,file_out):
    
    
    outfile = open(file_out, "w")    

    n = 1
    for polygon in polygonlist:
        npts = len(polygon)
        line = '%s   %s\n' % (str(npts), str(n))
        outfile.write(line)
        for row in polygon:
            line = "%s  %s \n" % (str(row[0]).ljust(14), str(row[1]).ljust(14))
            outfile.write(line)
            
        n += 1
        
        
        
        
        
        
        
def puysegurWorstCase():
    
    
    from osgeo import ogr
    from osgeo import osr
    import math
    from numpy import matrix

    #proj4 definition of local coordinate system
    #origin is in the middle of the bounding box
    proj4 = '+proj=tmerc +lat_0=-47.0 +lon_0=166.0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +k=1'   
    faultPatches = []
    epsgLL = 4326
    outfile = open('./test/puysegurWorstCase.dat', "w")    


    Puy1 = [-128350.0,-196020.0,-56440.0,115380.0,0.0,8000.0]
    Puy2 = [-83760.0,-151040.0,-17250.0,137250.0,8000.0,14800.0]
    Puy3 = [-36500.0,-72250.0,19500.0,168000.0,14800.0,20200.0]
    Puy4 = [-650.0,-17600.0,49000.0,194000.0,20200.0,25800.0]
    Puy5 = [41000.0,57750.0,72000.0,192000.0,25800.0,30000.0]
    faultPatches.append(Puy1)
    faultPatches.append(Puy2)
    faultPatches.append(Puy3)
    faultPatches.append(Puy4)
    faultPatches.append(Puy5)

    # set the spatial reference - of the input data
    coordsLL = osr.SpatialReference()
    coordsLL.ImportFromEPSG(epsgLL)
    # set the spatial reference of the fault patch (as defined by Hayes and Furlong 2010)   
    coordsPatch = osr.SpatialReference()
    coordsPatch.ImportFromProj4(proj4)
    
    driver = ogr.GetDriverByName('KML')
    kmlData = driver.CreateDataSource("./test/PuysegurWorstCasePatch2.kml")
    layer =  kmlData.CreateLayer("HayesPatch", coordsLL, ogr.wkbLineString)
    line = ogr.Geometry(type=ogr.wkbLineString)

    dip = 13.5
    dipRadians = (math.pi/180)*dip
    strike = 13.5
    strikeRadians = (math.pi/180)*strike
    
    for patch in faultPatches:
        widthSurface = (patch[5]-patch[4])/math.tan(dipRadians)
        widthRupture = (patch[5]-patch[4])/math.sin(dipRadians)
        xp = widthSurface*math.cos(-strikeRadians)
        yp = widthSurface*math.sin(-strikeRadians)
        lengthRupture = math.sqrt((patch[2]-patch[0])*(patch[2]-patch[0]) + (patch[3]-patch[1])*(patch[3]-patch[1]) )
        xc = (patch[0]+patch[2])/2 + xp/2
        yc = (patch[1]+patch[3])/2 + yp/2   
        dc = (patch[5]+patch[4])/2
        
        line = ogr.Geometry(type=ogr.wkbLineString)
        line.AddPoint(patch[0],patch[1])
        line.AddPoint(patch[0]+xp,patch[1]+yp)
        line.AddPoint(patch[2]+xp,patch[3]+yp)
        line.AddPoint(patch[2],patch[3])
        line.AddPoint(patch[0],patch[1])
  
        point = ogr.Geometry(type=ogr.wkbPoint) 
        point.SetPoint(0,xc,yc)
        point.AssignSpatialReference(coordsPatch)
        point.TransformTo(coordsLL)

        line.AssignSpatialReference(coordsPatch)
        line.TransformTo(coordsLL)
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(line)
        layer.CreateFeature(feature)
        feature.SetGeometry(point)
        layer.CreateFeature(feature)

        xc = float(point.GetX())
        yc = float(point.GetY())
        
        #line = "%s   %s   %s   %s   %s  %s  19.0    2.e10  2.e10   4.0   144.0\n" % (str(xc-15.0), str(yc), str(lengthRupture), \
        #                                       str(widthRupture),str(dc),str(dip))
        
        outfile.write(line)
    kmlData.Destroy()
    feature.Destroy()
    outfile.close()
    
    
def polygonList2KML(polygonlist=[]):

    from osgeo import ogr
    from osgeo import osr
    epsgLL = 4326
    epsgUTM = 32755
    # set the spatial reference - of the input data
    coordsLL = osr.SpatialReference()
    coordsLL.ImportFromEPSG(epsgLL)
    # set the spatial reference of the fault patch (as defined by Hayes and Furlong 2010)   
    coordsUTM = osr.SpatialReference()
    coordsUTM.ImportFromEPSG(epsgUTM)
    
    driver = ogr.GetDriverByName('KML')
    kmlData = driver.CreateDataSource("./test/NZTAS_splitsADDED_BC7_Redepped.kml")
    layer =  kmlData.CreateLayer("RicomBoundary", coordsLL, ogr.wkbLineString)

    for poly in polygonlist:
        line = ogr.Geometry(type=ogr.wkbLineString)
        for point in poly:
            line.AddPoint(point[0],point[1])

        line.AssignSpatialReference(coordsUTM)
        line.TransformTo(coordsLL)
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(line)
        layer.CreateFeature(feature)

    kmlData.Destroy()
    feature.Destroy()
    
    

def readSTSMaxFile(stsfile_in, stsfile_out, epsg_in = 32755, epsg_out = 4326):


    import csv
    from osgeo import ogr
    from osgeo import osr
    
    proj4 = '+proj=ortho +lat_0=-46.0 +lon_0=155.5 +x_0=0 +y_0=0'   

  

    infile = csv.reader(open(stsfile_in, 'r'), delimiter=',', quotechar=None)
    header = dict(((str, i) for i, str in enumerate(infile.next())))
    outfile = open(stsfile_out, "w") 
    
    # set the spatial reference - of the input data
    coords_in = osr.SpatialReference()
    coords_in.ImportFromEPSG(epsg_in)
    # set the spatial reference - of the output data    
    coords_out = osr.SpatialReference() 
    coords_out.ImportFromEPSG(epsg_out)    
    coordsRicom = osr.SpatialReference()
    coordsRicom.ImportFromProj4(proj4)       



    for row in infile:
        point = ogr.Geometry(type=ogr.wkbPoint)
        point.SetPoint(0, float(row[1]), float(row[2]))
        point.AssignSpatialReference(coords_in)
        point.TransformTo(coordsRicom)
        x = float(point.GetX())
        y = float(point.GetY())
        #z = int(point.GetZ())
        #nodes_out.append([x, y, float(row[2])])
        line = "%s  %s  %s %s\n" % (str(x).ljust(14), str(y).ljust(14), str(row[3]), str(row[0]))
        outfile.write(line)    
    
    
    outfile.close()


def ricomTs2anugaIC():
    #input: Ricom time series output (ETA, velocity) at a series of points (that define an ANUGA driving boundary)
    #output: ANUGA sts driving boundary file for each of the points

    import csv  
    tsdata = open("../input/anuga/tsdata.dat", "r").readlines() 
    tsUdata = open("../input/anuga/tsUdata.dat", "r").readlines()
    


    #URS order file is the indices and locations of the driving boundary (according to ANUGA)
    #index,longitude,latitude

    ursorder = csv.reader(open("../input/anuga/urs_order.csv", 'r'), delimiter=',', quotechar=None)
    ursorder.next()
    
    ricomProbedPts = csv.reader(open("../input/anuga/anugaBoundaryFileElementsProbeDepth.csv", 'r'), delimiter=',', quotechar=None)
    ricomProbedPts.next()

    ntsPts = int(tsdata[0].split()[0])
    ntsSteps = int(tsdata[0].split()[1])

    ntsUPts = int(tsUdata[0].split()[0])
    ntsUSteps = int(tsUdata[0].split()[1])
    
    eta = []
    velocity = []
    
    if (ntsUPts == ntsPts and ntsUSteps == ntsSteps):
        timeStep = 0    # timeStep number   
        i = 1           # row in tsfile
        j = 1           # row in tsUfile
        
        while (timeStep < ntsSteps):
            n = 0           # output point number    
            line = []       
            while n != ntsPts+1:
                row = tsdata[i].split()
                line.extend(row)
                n = n + len(row)
                i += 1
            
            eta.append(line)
            n = 0
            line = []
            while n !=  (ntsUPts*2)+1:
                row = tsUdata[j].split()
                line.extend(row)
                n = n + len(row)        
                j += 1
            
            velocity.append(line)
            timeStep += 1
    
    n = 0   #boudary point number
    #make a STS boundary file for each point in input time series data
    for index in ursorder:
        outfile = open('../output/anuga/ricomDrivingBoundaryforANUGA/sts_gauge_'+str(index[0])+'.csv', "w") 
        i = 0
        outfile.write('time, stage, xmomentum, ymomentum\n')
        
        point = ricomProbedPts.next()
        elevation =  float(point[2])
        
        print "Point: x = %s, y = %s, z = %s, elNum = %s\n" % (point[0],point[1],point[2],point[3]) 
        #in ANUGA xmomentum and ymomentum are "defined" as height X x_velocity and height X y_velocity
        
        while i < ntsSteps:
            if float(eta[i][0]) >= 5857.5:
                line  = "%s, %s, %s, %s\n" % (eta[i][0], eta[i][n+1], str(float(velocity[i][n*2+1])*( float(eta[i][n+1]) - elevation)),str( float(velocity[i][n*2+2]) * (float(eta[i][n+1]) - elevation)))
                outfile.write(line)
            i += 1
        
        n += 1
        outfile.close() 


def ply2kml(csvfile_in, kmlfile_out, epsg_in = 32755, epsg_out = 4326):


    import csv
    from osgeo import ogr
    from osgeo import osr
    
    infile = csv.reader(open(csvfile_in, 'r'), delimiter=',', quotechar=None)
    header = dict(((str, i) for i, str in enumerate(infile.next())))    
    
    epsgLL = 4326
    epsgUTM = 32755
    # set the spatial reference - of the input data
    coordsLL = osr.SpatialReference()
    coordsLL.ImportFromEPSG(epsgLL)
    # set the spatial reference of the fault patch (as defined by Hayes and Furlong 2010)   
    coordsUTM = osr.SpatialReference()
    coordsUTM.ImportFromEPSG(epsgUTM)
    
    driver = ogr.GetDriverByName('KML')
    kmlData = driver.CreateDataSource(kmlfile_out)
    layer =  kmlData.CreateLayer("PriorityAreasKB", coordsLL, ogr.wkbLineString) 
    
    plyValue = 0
    polygon_list = []
    polygon = []
    
    for point in infile:
        if point[2] != plyValue and plyValue != 0:
            polygon_list.append(polygon)
            polygon = []    
        plyValue = point[2]
        polygon.append(point)
        
    for polygon in polygon_list:
        
        line = ogr.Geometry(type=ogr.wkbLineString)
        
        for point in polygon:
            line.AddPoint(float(point[0]),float(point[1]))
        
        line.AssignSpatialReference(coordsUTM)
        line.TransformTo(coordsLL)
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(line)
        layer.CreateFeature(feature)

    feature.Destroy()
    kmlData.Destroy()

def ply2kml2(csvfile_in, kmlfile_out, epsg_in = 32755, epsg_out = 4326):


    import csv
    from osgeo import ogr
    from osgeo import osr
    
    infile = csv.reader(open(csvfile_in, 'r'), delimiter=',', quotechar=None)    
    #file_in = open(csvfile_in, "r").readlines() 

    
    epsgLL = 4326
    epsgUTM = 32755
    # set the spatial reference - of the input data
    coordsLL = osr.SpatialReference()
    coordsLL.ImportFromEPSG(epsgLL)
    # set the spatial reference of the fault patch (as defined by Hayes and Furlong 2010)   
    coordsUTM = osr.SpatialReference()
    coordsUTM.ImportFromEPSG(epsgUTM)
    
    driver = ogr.GetDriverByName('KML')
    kmlData = driver.CreateDataSource(kmlfile_out)
    layer =  kmlData.CreateLayer("PolyAREA", coordsLL, ogr.wkbLineString) 
    
    polygon = []
    line = ogr.Geometry(type=ogr.wkbLineString)

    for point in infile:  
        polygon.append(point)
        line.AddPoint(float(point[0]),float(point[1]))
        
    line.AssignSpatialReference(coordsUTM)
    line.TransformTo(coordsLL)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(line)
    layer.CreateFeature(feature)

    feature.Destroy()
    kmlData.Destroy()


def nod2csv(csvfile_in, kmlfile_out, csvfile_out, epsg_in = 32755, epsg_out = 4326):


    from osgeo import ogr
    from osgeo import osr

    infile = open(csvfile_in, "r").readlines() 
    outfile = open(csvfile_out, "w") 

    
    epsgLL = 4326
    epsgUTM = 32755
    # set the spatial reference - of the input data
    coordsLL = osr.SpatialReference()
    coordsLL.ImportFromEPSG(epsgLL)
    # set the spatial reference of the fault patch (as defined by Hayes and Furlong 2010)   
    coordsUTM = osr.SpatialReference()
    coordsUTM.ImportFromEPSG(epsgUTM)
    
    driver = ogr.GetDriverByName('KML')
    kmlData = driver.CreateDataSource(kmlfile_out)
    layer =  kmlData.CreateLayer("PolyAREA", coordsLL, ogr.wkbLineString) 
    
    polygon = []
    line = ogr.Geometry(type=ogr.wkbLineString)

    for row in infile:  
        point = row.split()
        line.AddPoint(float(point[0]),float(point[1]))
        linecsv = "%s,%s\n" % (str(point[0]),str(point[1]))
        outfile.write(linecsv)
        
    line.AssignSpatialReference(coordsUTM)
    line.TransformTo(coordsLL)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(line)
    layer.CreateFeature(feature)

    feature.Destroy()
    kmlData.Destroy()
    outfile.close()
    
    
def LL2LocalRicom(pointsLL,lat0=0,long0=0,latoff=0,longoff=0):
    
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
        
    
def convertXY (pointsIN,epsgIN,epsgOUT):
   
    from osgeo import ogr
    from osgeo import osr   
       
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
    
    
def fixSHPfile(input_shp, output_shp):
    
    from osgeo import ogr
    from osgeo import osr  
    import sys

    datasource_in = ogr.Open(input_shp)    
    if datasource_in is None:
        print "Open Failed.\n"
        sys.exit(1)
    
    #print datasource_in.GetLayerCount()
    #layer_in = datasource_in.GetLayerByName("footprints")
    layer_in = datasource_in.GetLayer(0)
    layer_in.ResetReading()
    geomType = layer_in.GetGeomType()
    
    if (geomType == ogr.wkbPolygon or geomType == ogr.wkbMultiPolygon):
        #input geometry is a polygon
        # Create new shapefile    
        driver = ogr.GetDriverByName('ESRI Shapefile')
        datasource_out = driver.CreateDataSource(output_shp)
        layer_out =  datasource_out.CreateLayer("footprints", layer_in.GetSpatialRef(), ogr.wkbPolygon)
    
        n = 0
        for feature in layer_in:
            geom = feature.GetGeometryRef()
            if geom is not None and geom.GetGeometryType() ==  ogr.wkbPolygon:
                
                count = geom.GetGeometryCount()
                n += 1
                if geom.IsValid():
                    layer_out.CreateFeature(feature)
                    
        datasource_out.Destroy()
    else:
        print "ERROR: Incorrect input geometry.  Building footprints must be a POLYGON.\n"
    
    
    feature.Destroy()
    datasource_in.Destroy()

            
            
def addPolygonListPG(polyList,dbname="kingstonbeach",user="tbone",tablename="polylist",epsg = "28355"):
    
    import psycopg2
    import sys
    
    #open the conection to the requested database
    line = "dbname='" + dbname + "' user='" + user +"'"
     
    try:
        conn = psycopg2.connect(line);
    except:
        print "Unable to connect to the database" 
        sys.exit()   

    cur = conn.cursor()


    #DROP table=tablename already exists
    cur.execute("DROP TABLE IF EXISTS %s;" % (tablename))
    #CREATE table=tablename
    cur.execute("CREATE TABLE %s (id int4)"  % (tablename))  
    cur.execute("SELECT AddGeometryColumn( '%s', 'polygon',%s,'POLYGON', 2 );" % (tablename,epsg))
  
    n = 1

    for poly in polyList:
        polySQL = "ST_GeomFromText('POLYGON((%s %s" % (poly[0][0],poly[0][1])
        poly.pop(0)
        for pt in poly:
            polySQL = polySQL + ", %s %s" % (pt[0],pt[1])

        polySQL = polySQL + "))',%s));" % (epsg)

        line = "INSERT INTO %s (id,polygon) VALUES (%s, " % (tablename,n)
        line = line + polySQL
        cur.execute(line)  
        n += 1  
        
        
    conn.commit()
    cur.close()
    conn.close()
        
    
def pg_geomList2Kml(geomList=[],kml_file="PostGIS.kml"):
#Writes a list of arbitrary PostGIS geometries to a KML file
#---------------------------------------------------------------------------------
#INPUT:     A list of PostGIS geometry entities as strings (i.e. POINT, POLYGON, etc.).
#OUTPUT:    Output KML file
#---------------------------------------------------------------------------------

    from osgeo import ogr
    from osgeo import osr 

    if(geomList != []):    
        
        polyList = []
        pointsList = []
        
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
                pt = geom[1][i+6:j]
                pt = pt.replace(',',' ').split()
                pointsList.append([float(pt[0]),float(pt[1])])
            else:
                print "Unhandled SQL geometry..."
    
        # set the spatial reference - of the polygon data
        coordsLL = osr.SpatialReference()
        coordsLL.ImportFromEPSG(4326)
        coordsUTM = osr.SpatialReference()
        coordsUTM.ImportFromEPSG(int(28355))
    
        #write the polygon data to KML
        driver = ogr.GetDriverByName('KML')
        kmlData = driver.CreateDataSource(kml_file)
        layerPolygons =  kmlData.CreateLayer("polygons", coordsLL, ogr.wkbLineString)
        layerPoints =  kmlData.CreateLayer("points", coordsLL, ogr.wkbPoint)
    
        featurePolygons = ogr.Feature(layerPolygons.GetLayerDefn())
        featurePoints = ogr.Feature(layerPoints.GetLayerDefn())
        
        #write the polygon data to the KML file
        for poly in polyList:
            line = ogr.Geometry(type=ogr.wkbLineString)
            for pt in poly:
                line.AddPoint(pt[0],pt[1])
                  
            line.AssignSpatialReference(coordsUTM)
            line.TransformTo(coordsLL)
            featurePolygons.SetGeometry(line)
            layerPolygons.CreateFeature(featurePolygons)
    
        #write the point data to the KML file
        for pt in pointsList:
            point = ogr.Geometry(type=ogr.wkbPoint) 
            point.SetPoint(0, pt[0], pt[1])
            point.AssignSpatialReference(coordsUTM)
            point.TransformTo(coordsLL)
            featurePoints.SetGeometry(point)
            layerPoints.CreateFeature(featurePoints)
        
        kmlData.Destroy()
        featurePolygons.Destroy()
        featurePoints.Destroy()
        
    else:
        print "No point in PostGIS Geometry List!!"
    

def pg_to_ricom(node_table="", tri_table="", ngh_file ="", tri_file ="",dbname="kingstonbeach",user="tbone",epsgOUT =4326):


 #write the ngh file header
    ngh_outfile.write("%s\n%s\n" % (str(len(vertices)),str(number_of_neighbours)))
    ngh_outfile.write("0.0     0.0     0.0     0.0\n")
    i = 1
    #write vertices
    for vert in vertices:
        line = "%s   %s    %s    %s     %s   " % (str(i), str(vert[0]), str(vert[1]), str(vert[2]), str(vert[3]))
        n = 0
        neighs = len(vert) - 4
        #construct the neighbours list
        while n < neighs:                   
            nline = "%s   " % (str(vert[4+n]))
            line = line + nline
            n += 1
        while n < number_of_neighbours:
            line = line + "0  "
            n += 1
        
        line = line + "\n"
        ngh_outfile.write(line) 
        i += 1
            
    ngh_outfile.close()
    
    for el in elements:
        el_outfile.write("%s  %s  %s  0   1\n" % (str(el[0]).ljust(8),str(el[1]).ljust(8),str(el[2]).ljust(8)))    


 

def shapefile_to_point(shp_filename,out_filename,epsgOUT):
    """
    This function imports a shapefile into a the postgis database
    
    IN: filename = shapefile filename - assume shapefile consists of LINESTRING only
        tablename = name of the table to be created in the database
        boundaryID - the id of the boundary polygon inside which the shapefile features are to be imported (i.e. the study area)
        
    OUT: function returns # of geometries imported if successful or 0 if unsuccessful
    
    """
    
    from osgeo import ogr
    from osgeo import osr 
    
    datasource_in = ogr.Open(shp_filename)
    outfile = open(out_filename, "w") 
    outfile.write("x y z\n")
    
    layer_in = datasource_in.GetLayer(0)
    layer_in.ResetReading()
    
    srs_in = layer_in.GetSpatialRef()       # the footprints spatial reference
    srs = osr.SpatialReference()            # the spatial reference of the PostGIS database (i.e. self.dbname)
    srs.ImportFromEPSG(epsgOUT)

    points = []
    for feature in layer_in:
        geom = feature.GetGeometryRef()
        geom.TransformTo(srs) 
        i = 0
        ptCount = geom.GetPointCount()
        depth = feature.GetFieldAsDouble('Depth')
        while i < ptCount:
            p = geom.GetPoint(i)
            outfile.write("%s %s %s\n" %(p[0], p[1], depth))
            i += 1        
        feature.Destroy()
        

    datasource_in.Destroy()

    outfile.close()


def ogr_LineString2Polygon(ogr_data_in,filename_out, epsgOUT):

	'''
	Given and ogr data source of closed linestring data (i.e. first pt == last pt)
	
	Convert to a polygon feature (written to filename_out)
	

	'''
	from osgeo import ogr
	from osgeo import osr 
	
	datasource_in = ogr.Open(ogr_data_in)
	layer_in = datasource_in.GetLayer(0)
	layer_in.ResetReading()
	
	srs_in = layer_in.GetSpatialRef()      
	srs_out = osr.SpatialReference()            
	srs_out.ImportFromEPSG(epsgOUT)
	
	#write the polygon data to KML
	driver = ogr.GetDriverByName('ESRI Shapefile')
	outfile = driver.CreateDataSource(filename_out)
	layerPolygons =  outfile.CreateLayer("Polygons", srs_out, ogr.wkbPolygon)	
	featurePolygons = ogr.Feature(layerPolygons.GetLayerDefn())
	
	#NOTE: To create a polygon you first need to create a ring
	#       then add it to the polygon
	for feature in layer_in:
		geom = feature.GetGeometryRef()
		geom.TransformTo(srs_out) 
		ptCount = geom.GetPointCount()
		firstPt = geom.GetPoint(0)
		lastPt = geom.GetPoint(ptCount-1)
		polygon = ogr.Geometry(type=ogr.wkbPolygon)			#create polygon
		ring = ogr.Geometry(type=ogr.wkbLinearRing)			#create ring
		
		#Only convert the feature to polygon if is a closed linestring
		if firstPt[0] == lastPt[0] and firstPt[1] == lastPt[1]:
			print "Valid LINESTRING"      
			i = 0
			p1 = geom.GetPoint(0)
			p2 = geom.GetPoint(1)
			p3 = geom.GetPoint(2)
			p4 = geom.GetPoint(3)
			p5 = geom.GetPoint(4)
			ring.AddPoint(p1[0],p1[1],4.5)
			ring.AddPoint(p2[0],p2[1],0)
			ring.AddPoint(p3[0],p3[1],0)
			ring.AddPoint(p4[0],p4[1],4.5)
			ring.AddPoint(p5[0],p5[1],4.5)
			'''
			while i < ptCount-1:
				p = geom.GetPoint(i)
				ring.AddPoint(p[0],p[1])
				i += 1
			'''
			ring.CloseRings()
			polygon.AddGeometry(ring)						#add ring to polygon	
			print polygon
			polygon.AssignSpatialReference(srs_out)
			featurePolygons.SetGeometry(polygon)
			layerPolygons.CreateFeature(featurePolygons) 
		feature.Destroy()
	
	  
	outfile.Destroy()
	featurePolygons.Destroy()

def ogr_LineString2Polygon(ogr_data_in,filename_out, epsgOUT):

	'''
	Given and ogr data source of closed linestring data (i.e. first pt == last pt)
	
	Convert to a polygon feature (written to filename_out)
	

	'''
	from osgeo import ogr
	from osgeo import osr 



	datasource_in = ogr.Open(ogr_data_in)
	layer_in = datasource_in.GetLayer(0)
	layer_in.ResetReading()
	
	srs_in = layer_in.GetSpatialRef()      
	srs_out = osr.SpatialReference()            
	srs_out.ImportFromEPSG(epsgOUT)
	
	#write the polygon data to KML
	driver = ogr.GetDriverByName('ESRI Shapefile')
	outfile = driver.CreateDataSource(filename_out)
	layerPolygons =  outfile.CreateLayer("Polygons", srs_out, ogr.wkbPolygon)	
	featurePolygons = ogr.Feature(layerPolygons.GetLayerDefn())
	
	#NOTE: To create a polygon you first need to create a ring
	#       then add it to the polygon
	for feature in layer_in:
		geom = feature.GetGeometryRef()
		geom.TransformTo(srs_out) 
		ptCount = geom.GetPointCount()
		print ptCount
		firstPt = geom.GetPoint(0)
		lastPt = geom.GetPoint(ptCount-1)
		polygon = ogr.Geometry(type=ogr.wkbPolygon)			#create polygon
		ring = ogr.Geometry(type=ogr.wkbLinearRing)			#create ring
		
		#Only convert the feature to polygon if is a closed linestring
		if firstPt[0] == lastPt[0] and firstPt[1] == lastPt[1]:
			print "Valid LINESTRING"      
			i = 0
			while i < ptCount-1:
				p = geom.GetPoint(i)
				ring.AddPoint(p[0],p[1])
				i += 1
			ring.CloseRings()
			polygon.AddGeometry(ring)						#add ring to polygon	
			print polygon
			polygon.AssignSpatialReference(srs_out)
			featurePolygons.SetGeometry(polygon)
			layerPolygons.CreateFeature(featurePolygons) 
		feature.Destroy()
	
	  
	outfile.Destroy()
	featurePolygons.Destroy()




        