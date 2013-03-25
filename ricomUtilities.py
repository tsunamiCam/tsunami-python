def extract_asciigrid(file_out, nod_file, ascii_filelist, reduction=10):

#Given a nod file defining a boundary and islands, extracts all the points
#contained in the input ascii grid files (i.e. ascii_filelist) that are inside
#the boundary and outside the islands


    from ricomGeospatialUtilities import asciiGrid_to_array
    from ricomGeospatialUtilities import NODfile_to_polygonlist
    from ricomGeospatialUtilities import NODfile_write
    from ricomGeospatialUtilities import pts_in_polygon
    from ricomGeospatialUtilities import get_polygon_bounds
    
    polygonlist = []
    pointsList = []
    bounds = []

    NODfile_to_polygonlist(file_in=nod_file, polygonlist=polygonlist)
    npoly = len(polygonlist)
    
    get_polygon_bounds(polygonlist[0],bounds=bounds)
    print "xmin = %s, xmax = %s, ymin = %s, ymax = %s" % (str(bounds[0]),str(bounds[1]),str(bounds[2]),str(bounds[3]))

#    bounds[0] = 526823.5
#    bounds[1] = 526825.5

    
    for asc_file in ascii_filelist:    

        print "Parsing: %s" % (asc_file)

        points = []
        n = 0
        asciiGrid_to_array(file_in=asc_file,bounds=bounds,points=points, reduction=reduction)
        
        print "Num Pts: %s" % (str(len(points)))

        while (n < npoly):
            pointsInside = []
            pointsOutside = []
            pts_in_polygon(polygon=polygonlist[n], pointsAll=points, pointsInside=pointsInside, pointsOutside=pointsOutside)
            if (n == 0): points = pointsInside
            else: points = pointsOutside 
            n += 1
            
        pointsList.append(points)    
    
    points = []
    NODfile_write(file_out=file_out, pointslist=pointsList, polygonlist=polygonlist)



def extract_fromXYZ(file_out, nod_file, xyz_filelist, reduction=10):

#Given a nod file defining a boundary and islands, extracts all the points
#contained in the input XYZ pts file (i.e. X,Y,Z) that are inside
#the boundary and outside the islands


#Input XYZ file (example):
#x,y,elevation
#623683.860551,5365560.66681,-81.69681
#623933.860551,5365560.66681,-81.99997
#624183.860551,5365560.66681,-82.0
#624433.860551,5365560.66681,-81.99982
#623183.860551,5365310.66681,-81.51879

#To be kind to memory:
#Find the region of the Polygon (nod file) and only extract points inside this


   

    from ricomGeospatialUtilities import asciiGrid_to_array
    from ricomGeospatialUtilities import NODfile_to_polygonlist
    from ricomGeospatialUtilities import NODfile_write
    from ricomGeospatialUtilities import pts_in_polygon
    from ricomGeospatialUtilities import get_polygon_bounds
    import csv

    polygonlist = []
    points = []
    pointsList = []
    bounds = []

    NODfile_to_polygonlist(file_in=nod_file, polygonlist=polygonlist)
    npoly = len(polygonlist)
    get_polygon_bounds(polygonlist[0],bounds=bounds)
    print "xmin = %s, xmax = %s, ymin = %s, ymax = %s\n" % (str(bounds[0]),str(bounds[1]),str(bounds[2]),str(bounds[3]))

    for xyz_file in xyz_filelist:    

        ptfile = csv.reader(open(xyz_file, 'r'), delimiter=',', quotechar=None)
        ptfile.next()
        line = ptfile.next()      
        xminf = float(line[0])
        xmaxf = float(line[0])
        yminf = float(line[1])
        ymaxf = float(line[1])        
        
        i = 0
        for point in ptfile:
            #check if the point is within the polygon bounds
            x = float(point[0])
            y = float(point[1])
            z = float(point[2])
            
#            if (x < xminf): xminf = x
#            if (x > xmaxf): xmaxf = x
#            if (y < yminf): yminf = y
#            if (y > ymaxf): ymaxf = y           
            
            if (((x > bounds[0]) and (x < bounds[1])) and ((y > bounds[2]) and (y < bounds[3]))):
                points.append([x,y,z])
                i += 1
        
        print "# inside bounds = %s\n" % (str(i))
#        print "xminf = %s, xmaxf = %s, yminf = %s, ymaxf = %s\n" % (str(xminf),str(xmaxf),str(yminf),str(ymaxf))

        n = 0        
        while (n < npoly):
            pointsInside = []
            pointsOutside = []
            pts_in_polygon(polygon=polygonlist[n], pointsAll=points, pointsInside=pointsInside, pointsOutside=pointsOutside)
            if (n == 0): points = pointsInside
            else: points = pointsOutside #
            n += 1
 
 
        #pts file from GA is a 10m grid
        #calculate a square boundary around the points
        #start in the bottom left hand corner
        boundary = []
        bounds_points = []
        get_polygon_bounds(points,bounds=bounds_points)
        
        vertexll = [(bounds_points[0]-10), (bounds_points[2]-10)]
        vertexlr = [(bounds_points[1]+10), (bounds_points[2]-10)]
        vertexur = [(bounds_points[1]+10), (bounds_points[3]+10)]
        vertexul = [(bounds_points[0]-10), (bounds_points[3]+10)]
        

        vx = vertexll[0]
        vy = vertexll[1]
        
        while (vx < vertexlr[0]):
            boundary.append([vx,vy])
            vx = vx + 10
        
        vx = vertexlr[0]
        vy = vertexlr[1]
        while (vy < vertexur[1]):
            boundary.append([vx,vy])
            vy = vy + 10

        vx = vertexur[0]
        vy = vertexur[1]

        while (vx > vertexul[0]):
            boundary.append([vx,vy])
            vx = vx - 10
                     
        vx = vertexul[0]
        vy = vertexul[1]
        while (vy > vertexll[1]):
            boundary.append([vx,vy])
            vy = vy - 10
        
        pointsList.append(points)    
        points = []
        
    boundarylist = []
    boundarylist.append(boundary)
    
    NODfile_write(file_out=file_out, pointslist=pointsList, polygonlist=boundarylist)

    
    
    
def extract_fromGEBCO(tecfile_out, nod_file, gebco_file):
    
  
    from ricomGeospatialUtilities import NODfile_to_polygonlist
    from ricomGeospatialUtilities import get_polygon_bounds  
    from ricomGeospatialUtilities import NODfile_write
    
    import csv
    
    #Input GEBCO file (example):
    #Longitude Latitude Depth
    #139.595833 -34.004167 55
    #139.604167 -34.004167 53
    #139.612500 -34.004167 49  
    
    ptfile = csv.reader(open(gebco_file, 'r'), delimiter=' ', quotechar=None)
    ptfile.next()
    line = ptfile.next()      

    polygonlist = []
    bounds = []
    points = []
    pointslist = []
    NODfile_to_polygonlist(file_in=nod_file, polygonlist=polygonlist)
    npoly = len(polygonlist)
    get_polygon_bounds(polygonlist[0],bounds=bounds)
    
    i = 0
    for point in ptfile:
        #check if the point is within the polygon bounds
        x = float(point[0])
        y = float(point[1])
        z = float(point[2])

        if (((x > bounds[0]) and (x < bounds[1])) and ((y > bounds[2]) and (y < bounds[3]))):
            points.append([x,y,z])
            i += 1

    pointslist.append(points)
    NODfile_write(file_out=tecfile_out, pointslist=pointslist, polygonlist=polygonlist)


def avergeDepthGEBCO(gebco_file):
    
  
    from ricomGeospatialUtilities import NODfile_to_polygonlist
    from ricomGeospatialUtilities import get_polygon_bounds  
    from ricomGeospatialUtilities import NODfile_write
    
    import csv
    
    #Input GEBCO file (example):
    #Longitude Latitude Depth
    #139.595833 -34.004167 55
    #139.604167 -34.004167 53
    #139.612500 -34.004167 49  
    
    ptfile = csv.reader(open(gebco_file, 'r'), delimiter=' ', quotechar=None)
    ptfile.next()
    ptfile.next()      

    points = []
    bounds = [-49.0,164.0,-45.0,168.0]
    
    i = 0
    depthTotal = 0.0
    for point in ptfile:
        #check if the point is within the polygon bounds
        x = float(point[0])
        y = float(point[1])
        z = float(point[2])

        if (((x > bounds[0]) and (x < bounds[1])) and ((y > bounds[2]) and (y < bounds[3]))):
            points.append([x,y,z])
            if(z<0):
                depthTotal = depthTotal + z
                i += 1

    depthAVG = depthTotal/float(i)
    print "AVG Depth = %s\n" % (str(depthAVG))
    
def tec2ngh(tec_in, ngh_in, ngh_out):
    ngh_file = open(ngh_in, "r").readlines() 
    tec_file = open(tec_in, "r").readlines() 
    outfile = open(ngh_out, "w")
    
    i = 0
    j = 0
    k = 0
    num_nodesNGH = 0
    num_nodesTEC = 0    

    for line in tec_file:
        dt = line.find("DT=")
        i = line.find("Nodes=")
        if i > 0:
            j = line.find(",")
            num_nodesTEC = int(line[i+6:j])
        if dt > 0:
            k +=1 
            break
        k+=1
    
    
    num_nodesNGH = int(ngh_file[2].split()[0])
    num_neighbours = int(ngh_file[3].split()[0]) 
   
    outfile.write(ngh_file[0])
    outfile.write(ngh_file[1])
    outfile.write(ngh_file[2])
    outfile.write(ngh_file[3]) 
    
    if num_nodesNGH != num_nodesTEC:
        print "Input Tecplot and NGH file are not compatible. TRY AGAIN"
        return 0
    
    i = 0
    while i < num_nodesTEC:
        l = ngh_file[i+4].split()
        z = tec_file[k+i].split()[2]
        l[4] = z
        newLine = ""
        newLine = "%s %s %s %s %s" % (l[0].ljust(10),l[1].ljust(16),l[2].ljust(16),l[3].ljust(4),l[4].ljust(20))
        j = 0
        while j < num_neighbours:
            newLine = newLine + " %s" % l[5+j].ljust(9)
            j += 1
            
        newLine = newLine + "\n"
        outfile.write(newLine)
        i += 1
            
        
def ngh2tec(ngh_in, tri_in, tec_out):
    
    ngh_file = open(ngh_in, "r").readlines() 
    el_file = open(tri_in, "r").readlines() 
    outfile = open(tec_out, "w")    

    nodes = []                   # List of [x, y] point coordinates  (# rows = number of points) - YES
    nodeCode = []           # List of point attributes (# rows = number of points) - YES
    triangles = []                  # List of triangles [point1 point2 point3]  - point indexes correspond to there position in points[]    - YES
    nodeNumber = []
    numNodes = int(ngh_file[2].split()[0])
    numNeighbours = int(ngh_file[3].split()[0]) 
    neighbours = []
    
    n = 0
    while n < numNodes:
        vertLine = ngh_file[n+4].split()
        node = [float(vertLine[1]), float(vertLine[2]), float(vertLine[4])]
        nodes.append(node)
        nodeNumber.append(int(vertLine[0]))
        
        i = 0
        nghLine = []
        while i < numNeighbours:               #initialise the nodeNeighbours list
            nghLine.append(int(vertLine[5+i]))
            i += 1
            
        neighbours.append(nghLine)
        nodeCode.append(vertLine[3])
        n += 1     

    # Triangles
    #--------------------------------------------------------------------------------
    for line in el_file:
        line_list = line.split() 
        #Ricom is BASE 1 - ANUGA is base 0
        triangles.append([int(line_list[0]), int(line_list[1]), int(line_list[2])]) 
    
    numTriangles = len(triangles)




    #Write the tecplot header
    line = "VARIABLES=\"X\" \"Y\" \"Z\" \"NodeNum\" \"NodeCode\"" 
    n = 0
    while n < numNeighbours:
        line = line + " \"n%s\"" % (str(n+1))
        n += 1
    
    line = line + "\n"
    outfile.write(line)
    line = "ZONE N=%s E=%s ET=TRIANGLE F=FEPOINT\n" % (str(numNodes),str(numTriangles))
    outfile.write(line)
    
    #Write the points to the tecplot file
    i = 0
    while i < numNodes:
        line = "%s %s %s %s %s" % (str(nodes[i][0]).ljust(15), str(nodes[i][1]).ljust(15), str(nodes[i][2]).ljust(15), str(nodeNumber[i]).ljust(8), str(nodeCode[i]).ljust(3))
        n = 0
        while n < numNeighbours:
            line = line + " %s" % (str(neighbours[i][n]).ljust(8))
            n += 1
            
        line = line + "\n"   
        outfile.write(line)  
        i += 1
        
    #write the triangles
    i = 0
    while i < numTriangles:
        line = "%s %s %s\n" % (str(triangles[i][0]).ljust(8), str(triangles[i][1]).ljust(8), str(triangles[i][2]).ljust(8))
        outfile.write(line)
        i += 1


#def ricom2sts(rcmTimeSeries,sts_in, sts_out, number_of_timesteps, number_of_points, description = ''):

def ricom2sts(sts_out,tsPointsFile,tsFile,tsUFile, epsgOUT = 32755, lat0 = 0.0,long0 = 0.0, starttime=0,endtime=0):

    
    #input: Ricom time series output (ETA, velocity) at a series of points (that define an ANUGA driving boundary)
    #       lat0,long0 - the offset of the grid of the input timeseries - see RCM file

    #output: ANUGA sts driving boundary file for each of the points
    from Scientific.IO.NetCDF import NetCDFFile     #@UnresolvedImport
    from Scientific.IO import NetCDF                #@UnresolvedImport
    import Numeric                                  #@UnresolvedImport
    from osgeo import ogr                           #@UnresolvedImport
    from osgeo import osr                           #@UnresolvedImport
    import csv
    from ricomGeospatialUtilities import ricomLocal2LL
    from ricomGeospatialUtilities import convertXY

    #tsPointsFile FORMAT: X,Y,Z,ElNumRcm,ElNumANUGA
    tsPoints = csv.reader(open(tsPointsFile, 'r'), delimiter=',', quotechar=None)
    tsPoints.next()
        
    points = []
    elevation = []
    permutationANUGA = []
    permutationRICOM = []

    for row in tsPoints:
        points.append([float(row[0]),float(row[1])])
        elevation.append(row[2])
        permutationRICOM.append(row[3])
        permutationANUGA.append(row[4])
        
    
    points = ricomLocal2LL(points,lat0=-46.0,long0=155.5,latoff=0,longoff=0) 
    points = convertXY(points,epsgIN=4326,epsgOUT=epsgOUT)
    
    x = []
    y = []
    
    for row in points:
        x.append(row[0])
        y.append(row[1])
    
    tsdata = open(tsFile, "r").readlines() 
    tsUdata = open(tsUFile, "r").readlines()
    ntsPts = int(tsdata[0].split()[0])
    ntsSteps = int(tsdata[0].split()[1])
    ntsUPts = int(tsUdata[0].split()[0])
    ntsUSteps = int(tsUdata[0].split()[1])
    tsdata.pop(0)
    tsUdata.pop(0)

    eta = []
    xmomentum = []
    ymomentum = []
    times = []
        
    if (ntsUPts == ntsPts and ntsUSteps == ntsSteps):
        timeStep = 0    # timeStep number   
        i = 0           # row in tsfile
        j = 0           # row in tsUfile
        while (timeStep < ntsSteps):
            n = 0           # output point number    
            line = []       
            while n != ntsPts+1:
                row = tsdata[i].split()
                line.extend(row)
                n = n + len(row)
                i += 1    
            
            currentTime =  float(line[0])
            if (currentTime >= starttime and currentTime <= endtime):
                times.append(float(line.pop(0)))
                eta.append(line)
            
                n = 0
                line = []
                while n !=  (ntsUPts*2)+1:
                    row = tsUdata[j].split()
                    line.extend(row)
                    n = n + len(row)        
                    j += 1 
                line.pop(0)
               
                n = 0
                xline = []
                yline = [] 
                while n < ntsUPts*2:
                    xline.append(float(line[n])*float(elevation[n/2]))
                    yline.append(float(line[n+1])*float(elevation[n/2]))
                    n += 2
                xmomentum.append(xline)
                ymomentum.append(yline)
            
            timeStep += 1
    
    time = []
    for row in times:
        time.append(row-starttime)
           
    deltaTime = times[1] - times[0]
    number_of_points = ntsPts
    number_of_timesteps = len(times)
    numbers_in_range = 2

    # NetCDF file definition
    stsout = NetCDFFile(sts_out, 'w')

    #write the header for the sts file
    stsout.institution = 'Australian Tsunami Research Centre'
    stsout.description = 'RiCOM time series to sts driving boundary file (netcdf) for ANUGA'
    stsout.revision_number = 'N/A'
    stsout.starttime = starttime
    #stsout.xllcorner = 426012.084591516
    #stsout.yllcorner = 5157592.30434623
    stsout.xllcorner = 0.0
    stsout.yllcorner = 0.0
    stsout.zone = 55
    stsout.false_easting = 500000
    stsout.false_northing = 10000000
    stsout.datum = 'wgs84'
    stsout.projection = 'UTM'
    stsout.units = 'm'

        
    # Dimension definitions
    stsout.createDimension('number_of_points', number_of_points)
    stsout.createDimension('number_of_timesteps', number_of_timesteps)
    stsout.createDimension('numbers_in_range', numbers_in_range)

    # Variable definitions
    stsout.createVariable('permutation', 'i', ('number_of_points',)).assignValue(permutationANUGA)
    stsout.createVariable('permutationRICOM', 'i', ('number_of_points',)).assignValue(permutationRICOM)
    stsout.createVariable('x', 'd', ('number_of_points',)).assignValue(x)
    stsout.createVariable('y', 'd', ('number_of_points',)).assignValue(y)
    stsout.createVariable('elevation', 'd', ('number_of_points',)).assignValue(elevation)
    stsout.createVariable('time', 'd', ('number_of_timesteps',)).assignValue(time)
    stsout.createVariable('stage', 'd', ('number_of_timesteps','number_of_points')).assignValue(eta)
    stsout.createVariable('stage_range', 'd', ('numbers_in_range',))
    stsout.createVariable('xmomentum', 'd', ('number_of_timesteps','number_of_points')).assignValue(xmomentum)
    stsout.createVariable('xmomentum_range', 'd', ('numbers_in_range',))
    stsout.createVariable('ymomentum', 'd', ('number_of_timesteps','number_of_points')).assignValue(ymomentum)
    stsout.createVariable('ymomentum_range', 'd', ('numbers_in_range',))

    stsout.close()
    
    
#def ngh2pts(ngh_in,pts_out,epsgIN = 32755, epsgOUT = 32755):
def xyz2pts():

    #input: a ngh file from RiCOM
    #output: a pts_file for ANUGA (this is used to interpolate onto the mesh)
    
    from Scientific.IO.NetCDF import NetCDFFile     
    from Scientific.IO import NetCDF
    import Numeric
    from osgeo import ogr
    from osgeo import osr
    import csv
    from ricomGeospatialUtilities import convertXY   
    from ricomGeospatialUtilities import pts_in_polygon   


    sts_out = '../output/anuga/ngh2pts/test.sts'


    # NetCDF file definition
    stsout = NetCDFFile(sts_out, 'w')

    xyz_in = '../input/anuga/ngh2pts/GEBCO_RICOM_in_bounding_polygon.txt'
    xyz_file = csv.reader(open(xyz_in, 'r'), delimiter=',', quotechar=None)
    xyz_file.next()    
    
    points = []
    elevation = []
    for row in xyz_file:
        points.append([row[0],row[1]])
        elevation.append(row[2])
        
    number_of_points = len(elevation)


  

    #write the header for the sts file
    stsout.institution = 'Australian Tsunami Research Centre'
    stsout.description = 'RiCOM time series to sts driving boundary file (netcdf) for ANUGA'
    stsout.xllcorner = 0.0
    stsout.yllcorner = 0.0
    stsout.zone = -1
    stsout.false_easting = 500000
    stsout.false_northing = 10000000
    stsout.datum = 'wgs84'
    stsout.projection = 'UTM'
    stsout.units = 'm'
        
    # Dimension definitions
    stsout.createDimension('number_of_points', number_of_points)
    stsout.createDimension('number_of_dimensions', 2)

    stsout.createVariable('points', 'd', ('number_of_points','number_of_dimensions')).assignValue(points)
    stsout.createVariable('elevation', 'd', ('number_of_points',)).assignValue(elevation)

    stsout.close()

#    nghfile = open('../input/anuga/ngh2pts/NZ_S1S2S3_RD_LidarTEST.xyz', "r").readlines() 
#
#    nghPoints = []
#    for row in nghfile:
#        point = row.split()
#        nghPoints.append([float(point[1]),float(point[2]),float(point[4])])
#    
#    bounds_file = '../input/anuga/ngh2pts/bounding_polygon.csv'
#    polyfile = csv.reader(open(bounds_file, 'r'), delimiter=',', quotechar=None)
#    outfile = open('../input/anuga/ngh2pts/GEBCO_RICOM_in_bounding_polygon.csv', "w")    
#
#    boundingPoly = []
#    for row in polyfile:
#        boundingPoly.append([float(row[0]),float(row[1])])
#
#    
#    pointsInside = []
#    pts_in_polygon(boundingPoly, pointsAll=nghPoints, pointsInside=pointsInside)
#
#      
#    for row in pointsInside:
#        line = '%s,%s,%s\n' % (str(row[0]),str(row[1]),str(row[2]))
#        outfile.write(line)
#        
#        
#        
#    gebco_file = '/Users/tbone/Documents/PhD/RiCOM/Grids/Data/GEBCO/AUSTAS30sec_UTM55.asc'
#    ptfile = csv.reader(open(gebco_file, 'r'), delimiter=' ', quotechar=None)
#    ptfile.next()    
#
#    points = []
#    for row in ptfile:
#        points.append([float(row[0]),float(row[1]),float(row[2])])
#    
#    pointsInside = []
#    pts_in_polygon(boundingPoly, pointsAll=points, pointsInside=pointsInside)        
#        
#    for row in pointsInside:
#        line = '%s,%s,%s\n' % (str(row[0]),str(row[1]),str(row[2]))
#        outfile.write(line)
#        
#    outfile.close()
#    

#input grid is in ESRI ASCII grid format
def GEBCOtoTEC(gebco_in,tec_out):
    
    
    outfile = open(tec_out, "w")    

    ascfile = open(gebco_in, "r").readlines() 
    nx = int(ascfile[0].split()[1])
    ny = int(ascfile[1].split()[1])
    xllcenter = float(ascfile[2].split()[1])
    yllcenter = float(ascfile[3].split()[1])
    cell_size = float(ascfile[4].split()[1])
    
    yulcenter = yllcenter + (ny-1)*cell_size
    i = 0
    j = 0
    
    #remove the header
    ascfile.pop(0)
    ascfile.pop(0)
    ascfile.pop(0)
    ascfile.pop(0)
    ascfile.pop(0)
  
    numVolumes = (nx-1) * (ny-1) * 2
    yulcenter = yllcenter + (ny-1)*cell_size

    #Write the tecplot header
    line = "VARIABLES=\"X\" \"Y\" \"eta\"\n" 
    outfile.write(line)
    line = "ZONE N=%s E=%s ET=TRIANGLE F=FEPOINT\n" % (str(nx*ny),str(numVolumes))
    outfile.write(line)
    
    #Write the points to the tecplot file
    i = 0
    j = 0
    
    
    
    while j < ny:
        i = 0
        y = yulcenter - j*cell_size
        depths = ascfile[j].split()
        while i < nx:
            z = -float(depths[i])
            x = xllcenter + i*cell_size
            outfile.write("%s %s %s\n" % (str(x).ljust(11),str(y).ljust(11),str(z)))
            i += 1 
        j += 1
    
    #write the volumes to the tecplot file
    j = 0
    while j < (ny - 1):         #iterate through the rows of the grid
        i = 0
        while (i < nx - 1):         #iterate through the columns of the grid
            # points of the square element relating to the current grid point
            #  -   -   -   -   -   -   -
            #  -   -   -   -   p1  p4  -
            #  -   -   -   -   p2  p3  -
            #  -   -   -   -   -   -   -
            p1 = j*nx + i
            p2 = (j + 1)*nx + i 
            p3 = (j + 1)*nx + i + 1
            p4 = j*nx + i + 1
            #Each point has two elements (volumes) associated with it
            #add one to each point index because ANUGA appears to base 0.  Tecplot doesn't accepts 0 indices
            v1 = [p1+1, p2+1, p3+1]
            v2 = [p1+1, p3+1, p4+1]
            line1 = "%s %s %s\n" % (str(v1[0]).ljust(8), str(v1[1]).ljust(8), str(v1[2]))
            line2 = "%s %s %s\n" % (str(v2[0]).ljust(8), str(v2[1]).ljust(8), str(v2[2]))
            outfile.write(line1)
            outfile.write(line2)
            i += 1
        j += 1
    
    
    outfile.close()

#input grid is in ESRI ASCII grid format
def GEBCOtoMOST(gebco_in,most_out):
    
    
    outfile = open(most_out, "w")    

    ascfile = open(gebco_in, "r").readlines() 
    ncols = int(ascfile[0].split()[1])
    nrows = int(ascfile[1].split()[1])
    xllcenter = float(ascfile[2].split()[1])
    yllcenter = float(ascfile[3].split()[1])
    cell_size = float(ascfile[4].split()[1])
    
    yulcenter = yllcenter + (nrows-1)*cell_size
    i = 0
    j = 0
    
    #remove the header
    ascfile.pop(0)
    ascfile.pop(0)
    ascfile.pop(0)
    ascfile.pop(0)
    ascfile.pop(0)
    
#    ascfile.reverse()

    cam = 1
    
    outfile.write('%s\n' % str(ncols))
    outfile.write('%s\n' % str(nrows))
    
    while j < ncols:
        outfile.write('%s\n' % str(xllcenter+j*cell_size))   
        j += 1
    
    while i < nrows:
        outfile.write('%s\n' % str(yulcenter-i*cell_size))
        i += 1
        
    for row in ascfile:
        row = row.split()
        for pt in row:
            outfile.write('%s   ' % str(-float(pt))) 
            #outfile.write('%s   ' % pt)           
          
        outfile.write('\n')
    
    outfile.close()
        
    
#    while n < nrows:
#        i = 0
#        row = ascfile[n+6].split()
#        while i < ncols:
#            #start in the top left corner (i = 0, j = nrows)
#            x = float(xllcenter + i*cell_size)
#            y = float(yllcenter + (j-1)*cell_size)
#            z = round(float(row[i]),4)
#
#            points.append([x,y,z])
#            i = i + 1
#        j = j - 1
#        n += 1    
#    
#    
#
#    ptfile = csv.reader(open(gebco_in, 'r'), delimiter=' ', quotechar=None)
#    outfile = open(most_out, "w")    
#
#    i = 0
#    while i<12:         #expecting a 12 line header for the GEBCO ascii file
#        ptfile.next()
#        i+=1
#    
#    
#    for pt in ptfile:
#        x = float(pt[0])
#        y = float(pt[1])
#        z = float(pt[2])
#
#
#


def tec_convertCoords(tec_in,epsgOUT = 4326, lat0=0,long0=0,latoff=0,longoff=0):
	'''
	Convert the coordinates of a Fulldomain output file from Ricom
	
	Ricom local Coords -> epsgOUT
	
	'''

	from ricomGeospatialUtilities import ricomLocal2LL
	from ricomGeospatialUtilities import convertXY
	
	# VARIABLES="X" "Y" "Z" "ETA" "U" "V" 
	#ZONE N= 395116 E= 783117 ET=TRIANGLE F=FEPOINT
	#  SOLUTIONTIME=  0.000000000000000E+000
	tec_file = open(tec_in, "r").readlines()
	outfile = open("Fulldomain2452.dat", "w")
	i = 0
	j = 0
	k = 0
	num_nodesNGH = 0
	num_nodesTEC = 0
	
	for line in tec_file:

		i = line.find("N=")
		j = line.find("E=")
		if i > 0:
			num_nodesTEC = int(line[i+3:j])
			k +=2 
			break
		k+=1
	
	i = 0
	
	while i < num_nodesTEC:
		
		line = tec_file[k].split()
		
		xy = []
		xy.append([float(line[0]),float(line[1])])
		xyLL = ricomLocal2LL(xy,37.838428, 143.22166)
		xy2452 = convertXY(xyLL,4326,2452)
		line[0] = xy2452[0][0]
		line[1] = xy2452[0][1]
		line2 = "%s    %s    %s    %s    %s    %s\n" % (line[0],line[1],line[2],line[3],line[4],line[5])
		
		tec_file[k] = line2	
		k+=1
		i+=1
		
	for line in tec_file:
		outfile.write(line) 
		
def tec2gdal(tec_in, asc_out, xmin, xmax, ymin, ymax, step):
    
	'''
	tec_in - the z values of a tecplot rectangular grid 
	Write Data File -> Field Data, Pt, ASCII, Rectangular Zone, Z
	
	Tecplot Header:
		TITLE     = ""
		VARIABLES = "Z"
		ZONE T="Rectangular zone"
		STRANDID=0, SOLUTIONTIME=0
		I=3000, J=3950, K=1, ZONETYPE=Ordered
		DATAPACKING=POINT
		DT=(SINGLE )
		
	ASCII Grid Header:
	
		ncols 		1989
		nrows        926
		xllcorner    140.889965281893
		yllcorner    38.150035995522
		cellsize     0.000040063538
		NODATA_value -3.4028234663852885981e+38
	'''
	tecIn = open(tec_in, "r").readlines() 
	outfile = open(asc_out, "w")    
	
	
	line = tecIn[4]
	i = line.find('I=')
	j = line.find(', J=')
	k = line.find('J=')
	l = line.find(', K=')
	
	numX = int(line[i+2:j])
	numY = int(line[k+2:l])
	
	
	depths = []
	x = []
	y = []
	
	points = []
	numPts = numX*numY
	print numX
	print numY
	print numPts
	
	outfile.write("ncols         %s\n" % numX)
	outfile.write("nrows         %s\n" % numY)
	outfile.write("xllcorner     %s\n" % xmin)
	outfile.write("yllcorner     %s\n" % ymin)
	outfile.write("cellsize      %s\n" % step)
	outfile.write("NODATA_value  %s\n" % -9999) 
	i=0
	j=0
	print "hello!!! =  %s" % tecIn[7+i+(j*i)]
	
	
	asc_lines = []
	while j < numY:
		line = ""
		i=0
		while i < numX:
			z = float(tecIn[7+i+(j*numX)].split()[0])
			line = line + " %s" % z
			i += 1
	
		asc_lines.append(line)
		j+=1
		
	asc_lines.reverse()
	for line in asc_lines:
		outfile.write(line + "\n")

	outfile.close()



def tec2mostGRD(tec_in, most_out):
    
    from ricomGeospatialUtilities import convertXY

    
    tecIn = open(tec_in, "r").readlines() 
    outfile = open(most_out, "w")    
    
    
    line = tecIn[6]
    i = line.find('I=')
    j = line.find(', J=')
    k = line.find('J=')
    l = line.find(', K=')
   
    numX = int(line[i+2:j])
    numY = int(line[k+2:l])

  
    depths = []
    x = []
    y = []
    
    points = []
    numPts = numX*numY
    n=0
    while n < numPts:
        line = tecIn[9+n].split()
        points.append([float(line[0]),float(line[1]),float(line[2])])
        n+=1
    #points = convertXY(points,32755,4326)
    
    j = 0       #data starts on line # 9
    i = 0       #data ordered in lines along x

    while j < numY:
        i = 0
        line = []
        while i < numX:
            pt = points[(j*(numX))+i]
            if j == 0:
                x.append(pt[0])
            
            line.append(pt[2])
            i+=1

        y.append(pt[1])
        depths.append(line)
        j+=1
        
    
    
    y.reverse()
    depths.reverse()    
    
    line = "%i %i\n" % (numX,numY)
    outfile.write(line)
    for pt in x:
        outfile.write("%f\n" % (pt))

    for pt in y:
        outfile.write("%f\n" % (pt))
        
    for line in depths:
        for pt in line:
            outfile.write("%f   " % (pt))
        outfile.write("\n")

    outfile.close()


def postGISTest(database='kingstonbeach'):
    
    import psycopg2
    import sys
    from ricomGeospatialUtilities import convertXY

    
    try:
        conn = psycopg2.connect("dbname='kingstonbeach' user='tbone'");
    except:
        print "Unable to connect to the database" 
        sys.exit()   

    cur = conn.cursor()

   #cur.execute("SELECT ST_AsText(points) AS geom FROM model_pts WHERE ST_Distance(points, ST_GeomFromText('POINT(526696.366853242 5241526.92119538)', 32755)) < 200;")
    
    
    
#    cur.execute("SELECT m.ID, f.gid, ST_AsText(m.points) FROM footprints AS f, model_pts AS m \
#                WHERE ST_Distance(m.points,f.the_geom) < 100 ORDER BY f.gid;")

#
#    cur.execute("SELECT ID, etamax, z FROM triangles WHERE z > 0 AND etamax > 0.5 ORDER BY z;")
#    footprints = cur.fetchall()
#    for row in footprints:
#        print "%s, %s, %s" % (str(row[0]),str(row[1]),str(row[2]))

#    cur.execute("SELECT t.ID, ST_AsText(t.element), ST_Distance(t.element,ST_GeomFromText(%s,32755)) as distance FROM triangles AS t \
#                WHERE ST_Distance(t.element,ST_GeomFromText(%s,32755)) < 40 ORDER BY distance;", (footprint[0],footprint[0]))
 
       
        
#    cur.execute("SELECT ele.id, ST_AsText(ele.element), f.gid as order_id, ele.v1, ele.v2, ele.v3 as dist FROM buildingele AS ele, footprints AS f \
#                WHERE ST_Intersects(ele.element,f.the_geom) ORDER BY order_id;")

    cur.execute("SELECT ST_AsText(f.the_geom), ST_AsText(ST_PointOnSurface(f.the_geom)) FROM footprints AS f, boundingpoly AS p WHERE ST_Contains(p.polygon,f.the_geom);")
    
    footprints = cur.fetchall()
        
    #cur.execute("SELECT ST_AsText(ST_ClosestPoint('POINT(526696.366853242 5241526.92119538)',points) FROM model_pts;")

    #cur.execute("SELECT ST_AsText(points) AS geom FROM model_pts WHERE ST_ClosestPoint(points, ST_GeomFromText('POINT(526696.366853242 5241526.92119538)', 32755));")


    rows = cur.fetchall()
    n = 0
    print "\nShow me the result:\n"
    for row in rows:
        n+=1
        print "%s, %s, %s, %s, %s" % (str(row[2]), str(row[0]), str(row[3]), str(row[4]), str(row[5]))


    print "\n%s" % (str(n))


#    print "Points index...\n"
#
#    cur.execute("CREATE INDEX pts_idx ON model_pts USING GIST (points);")
#
#    print "Element index...\n"
#
#
#    cur.execute("CREATE INDEX tri_idx ON triangles USING GIST (element);")
#
#    print "Footprints index...\n"
#
#    
#    cur.execute("CREATE INDEX footprint_idx ON footprints USING GIST (the_geom);")


#    cur.execute("CREATE TABLE model_pts (ID int4, z float);")
#    cur.execute("SELECT AddGeometryColumn( 'model_pts', 'points', 32755, 'POINT', 2);")



    # cur.execute("""INSERT INTO some_table (an_int, a_date, a_string) VALUES (%s, %s, %s);""", (10, datetime.date(2005, 11, 18), "O'Reilly"))



    #cur.execute("INSERT INTO roads_geom (ID, NAME, geom) VALUES (1,'Test Point',ST_GeomFromText('POINT(191232 243118)',32755));")
    #cur.execute("DELETE FROM model_pts WHERE ID=2;")
    #cur.execute("CREATE TABLE roads_geom ( ID int4, NAME varchar(25) )")
    #cur.execute("SELECT AddGeometryColumn( 'roads_geom', 'geom', 32755, 'POINT', 2)")
    
    
    conn.commit()
    cur.close()
    conn.close()
    
    
def pg_addRicomGrid(ngh_in,tri_in,tec_in, database='camtest'):
    
    
    
    import psycopg2
    import sys
    from ricomGeospatialUtilities import convertXY




     
    ngh_file = open(ngh_in, "r").readlines() 
    tri_file = open(tri_in, "r").readlines() 
    
    triangles = []                  # List of triangles [point1 point2 point3]  - point indexes correspond to there position in points[]    - YES
    numNodes = int(ngh_file[0].split()[0])

    nodes = []                
    n = 0
    while n < numNodes:
        vertLine = ngh_file[n+3].split()
        node = [float(vertLine[1]), float(vertLine[2]), float(vertLine[4])]
        nodes.append(node)
        n += 1     

    nodesUTM = convertXY(nodes,4326,32755)



    # Triangles
    #--------------------------------------------------------------------------------
    for line in tri_file:
        line_list = line.split() 
        triangles.append([int(line_list[0]), int(line_list[1]), int(line_list[2])]) 
    
    
    print "Reading Tecplot file (ETAMAX)...\n" 


    # Assuming the following tecplot header:
    # ---------------------------------------
    # VARIABLES="X" "Y" "Z" "ETAMAX" 
    # ZONE N=      509544  E=      248605
    # ZONETYPE=FETRIANGLE DATAPACKING=BLOCK
    # VARLOCATION=([4]=CELLCENTERED) 

    tec_file = open(tec_in, "r").readlines() 

    ETAMAX = [] 
    line = tec_file[1]
    i = line.find('N=')
    j = line.find('E=')
   
    numNodes = int(line[i+2:j])
    numElements = int(line[j+2:])
    
    #X, Y, Z are node values, ETAMAX is an element value
    
    linesNodes  = int(numNodes/6)
    if (numNodes%6 > 0):
        linesNodes += 1
    linesElements = int(numElements/6) 
    if (linesElements%6 > 0):
        linesElements += 1
    
    n = 0    
    while n < linesElements:
        lineEta = tec_file[4+n+linesNodes*3].split()
        
        for pt in lineEta:
            ETAMAX.append(float(pt))
            
        n+=1
        
    try:
        conn = psycopg2.connect("dbname='camtest' user='tbone'");
    except:
        print "Unable to connect to the database" 
        sys.exit()   

    cur = conn.cursor()


    print "Creating new tables..." 


    cur.execute("DROP TABLE triangles;")
    cur.execute("DROP TABLE model_pts;")

    cur.execute("CREATE TABLE model_pts (ID int4, z float);")
    cur.execute("SELECT AddGeometryColumn( 'model_pts', 'points', 32755, 'POINT', 2);")

    cur.execute("CREATE TABLE triangles (ID int4, v1 int4, v2 int4, v3 int4, etamax float, z float);")
    cur.execute("SELECT AddGeometryColumn( 'triangles', 'element', 32755, 'POLYGON', 2);")


    print "Inserting model points..." 



    n = 0
    for row in nodesUTM:
        cur.execute("INSERT INTO model_pts (ID, points, z) VALUES (%s,ST_GeomFromText('POINT(%s %s)',32755),%s);", (n+1, row[0], row[1], row[2]))
        n+=1
    
    
    print "Inserting triangles..." 
    
    
    n = 0
    for row in triangles:
        zAvg = (nodesUTM[row[0]-1][2] + nodesUTM[row[1]-1][2] + nodesUTM[row[2]-1][2]) / 3
        cur.execute("INSERT INTO triangles (ID, v1, v2, v3, etamax, z, element) VALUES (%s,%s,%s,%s,%s,%s,ST_GeomFromText('POLYGON((%s %s, %s %s, %s %s, %s %s))',32755));", \
                            (n+1, row[0], row[1], row[2],ETAMAX[n],zAvg, nodesUTM[row[0]-1][0],nodesUTM[row[0]-1][1],nodesUTM[row[1]-1][0],nodesUTM[row[1]-1][1], \
                             nodesUTM[row[2]-1][0],nodesUTM[row[2]-1][1],nodesUTM[row[0]-1][0],nodesUTM[row[0]-1][1]))   
        n += 1

    
    
    
   
    conn.commit()
    cur.close()
    conn.close()    
    
    
    
def tecETAMAX_to_postgis(tec_in, database = 'camtest'):

    import psycopg2
    import sys


    print "Reading Tecplot file (ETAMAX)..." 


    # Assuming the following tecplot header:
    # ---------------------------------------
    # VARIABLES="X" "Y" "Z" "ETAMAX" 
    # ZONE N=      509544  E=      248605
    # ZONETYPE=FETRIANGLE DATAPACKING=BLOCK
    # VARLOCATION=([4]=CELLCENTERED) 

    tec_file = open(tec_in, "r").readlines() 

    X = []
    Y = []
    Z = []
    ETAMAX = []
    elements = []
 
    line = tec_file[1]
    i = line.find('N=')
    j = line.find('E=')
   
    numNodes = int(line[i+2:j])
    numElements = int(line[j+2:])
    
    #X, Y, Z are node values, ETAMAX is an element value
    
    linesNodes  = int(numNodes/6)
    if (numNodes%6 > 0):
        linesNodes += 1
    linesElements = int(numElements/6) 
    if (linesElements%6 > 0):
        linesElements += 1
        
    n = 0
    while n < linesNodes:
        lineX = tec_file[4+n].split()
        lineY = tec_file[4+n+linesNodes].split()
        lineZ = tec_file[4+n+linesNodes*2].split()
        
        for pt in lineX:
            X.append(float(pt))
        for pt in lineY:
            Y.append(float(pt))
        for pt in lineZ:
            Z.append(float(pt))
            
        n+=1
    
    n = 0    
    while n < linesElements:
        lineEta = tec_file[4+n+linesNodes*3].split()
        
        for pt in lineEta:
            ETAMAX.append(float(pt))
            
        n+=1
        
    n = 0
    while n < numElements:
        lineEle = tec_file[4+n+linesNodes*3+linesElements].split()
        elements.append([int(lineEle[0]),int(lineEle[1]),int(lineEle[2])])
        n+=1

    
    try:
        conn = psycopg2.connect("dbname='camtest' user='tbone'");
    except:
        print "Unable to connect to the database" 
        sys.exit()   

    cur = conn.cursor()


    cur.execute("ALTER TABLE triangles ADD COLUMN etamax float;")


    print "Inserting triangles Data..." 



    n = 1
    for row in ETAMAX:
        sql =  "INSERT INTO triangles (etamax) VALUES (%s);" % (str(row))
        cur.execute(sql)
        
        if n%1000 == 0:
            print "%s  of  %s,   row = %s" % (str(n), str(numElements), str(row))

        n+=1
      
      
    print "Number of ETA points = %s" % (str(n))
    
    conn.commit()
    cur.close()
    conn.close()  


def pg_index(database = 'camtest'):

    import psycopg2
    import sys
 
    try:
        conn = psycopg2.connect("dbname='camtest' user='tbone'");
    except:
        print "Unable to connect to the database" 
        sys.exit()   

    cur = conn.cursor()
      
    print "Building points index..."
    cur.execute("CREATE INDEX pts_idx ON model_pts USING GIST (points);")
    print "Building triangles index...\n"
    cur.execute("CREATE INDEX tri_idx ON triangles USING GIST (element);")
      
    conn.commit()
    cur.close()
    conn.close()  
    
def netcdfTEST():
    
    #output: ANUGA sts driving boundary file for each of the points
    from Scientific.IO.NetCDF import NetCDFFile     
    from Scientific.IO import NetCDF
    
    netcdf_in = NetCDFFile("../input/netcdf/test.nc", 'r')
    
    
    x = netcdf_in.variables['x'][:]
    points = netcdf_in.variables['points'][:]
    elements = netcdf_in.variables['elements'][:]
    
    line = [points[0][0],points[0][1],points[0][2]]
    
    number_of_points = netcdf_in.dimensions['number_of_points']
    number_of_vertices = netcdf_in.dimensions['number_of_verticies']
    number_of_elements = netcdf_in.dimensions['number_of_elements']
    
    cam = 0
    
    
    
#create .POLY file from postGIS footprint table
#to be used with the Triangle utility
def createPOLY(holes=1):
    
    import psycopg2
    import sys
        
    try:
        conn = psycopg2.connect("dbname='kingstonbeach' user='tbone'");
    except:
        print "Unable to connect to the database" 
        sys.exit()   

    cur = conn.cursor()
    
    cur.execute("SELECT ST_AsText(f.the_geom), ST_AsText(ST_PointOnSurface(f.the_geom)) FROM footprints AS f, boundingpoly AS p WHERE ST_Contains(p.polygon,f.the_geom);")
    
    footprints = cur.fetchall()

    polygonList = []
    pointInsideList = []
    #extract the polygon verices for each footprint    
    for poly in footprints:
        i = poly[1].find('POINT(')
        j = poly[1].find(')')
        pointInside = poly[1][i+6:j]
        pointInside = pointInside.replace(',',' ').split()
        pointInsideList.append(pointInside)
        i = poly[0].find('POLYGON((')
        j = poly[0].find('))')
        points = poly[0][i+9:j]
        points = points.replace(',',' ').split()
        points.pop()            #remove repeated last item in polygon 
        points.pop()
        i = 0
        polygon = []
        while i < len(points):
            polygon.append([points[i], points[i+1]])
            i += 2
        polygonList.append(polygon)
        

    ply_file = open("../output/SQL/KBBuildingBoundingPoly2.ply", "r").readlines() 
    boundingPoly = []
    for row in ply_file:
        row =  row.split()
        boundingPoly.append(row)
    

    #remove repeated last item in polygon
    boundingPoly.pop()               
    
    outfile = open("../output/SQL/kingstonbeach5.poly", "w")    
    
    
    numVerticies = len(boundingPoly)
    
    for poly in polygonList:
        numVerticies = numVerticies + len(poly)
    
    numHoles = len(pointInsideList)


    outfile.write("%s 2 0 1\n" % (str(numVerticies)))
    
    
    segments = []
    #write the vertices in the footprint and bounding poly lists
    i = 1
    for pt in boundingPoly:
        ptYshift = float(pt[1]) - 5000000
        outfile.write("%s   %s    %s   1\n" % (str(i), pt[0], str(ptYshift)))
        if i != len(boundingPoly):
            segments.append([i,i+1,1])
            
        else:
            index = i - len(boundingPoly) + 1
            segments.append([i,index, 1])
        
        i += 1

    polyIndex = 2
    for poly in polygonList:
        lengthPoly = len(poly)
        j = 1
        for pt in poly:
            ptYshift = float(pt[1]) - 5000000
            outfile.write("%s   %s    %s    %s\n" % (str(i), pt[0], str(ptYshift), str(polyIndex)))
           
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
    if(holes == 1):
        outfile.write("%s\n" % (str(len(pointInsideList))))
        i = 1
        for pt in pointInsideList:
            ptYshift = float(pt[1]) - 5000000
            outfile.write("%s   %s   %s\n" % (str(i), str(pt[0]), str(ptYshift)))
            i += 1
    else:
        outfile.write("0\n")


    
    outfile.close()
    
    lengthInside = len(pointInsideList)
    lengthPoly = len(polygonList)

def Triangle2Tecplot():
    
    node_file = open("../output/SQL/KingstonBeach5/kingstonbeach5.1.node", "r").readlines() 
    ele_file = open("../output/SQL/KingstonBeach5/kingstonbeach5.1.ele", "r").readlines() 
    outfile = open("../output/SQL/KingstonBeach5/kingstonbeach5.1.dat", "w")    

    nodeHeader = node_file[0].split()
    node_file.pop(0)
    node_file.pop()
    
    eleHeader = ele_file[0].split()
    ele_file.pop(0)
    ele_file.pop()
     
    #Write the tecplot header
    line = "VARIABLES=\"X\" \"Y\" \"Z\"\n" 
    outfile.write(line)
    line = "ZONE N=%s E=%s ET=TRIANGLE F=FEPOINT\n" % (nodeHeader[0],eleHeader[0])
    outfile.write(line)
    
    #Write the points to the tecplot file
    for pt in node_file:
        pt = pt.split()
        y = float(pt[2])
        outfile.write("%s     %s     0.0\n" % (pt[1],str(y)))
    

    #Write the elements to the tecplot file
    for pt in ele_file:
        pt = pt.split()
        outfile.write("%s     %s     %s\n" % (pt[1],pt[2],pt[3]))    



def TriangleExtractInnerBoundaries(boundary_poly,filename, epsgBP, triangle_poly_file, triangle_node_file,epsgTRI, holes = 1):
	'''
	
	'''
	
	
	node_file = open(triangle_node_file, "r").readlines() 
	poly_file = open(triangle_poly_file, "r").readlines() 
	bp_file = open(boundary_poly, "r").readlines() 
	
	
	node_outfile = open(filename + ".node", "w")    
	poly_outfile = open(filename + ".poly", "w")    
	
	boundingPoly = []
	for row in bp_file:
		row =  row.split()
		boundingPoly.append(row)  
	
	boundingPoly = convert_points(boundingPoly, epsgBP,epsgTRI)
	boundingPoly.pop()		#remove repeated last item in the bounding polygon
	
	number_of_bounding_poly_pts = len(boundingPoly)

	'''
	1) Read the triangle .node file and extract all internal boundary nodes (i.e. code > 2)
	2) Create an array that maps the existing node code the new node code (i.e. with the new
		outer boundary and no internal nodes)
	'''
	nodeHeader = node_file[0].split()
	numNodes = int(nodeHeader[0])
	n = 0			#index for the line number of the node file
	nodesInside = []
	indexMap = []	#maps the old node indices to the new indices
	i = 0			#index for the internal nodes

	while n < numNodes:
		node = node_file[n+1].split()
		x = float(node[1])
		y = float(node[2])
		code = int(node[3])
		index =  int(node[0])
		
		if code > 1:	#(i.e. an internal boundary node)
			newIndex = i+1+number_of_bounding_poly_pts
			nodesInside.append([newIndex,x,y,code])
			indexMap.append(newIndex)
			i+=1
		else:
			indexMap.append(0)
			
		n+=1
		
		
	number_of_internal_boundary_nodes = i
	number_of_nodes = number_of_internal_boundary_nodes + number_of_bounding_poly_pts
	
	
	segmentsInside = []
	polyHeader = poly_file[0].split()
	i_line = int(polyHeader[0]) + 1
	numPolySegments = int(poly_file[i_line].split()[0])
	i_line += 1
	i = 0
	while i < numPolySegments:
		seg = poly_file[i_line].split()
		n1 = int(seg[1])
		n2 = int(seg[2])
		code = int(seg[3])
		if code > 1:
			n1_new = indexMap[n1-1]
			n2_new = indexMap[n2-1]
			segmentsInside.append([n1_new,n2_new,code])
		i_line += 1
		i+=1

	'''
	write the Bounding Polygon vertices to the output .node file and the segments 
	to the output .poly file
	'''
	#write the .node file header
	node_outfile.write("%s 2 0 1\n" % (str(number_of_nodes)))
	i = 1
	segmentsBoundingPoly = []
	for pt in boundingPoly:
		node_outfile.write("%s   %s    %s   1\n" % (str(i), pt[0], pt[1]))
		if i != len(boundingPoly):
			segmentsBoundingPoly.append([i,i+1,1])
			
		else:
			index = i - len(boundingPoly) + 1
			segmentsBoundingPoly.append([i,index, 1])
		
		i += 1
		
	'''
	Write the internal nodes to the .node file
	'''
	for pt in nodesInside:
		node_outfile.write("%s   %s    %s   %s\n" % (pt[0], pt[1], pt[2], pt[3]))
	
	'''
	Write the .poly file header - nodes defined in the .node file
	'''
	number_of_segments = len(segmentsInside) + len(segmentsBoundingPoly)
	poly_outfile.write("0 2 0 1\n")
	poly_outfile.write("%s 1\n" % number_of_segments)
	
	'''
	Write the segments to the .poly file
	'''
	i = 1
	#write the boundary segments
	for seg in segmentsBoundingPoly:
		poly_outfile.write("%s   %s   %s   %s\n" % (str(i), str(seg[0]), str(seg[1]), str(seg[2])))
		i += 1
		
	#write the internal (buildings) segments
	
	for seg in segmentsInside:
		poly_outfile.write("%s   %s   %s   %s\n" % (str(i), str(seg[0]), str(seg[1]), str(seg[2])))
		i += 1
	
	number_of_holes = int(poly_file[i_line].split()[0])
	#write number of holes to .poly file
	poly_outfile.write("%s\n" % number_of_holes)
 	i_line += 1
 	i = 0
 	while i <  number_of_holes:
 		poly_outfile.write(poly_file[i_line])
 		i_line += 1
 		i += 1
 	
 	
	node_outfile.close()
	poly_outfile.close()
	
		
	
#convert the output of the Triangle utility to Ricom ngh/el files

def Triangle2RiCOM(node_in, ele_in, poly_in, edge_in, ngh_out, el_out,nod_out, holes = 1):


    node_file = open(node_in, "r").readlines() 
    ele_file = open(ele_in, "r").readlines() 
    poly_file = open(poly_in, "r").readlines() 
    edge_file = open(edge_in, "r").readlines() 
    
    
    ngh_outfile = open(ngh_out, "w")    
    el_outfile = open(el_out, "w")    

    nod_outfile = open(nod_out, "w")    


    nodeHeader = node_file[0].split()
    node_file.pop(0)
    node_file.pop()
    
    eleHeader = ele_file[0].split()
    ele_file.pop(0)
    ele_file.pop()

    polyHeader = poly_file[1].split()
    poly_file.pop(0)
    poly_file.pop(0)
    poly_file.pop()
    
    edgeHeader = edge_file[0].split()
    edge_file.pop(0)
    edge_file.pop()
    
    
    vertices = []
    elements = []
    outerSegments = []
    outerPoly = []
    #get the vertices
    #each line of vertices array - x,y,NodeCode,z
    for vert in node_file:
        vert = vert.split()
        vertices.append([float(vert[1]),float(vert[2]), int(0), float(0.0)])
    
    #get the elements  
    for ele in ele_file:
        ele = ele.split()
        elements.append([ele[1],ele[2],ele[3]])
    
    #get the neighbours
    for edge in edge_file:
        edge = edge.split()
        vertices[int(edge[1])-1].append(int(edge[2]))
        vertices[int(edge[2])-1].append(int(edge[1]))
        
    
    number_of_neighbours = 0
    cam = 0
    for vert in vertices:
        neighs = len(vert) - 4
        if neighs > number_of_neighbours:
            number_of_neighbours = neighs
    
    
    #write the node codes to the vertices list
    i = 0
    numSegments = int(polyHeader[0])
    while i < numSegments:
        seg = poly_file[i].split()
        code = int(seg[3])
        if code == 1:     #outer boundary
            outerSegments.append([int(seg[1]),int(seg[2])])
            nodeCode = 1
        elif code > 1 and holes == 1:              #inner boundary
            nodeCode = 2
        else:
            nodeCode = 0
            
        vertices[int(seg[1])-1][2] = nodeCode
        vertices[int(seg[2])-1][2] = nodeCode
        i += 1


    #Construct the outer polygon (for use with griding in Ricom)
    #Assuming that all the segments in outerSegments are out of order - need to create an ordered polygon        

    numBoundaryPts = len(outerSegments)
    n1 = outerSegments[0][0]
    n2 = outerSegments[0][1]
    v1 = vertices[n1-1]             #Starting segment
    v2 = vertices[n2-1]
    outerPoly.append([v1[0],v1[1]])
    
    i = 2
    while i < numBoundaryPts:
        outerPoly.append([v2[0],v2[1]])
        neighs = len(v2) - 4
        n = 0
        while n < neighs:
            n3 = v2[4+n]
            v3 = vertices[n3-1]
            if (v3[2] == 1) and (n3 != n1):
                n1 = n2
                v2 = v3
                n2 = n3
                break
            else:
                n += 1
        i += 1
    
     
    #write the outer boundary polygon as a nodefile 
    nod_outfile.write("%s\n1     0\n%s\n" % (str(len(outerPoly)+1),str(len(outerPoly))))        
    for vert in outerPoly:
        nod_outfile.write("%s       %s       0.0\n" % (str(vert[0]),str(vert[1])))
        
    nod_outfile.write("0")
    
    nod_outfile.close()
    
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
    
    
def convertNGHcoords(ngh_in, ngh_out, epsgIN = 28355, epsgOUT = 4326, xShift=0,yShift=0):


    from ricomGeospatialUtilities import convertXY
    
    ngh_fileIn = open(ngh_in, "r").readlines() 
    ngh_fileOut = open(ngh_out, "w")


    numNodes = int(ngh_fileIn[2].split()[0])
    numNeighs = int(ngh_fileIn[3].split()[0])
    
    ngh_fileOut.write(ngh_fileIn[0])
    ngh_fileOut.write(ngh_fileIn[1])
    ngh_fileOut.write(ngh_fileIn[2])
    ngh_fileOut.write(ngh_fileIn[3])

    ngh_fileIn.pop(0)
    ngh_fileIn.pop(0)
    ngh_fileIn.pop(0)
    ngh_fileIn.pop(0)

    points = []
    for vert in ngh_fileIn:
        vert = vert.split()
        points.append([float(vert[1]),float(vert[2]),float(vert[4])])
    
    
    cam = len(ngh_fileIn)
   
    if(epsgIN != epsgOUT):    
        points = convertXY(points,epsgIN,epsgOUT)
    
    
        #write the ngh file header
#    ngh_fileOut.write("%s\n%s\n" % (str(numNodes),str(numNeighs)))
 #   ngh_fileOut.write("0.0     0.0     0.0     0.0\n")
    i = 0
    #write vertices
    for pt in ngh_fileIn:
        pt = pt.split()
        code = int(pt[3])
        if (code != 0) and (code != 1) and (code != 2) and (code != 7):            #remove bad codes (Trigrid seems to add weird numbers)
            pt[3] = '0' 
        
        line = "%s %s %s %s %s   " % (str(pt[0]).ljust(8), str(points[i][0]+xShift).ljust(16), str(points[i][1]+yShift).ljust(16), str(pt[3]).ljust(3), str(pt[4]).ljust(12))
        n = 0
        #construct the neighbours list
        while n < numNeighs:                   
            nline = "%s" % (str(pt[5+n])).ljust(8)
            line = line + nline
            n += 1
        
        line = line + "\n"
        ngh_fileOut.write(line) 
        i += 1

    ngh_fileOut.close()
    
    
def addPolygonPG(ply_in="../output/SQL/KBBuildingBoundingPoly2.ply",dbname="kingstonbeach",user="tbone",tablename="boundingpoly",epsg = "28355"):
    
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
    
    ply_file = open(ply_in, "r").readlines() 
    boundingPoly = []
    for row in ply_file:
        row =  row.split()
        boundingPoly.append(row)    
        
    line = "DROP TABLE " + tablename + ";"
    cur.execute(line)
    line = "CREATE TABLE " + tablename + " (ID int4);"
    cur.execute(line)
    line = "SELECT AddGeometryColumn( '" + tablename + "', 'polygon', " + epsg + ", 'POLYGON', 2);"
    cur.execute(line)
    
    polygonSQL = "ST_GeomFromText('POLYGON(("
       
    numPts = len(boundingPoly)
    i = 0
    while i < (numPts - 1):
        ptText = "%s %s, " % (str(boundingPoly[i][0]),str(boundingPoly[i][1]))
        polygonSQL = polygonSQL + ptText
        i += 1
        
    
    ptText = "%s %s))',%s));" % (str(boundingPoly[i][0]),str(boundingPoly[i][1]),epsg)
    polygonSQL = polygonSQL + ptText

    ptText = "INSERT INTO " + tablename + " (ID,polygon) VALUES (1,"
    polygonSQL = ptText + polygonSQL
    cur.execute(polygonSQL)    
        
        
    conn.commit()
    cur.close()
    conn.close()   
    
def addNGHRTR_to_PG(ngh_in,rtr_in=[], epsgIN=28355,dbname="kingstonbeach",user="tbone",tablename="ngh",tablenameTri="tri", shiftX=0.0, shiftY=0.0):


    import psycopg2
    import sys
    
    #open the conection to the requested database     
    line = "dbname='%s' user='%s'" % (dbname,user)
    try:
        conn = psycopg2.connect(line);
    except:
        print "Unable to connect to the database" 
        sys.exit()   
    
    cur = conn.cursor()


    ngh_fileIn = open(ngh_in, "r").readlines() 

    numNodes = int(ngh_fileIn[0].split()[0])
    numNeighs = int(ngh_fileIn[1].split()[0])
    ngh_fileIn.pop(0)
    ngh_fileIn.pop(0)
    ngh_fileIn.pop(0)

    points = []
    for vert in ngh_fileIn:
        vert = vert.split()
        points.append([float(vert[1]),float(vert[2]),float(vert[4])])

    #DROP table=tablename already exists
    cur.execute("DROP TABLE IF EXISTS %s;" % (tablename))
  
  
    #CREATE table=tablename
    #-----------------------------------  
    line = "CREATE TABLE %s (node_number int4, x float, y float, z float, node_code int4, number_of_neighbours int4"  % (tablename) 
    n = 0
    while n < numNeighs:                   
        nline = ", n%s int4" % (str(n+1))
        line = line + nline
        n += 1      
    line = line + ");"
    cur.execute(line)
    
    #add a POINT geometry column to hold define the [x,y] point
    cur.execute("SELECT AddGeometryColumn('%s', 'points', %s, 'POINT', 2);" % (tablename,epsgIN))    
    points = []
   
    print "Adding points to PostGIS table..."

    for pt in ngh_fileIn:
        
        pt = pt.split()       
        #add NGH file data into the PostGIS table
        line = "INSERT INTO %s (node_number, x, y, z, node_code, number_of_neighbours" % (tablename)
        n = 0
        while n < numNeighs:                   
            nline = ", n%s" % (str(n+1))
            line = line + nline
            n += 1 
        x = float(pt[1]) + shiftX
        y = float(pt[2]) + shiftY
        
        points.append([x,y])
        
        line = line +", points) VALUES (%s, %s, %s, %s, %s, %s" % (int(pt[0]),x,y,float(pt[4]),int(pt[3]),numNeighs)
        n = 0
        while n < numNeighs:                   
            nline = ", %s" % (int(pt[n+5]))
            line = line + nline
            n += 1 
            
        line = line + ", ST_GeomFromText('POINT(%s %s)',%s));" % (x, y, epsgIN)
        cur.execute(line)
    
    conn.commit()

    print "Building points index..."
    cur.execute("CREATE INDEX %s_idx ON %s USING GIST (points);" % (tablename,tablename))

    conn.commit()


    if(rtr_in != []):
        #Write the triangles
        #--------------------------------------------------------------------------------
        print "Adding triangles to PostGIS table..."
    
        rtr_file = open(rtr_in, "r").readlines() 
        numTriangles = len(rtr_file)
    
        triangles = []
        for tri in rtr_file:
            tri = tri.split()
            triangles.append([int(tri[0]),int(tri[1]),int(tri[2]), int(tri[4])])
    
        #DROP table=tablename already exists
        cur.execute("DROP TABLE IF EXISTS %s;" % (tablenameTri))
        cur.execute("CREATE TABLE %s (tri_number int4, v1 int4, v2 int4, v3 int4, tri_code int4);"  % (tablenameTri))    
        cur.execute("SELECT AddGeometryColumn('%s', 'triangles', %s, 'POLYGON', 2);" % (tablenameTri,epsgIN))    
       
        n = 1
        for tri in triangles:
            v1 = tri[0] - 1
            v2 = tri[1] - 1
            v3 = tri[2] - 1
            polySQL = "ST_GeomFromText('POLYGON((%s %s,%s %s,%s %s,%s %s))',%s));" % (points[v1][0],points[v1][1],points[v2][0], \
                                                                                      points[v2][1],points[v3][0],points[v3][1], \
                                                                                      points[v1][0],points[v1][1],epsgIN)
            line = "INSERT INTO %s (tri_number, v1, v2, v3, tri_code, triangles) VALUES (%s, %s , %s, %s, %s, " % (tablenameTri, n, v1+1, v2+1, v3+1, tri[3])
            line = line + polySQL
            cur.execute(line)  
            n += 1      
        
        
        conn.commit()
        print "Building triangles index..."
        cur.execute("CREATE INDEX %s_idx ON %s USING GIST (triangles);" % (tablenameTri,tablenameTri))
        conn.commit()

    
    cur.close()
    conn.close()
    
     
def offsetPOLYGON(polyList=[],offset=-0.1,dbname="kingstonbeach",user="tbone",tablename="footprints",geometry="the_geom",epsg = "28355" ):
    #Assumes that the input polygon is closed (i.e. the endpoint is the same as the start)
    #Assumes that input polygon is Clockwise defines (is this the standard??)
    #polyIn is a list of x,y points
    import math
  
    import psycopg2
    import sys
    
    from osgeo import ogr
    from osgeo import osr 
    
    from ricomGeospatialUtilities import addPolygonListPG
    
    
    #open the conection to the requested database     
    line = "dbname='%s' user='%s'" % (dbname,user)
    try:
        conn = psycopg2.connect(line);
    except:
        print "Unable to connect to the database" 
        sys.exit()   
    
    cur = conn.cursor()
  
    cur.execute("SELECT ST_AsText(%s) FROM %s;" % (geometry,tablename) )
    polySQL = cur.fetchall()
    polyList = []
    for poly in polySQL:
        #extract the polygon verices from postgis polygon    
        poly = poly[0]
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
    
    cam = len(polyList)
    
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
                
            #translate vector to the orgin (the point at index = n)
            d1 = [(poly[i0][0] - poly[i1][0]),(poly[i0][1] - poly[i1][1])]
            d2 = [(poly[i2][0] - poly[i1][0]),(poly[i2][1] - poly[i1][1])]
            #make all vectors unit vectors
            m1 = math.sqrt(d1[0]*d1[0]+d1[1]*d1[1])
            m2 = math.sqrt(d2[0]*d2[0]+d2[1]*d2[1])
            u1 = [d1[0]/m1, d1[1]/m1]
            u2 = [d2[0]/m2, d2[1]/m2]
            
            a = (math.atan2(d2[1], d2[0])-math.atan2(d1[1], d1[0]))
            if a > 0:
                a = (-2*math.pi+a)/2
            else:
                a = a/2
            
            R = [[math.cos(a), -math.sin(a)],[math.sin(a),math.cos(a)]]
            nx = u1[0]*R[0][0] + u1[1]*R[0][1]
            ny = u1[0]*R[1][0] + u1[1]*R[1][1]
    
            polyOffset.append([poly[n][0]+nx*offset,poly[n][1]+ny*offset])            
            n+=1
            
        polyOffsetList.append(polyOffset)

    
    # set the spatial reference - of the polygon data
    coordsLL = osr.SpatialReference()
    coordsLL.ImportFromEPSG(4326)
    coordsUTM = osr.SpatialReference()
    coordsUTM.ImportFromEPSG(int(epsg))

    #write the polygon data to KML
    driver = ogr.GetDriverByName('KML')
    kmlData = driver.CreateDataSource("../output/SQL/KingstonBeach5/FootprintOffset.kml")
    layer =  kmlData.CreateLayer("FootprintOffset", coordsLL, ogr.wkbLineString)
    feature = ogr.Feature(layer.GetLayerDefn())

    for poly in polyOffsetList:
        line = ogr.Geometry(type=ogr.wkbLineString)
        for pt in poly:
            line.AddPoint(pt[0],pt[1])
              
        line.AssignSpatialReference(coordsUTM)
        line.TransformTo(coordsLL)
        feature.SetGeometry(line)
        layer.CreateFeature(feature)
    
    for poly in polyList:
        line = ogr.Geometry(type=ogr.wkbLineString)
        for pt in poly:
            line.AddPoint(pt[0],pt[1])
              
        line.AssignSpatialReference(coordsUTM)
        line.TransformTo(coordsLL)
        feature.SetGeometry(line)
        layer.CreateFeature(feature)
    
    poly_tablename = tablename + "_offset"    
    addPolygonListPG(polyOffsetList,dbname=dbname,user=user,tablename=poly_tablename,epsg = epsg)
    
    kmlData.Destroy()
    feature.Destroy()

    cur.close()
    conn.close()
    
def test():
    import psycopg2
    import sys
    from osgeo import ogr
    from osgeo import osr 
    
    from ricomGeospatialUtilities import pg_geomList2Kml
    
    dbname = "kingstonbeach"
    user = "tbone"

    #open the conection to the requested database     
    line = "dbname='%s' user='%s'" % (dbname,user)
    try:
        conn = psycopg2.connect(line);
    except:
        print "Unable to connect to the database" 
        sys.exit()   
    
    cur = conn.cursor()
    #cur.execute("SELECT f.id, t.tri_number, ST_AsText(t.triangles) FROM tri_holes AS t, footprints_offset AS f WHERE ST_Intersects(t.triangles,f.polygon) ORDER BY f.id;" )
    #cur.execute("SELECT f.id, n.node_number, ST_AsText(n.points) FROM ngh_holes AS n, footprints_offset AS f WHERE ST_Contains(f.polygon, n.points) ORDER BY f.id;" )
    
    cur.execute("SELECT f.id, t.tri_number, ST_AsText(t.triangles) FROM tri AS t, footprints_offset AS f, boundingpoly AS p \
                WHERE (ST_Intersects(f.polygon,t.triangles) AND ST_Contains(p.polygon,f.polygon)) ORDER BY f.id;" )
    #cur.execute("SELECT f.id, t.tri_number, ST_AsText(t.triangles) FROM tri AS t, footprints_offset AS f WHERE ST_Contains(f.polygon, t.triangles) ORDER BY f.id;" )

    #cur.execute("SELECT f.id, n.node_number, ST_AsText(n.points) FROM ngh AS n, footprints_offset AS f WHERE ST_Contains(f.polygon, n.points) ORDER BY f.id;" )
    
    pointsReturned = cur.fetchall()
    
    geomList = []
    for row in pointsReturned:
        geomList.append(row[2])
    
#    n = 0
#    print "\nShow me the result:\n" 
#    for row in pointsReturned:
#        print "%s, %s, %s" % (str(row[0]), str(row[1]),str(row[2]))
#        n += 1
#    conn.commit()

    print "%s" % str(len(pointsReturned))
    
    cur.close()
    conn.close()


    pg_geomList2Kml (geomList=geomList,kml_file="../output/SQL/FootprintTri.kml")

    
def convert_points(pointsIN,epsgIN,epsgOUT):
       
    """
    Given a set of points ([x,y,z] or [x,y]) convert from espgIN to self.epsg
    
    """
    
    from osgeo import ogr #@UnresolvedImport
    from osgeo import osr #@UnresolvedImport 
    
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
    


        