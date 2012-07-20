import sys
from optparse import OptionParser
import os

#from osgeo import ogr
#from osgeo import osr
from osgeo import gdal
from osgeo.gdalconst import *
import struct
from numpy import *
'''
def main():
	import scipy 
	import scipy.interpolate 
	
	# the two axes 
	x = scipy.linspace(-1.0,1.0,10) 
	y = x 
	
	# make some pretend data 
	gridy, gridx = scipy.meshgrid(x,y) 
	z = scipy.sin(gridx)*scipy.sin(gridy) 
	
	# create a spline interpolator 
	spl = scipy.interpolate.RectBivariateSpline(x,y,z) 
	
	# make some new axes to interpolate to 
	nx = scipy.linspace(-1.0,1.0,100) 
	ny = nx 
	
	# evaluate 
	nz = spl(nx, ny) 
	
	# with matplotlib, compare: 
	import pylab 
	pylab.matshow(z) 
	pylab.matshow(nz) 
'''

def main():

	#Using the optparse module to parse the commandline arguements (see manual for more information)

	parser = OptionParser(usage="usage: %prog -i input_file -o output_file", version="%prog 1.0")    
	
	parser.add_option("-r", "--inputraster",
					  dest="raster_in",
					  help="The name of the input raster file")
	parser.add_option("-n", "--inputngh",
					  dest="ngh_in",
					  help="The name of the input ngh file")
	parser.add_option("-o", "--outputngh",
					  dest="ngh_out",
					  #default='output.asc',
					  help="The Name of the output ngh file")
	parser.add_option("-z", "--zscale",
					  dest="z_scale",
					  #default='output.asc',
					  help="Scale factor to apply to the z values of the raster grid")	
	
	
	(opts, args) = parser.parse_args()
	dataset = gdal.Open(opts.raster_in, GA_ReadOnly)
	if dataset is None:
		print "No gdal dataset named: %s" % opts.raster_in
		sys.exit()


	#Multiply the z value by a scaler.  NOTE: Japan500m bathy depth is positive
	z_scale = float(opts.z_scale)
	
	#read the input NGH file
	ngh_fileIn = open(opts.ngh_in, "r").readlines() 
	ngh_fileOut = open(opts.ngh_out, "w")
	numNodes = int(ngh_fileIn[2].split()[0])
	numNeighs = int(ngh_fileIn[3].split()[0])

	
	header = []
	header.append(ngh_fileIn[0])
	header.append(ngh_fileIn[1])
	header.append(ngh_fileIn[2])
	header.append(ngh_fileIn[3])
	ngh_fileIn.pop(0)
	ngh_fileIn.pop(0)
	ngh_fileIn.pop(0)
	ngh_fileIn.pop(0)

	#extract the points from the input NGH file
	points = []
	for vert in ngh_fileIn:
		vert = vert.split()
		points.append([float(vert[1]),float(vert[2]),float(vert[4])])
	
	


    
	#Tutorial - get info about the dataset
	print 'Driver: ', dataset.GetDriver().ShortName,'/', dataset.GetDriver().LongName
	print 'Size is ',dataset.RasterXSize,'x',dataset.RasterYSize,'x',dataset.RasterCount
	print 'Projection is ',dataset.GetProjection()
	
	T = dataset.GetGeoTransform()
	if not T is None:
		print 'Top-Left = (',T[0], ',',T[3],')'
		print 'Pixel Size = (',T[1], ',',T[5],')'
		
	#Tutorial get the raster band 
	band = dataset.GetRasterBand(1)
	print 'Band Type=',gdal.GetDataTypeName(band.DataType)

	min = band.GetMinimum()
	max = band.GetMaximum()
	if min is None or max is None:
		(min,max) = band.ComputeRasterMinMax(1)
		print 'Min=%.3f, Max=%.3f' % (min,max)



	xlr =  T[0] + band.XSize*T[1] + band.YSize*T[2] 
	ylr =  T[3] + band.XSize*T[4] + band.YSize*T[5]
	
	xul =  T[0]
	yul =  T[3]
	
	print "lower right x = %s, y = %s" % (xlr, ylr)

			
	#create list for each of the x coords (P - pixels), y coords (L - lines) and z values (elevation)
	
	x = arange(0,dataset.RasterXSize,1,dtype=int)
	y = arange(0,dataset.RasterYSize,1,dtype=int)

	iL = 0
	iP = 0
	z = []
	while iL < dataset.RasterYSize:
		#read one line at a time and insert into the z array
		scanline = band.ReadRaster( 0, iL, dataset.RasterXSize, 1,  dataset.RasterXSize, 1, GDT_Float32 )
		tuple_of_floats = struct.unpack('f' * (dataset.RasterXSize), scanline)
		z.append(tuple_of_floats)
		
		iL += 1
		#print iL

	z = array(z).transpose()
	import scipy 
	import scipy.interpolate 

	# create a spline interpolator 
	#spl = scipy.interpolate.RectBivariateSpline(x,y,z,kx = 2, ky = 2) 
	spl = scipy.interpolate.RectBivariateSpline(x,y,z,kx = 2, ky = 2) 

	#f = scipy.interpolate.interp2d(x,y,z, kind='linear') 
	
	num_pts = len(points)
	n = 0
	i = 0
	while n < num_pts:
		pt = points[n]
		#print "PT = %s" % pt
		
		#if the point is inside the bounds of the raster - interpolate the z from the raster
		if (pt[0] >= xul and pt[0] <= xlr) and (pt[1] >= ylr and pt[1] <= yul):
		
			#transformation from pixel coords to geo-coords (P = pixel, L = line, Xp = x coord of pt, Yp = ycoord of pt)
			#Xp = T[0] + P*T[1] + L*T[2];
			#Yp = T[3] + P*T[4] + L*T[5];
			L = (T[4]*(pt[0] - T[0]) + T[1]*(T[3] - pt[1])) / (T[3]*T[4] - T[5]*T[1])
			P = (T[5]*(pt[0] - T[0]) + T[2]*(T[3] - pt[1])) / (T[1]*T[5] - T[4]*T[2])
		

		
			#L = int(L)
			#P = int(P)
	
			#print "L = %s, P = %s" % ((L), (P))
			
			
			#old interpolation is simply nearest neighbour point.
			scanline = band.ReadRaster( int(P), int(L),1, 1, 1, 1, GDT_Float32 )
			tuple_of_floats = struct.unpack('f' * 1, scanline)
			
			#find the depth of the nearest neighbour
			zNN = float(tuple_of_floats[0])
			
			
			if zNN != 0:	#interpolate as long as the nearest neighbour is not NoDATA
				#new interpolation method uses a spline fit method
				zInterp = spl(P,L)[0][0]
				#zInterp = f(P,L)[0][0]
				#print zInterp
				#if zInterp < 0.0:
				points[n][2] = zInterp*z_scale
				i += 1
		#If outside - set the point elevation to 0
		else:
			points[n][2] = 0.0
		n +=1


	print "Number of points inside = %s" % i
	
	#Write the new z values to the output NGH file
	#----------------------------------------------------------------------------------------------
	#write the ngh file header
	ngh_fileOut.write(header[0])	
	ngh_fileOut.write(header[1])	
	ngh_fileOut.write(header[2])	
	ngh_fileOut.write(header[3])	

	i = 0
	#write vertices
	for pt in ngh_fileIn:
		pt = pt.split()
		code = int(pt[3])
		if (code != 0) and (code != 1) and (code != 2) and (code != 7):            #remove bad codes (Trigrid seems to add weird numbers)
			pt[3] = '0' 
		
		line = "%s %s %s %s %s   " % (str(pt[0]).ljust(8), str(points[i][0]).ljust(16), str(points[i][1]).ljust(16), str(pt[3]).ljust(3), str(points[i][2]).ljust(12))
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

if __name__ == '__main__':
    main()