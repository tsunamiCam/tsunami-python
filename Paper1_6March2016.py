#Import the required modules


import psycopg2 #@UnresolvedImport
#import errorcodes

from netCDF4 import Dataset #@UnresolvedImport
import time
import os
import sys
from osgeo import ogr #@UnresolvedImport
from osgeo import osr #@UnresolvedImport




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


class FitProbit:

    def __init__(self,dbname ):



        #connect to the POSTGIS database
        try:           #open connection to the requested database
            self.conn = psycopg2.connect(database=dbname, user="tbone");
        except Exception, e:
            if e[0].find("could not connect to server") >= 0:
                print "Postgis database server not running. Please start postgres\n"
                sys.exit()

        self.cur =self.conn.cursor()




    def __del__ (self):
        """
        Class deconstructor

        """
        self.cur.close()
        self.conn.close()


    def cdf_probit_ln(self,x,b0,b1):
        x = b1*log(x) + b0
        Prob = 0.5*(1 + erf(x/sqrt(2)))
        return Prob

    #The CDF for PROBIT
    def cdf_probit(self,x,b0,b1):

        x = b1*x + b0
        Prob = 0.5*(1 + erf(x/sqrt(2)))
        return Prob

    #Given the fitted PROBIT LN parameters (b0 and b1), solve for mu (i.e. the standard mean)
    def get_mu_from_cdf_probit_ln(self,b0,b1):

        #The standard mean of a CDF is were probability = 0.5
        mu = (sqrt(2)*erfinv(0) - b0)/b1
        x = 5
        sig = (log(x) - mu) / ( sqrt(2) * erfinv(erf((b1*log(x) + b0)/sqrt(2))) )

        return sig,mu


    def create_fragility_probit_ln_NEW(self,tablename="buildings_results_all_6march2016", column = "r1",where="",type = 1,log_depth=True):

        if type == 1:
            self.cur.execute("\
            SELECT  \
                id, \
                %s, \
                CASE WHEN damage IN (0,6,7) THEN 1 ELSE 0 END \
            FROM %s \
            WHERE \
               %s > 0" % (column,tablename,column))
        elif type == 2:
            self.cur.execute("\
            SELECT  \
                id, \
                %s, \
                CASE WHEN damage IN (6,7) THEN 1 ELSE 0 END \
            FROM %s \
            WHERE \
                damage NOT IN (0,8) AND %s > 0 AND (visibility_before in (1,2) or visibility_after in (1,2));" % (column,tablename, column))
        else:
            print "A valid type has not been provided.  Choose type = 1 to use all data in the database and type = 2 to use only the TIB data."
            return


        damaged_bldgs = self.cur.fetchall()

        depth = []
        phiinv = []
        histogram = []



        for b in damaged_bldgs:
            if b[1] != None and b[1] > 0:
                depth.append([b[1],b[2]])


        depth.sort()
        #print depth
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
            binary_damage.append(damage)
            i+=1

        import pysal
        from pysal.spreg.probit import Probit

        depths_ln = []
        for d in binary_depths:
            depths_ln.append(log(d))


        #x = np.array(binary_depths, ndmin=2).T

        x_ln = np.array(depths_ln, ndmin=2).T
        x = np.array(binary_depths, ndmin=2).T
        y = np.array(binary_damage, ndmin=2).T



        depths1 = []
        depths0 = []
        i=0
        while i < num_bldgs:
            if depth[i][1] ==  1:
                depths1.append(depth[i][0])
            else:
                depths0.append(depth[i][0])
            i+=1

        #print mean(np.asarray(depths1)),std(np.asarray(depths1))
        #print mean(np.asarray(depths0)),std(np.asarray(depths0))



        #log depth
        model = Probit(y,x_ln)
        betas = np.around(model.betas, decimals=6)
        b0_ln = betas[0][0]
        b1_ln = betas[1][0]

        #Don't log depths
        model = Probit(y,x)
        betas = np.around(model.betas, decimals=6)
        b0 = betas[0][0]
        b1 = betas[1][0]

        #print model.summary
        i=0

        damage_class_count = 0
        while i< len(binary_damage):


            if binary_damage[i] == 0:
                damage_class_count+=1
                binary_damage[i] = 0.02*1
            else:
                binary_damage[i] = 0.98*1
            i+=1


        if (log_depth == True):
            sig,mu = self.get_mu_from_cdf_probit_ln(b0_ln,b1_ln)
            print mu, sig, exp(mu), len(x_ln)
            return b0_ln,b1_ln,binary_depths,binary_damage

        else:
            return b0,b1,binary_depths,binary_damage


    def create_fragility_binned_ln(self,tablename="buildings_results_all_6march2016", column = "r1",where="",weight = True,log_depth=True,bin_size = 0):

        self.cur.execute("SELECT id,%s,D5 FROM %s" % (column,tablename))
        damaged_bldgs = self.cur.fetchall()

        #print "SUPPASRI: Num. of bldgs. with damage %s = %s" % (damage_class,len(damaged_bldgs))

        depth = []
        phiinv = []
        for b in damaged_bldgs:
            if b[1] != None and b[1] > 0:
                depth.append([b[1],b[2]])

        depth.sort()
        print len(depth)
        num_bldgs = len(depth)
        depth_class1 = []
        i = 0

        start = depth[0][0]
        end = depth[0][0] + 1.0

        p = []
        depth_bins = []
        depths_in_bin = []
        damage_in_bin = []
        #bin_size = 8
        count = 0
        count_class = 0

        #according to Porter 2007 - bin size should be SQRT(num_bldgs) rounded up to the nearest integer value

        if bin_size == 0:
            bin_size = int(sqrt(num_bldgs))
            if np.mod(float(sqrt(float(num_bldgs))),1) > 0:
                bin_size = bin_size + 1

        j = 0
        start = 0





        #Fill the bins with counts of [survived, destroyed] specimens
        histogram_bins = []

        #print "SUPPASRI: Bin Size = %s" % bin_size
        while i < num_bldgs:
            damage = depth[i][1]
            fd = depth[i][0]
            count += 1
            depths_in_bin.append(fd)
            damage_in_bin.append(damage)

            if damage == 1:
                count_class += 1

            if count == bin_size:

                j += 1

                prob = float(count_class)/float(count)
                start = i+1
                #if prob == 0 and j == 1:
                if prob == 0:
                #if prob == 0:

                    average_depth_bin = 0
                    for dep in depths_in_bin:
                        average_depth_bin = average_depth_bin + dep
                    average_depth_bin = average_depth_bin/len(depths_in_bin)
                    depth_bins.append(average_depth_bin)
                    p.append(0.01)
                    #p.append(1.0/float(bin_size+1))

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

                #[survived, destroyed]
                histogram_bins.append([bin_size-count_class, count_class])
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
                #p.append(0.01)
                p.append(1.0/float(bin_size+1))


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

            histogram_bins.append([count-count_class, count_class])

        def phiinv(p):
            y = sqrt(2)*erfinv(2*p - 1)
            return y

        def fragility(x,mu,sig):
            Prob = 0.5*(1 + erf((x - mu)/(sqrt(2)*sig)))
            return Prob

        pinv = []
        p100 = []
        depth_binsln = []
        for d in depth_bins:
            depth_binsln.append(log(d))

        for prob in p:
            pinv.append(phiinv(prob))
            p100.append(prob*100)


        from scipy.stats import linregress
        #r^2 = r_value^2

        '''
        if damage_class == 2:
            i = 0
            while i<20:
                pinv.pop()
                depth_binsln.pop()
                histogram_bins.pop()
                depth_bins.pop()
                p.pop()
                i+=1
        '''

        slope, intercept, r_value, p_value, std_err = linregress(pinv,depth_bins)
        slopeln, interceptln, r_valueln, p_valueln, std_errln = linregress(pinv,depth_binsln)

        sig = slope
        mu = intercept
        r2 = r_value * r_value
        #fit =  polyfit(pinv,depth_binsln,1,full=True)
        #sig,mu = fit[0]

        sig_ln = slopeln
        mu_ln = interceptln
        r2_ln = r_valueln * r_valueln
        #print "mu_ln = %s, sig_ln = %s, r2_ln = %s" % (mu_ln,sig_ln,r2_ln)
        #fit =  polyfit(pinv,depth_binsln,1,full=True)
        #sig,mu = fit[0]


        bin = depth_bins
        #prob = p100
        prob = p

        survived = []
        destroyed = []

        for b in histogram_bins:
            survived.append(b[0])
            destroyed.append(b[1])


        if (log_depth == True):
            #return sig_ln,mu_ln,bin,prob, survived, destroyed,bin_size
            return sig_ln,mu_ln,bin,prob, survived, destroyed,bin_size

        else:
            return sig,mu,bin,prob, survived, destroyed, r2
