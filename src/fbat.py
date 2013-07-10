
# FBAT - ISB implementation.
# david l gibbs
# dgibbs@systemsbiology.org
# June 20, 2013


# Goals:
# single marker fbat test
# rare-variant region-based test
#
 
# Input:
# tped and tfam files as specified by plink.  No recoding needed.
# for now: trios only.
# and additive coding for X .. i.e. the number of A alleles  
#    also expected to be biallelic
 
# Output:
# tab separated table of test results, p-values


# the test statistic is defined as:
# U = sum_k U_k
# where U_k = sum_ij (Y_k - mu)(X_ijk - E(X_ijk) = T_ij * X_ijk - E(X_ijk) * T_ij
# k indexes the pedigrees
# j is the offspring of the ith nuclear family
#
# Var(U) = sum_i (U_i)^2 for one nuclear family per pedigree
# where U_i is the i-th families contribution to U

import sys
import math
import singleTest
import rareTest

class FBAT:
    """ the fbat test object """
    
    def __init__(self):
        """ any initial tasks """
        self.tfamfile = ""
        self.tpedfile = ""
        self.freqcutoff = 0
        self.tfam = []   # list of the family info, first six columns of ped
        self.tped = []   # list of the snps .. each individual has two alleles [A11, A21, A12, A22, A13, A23, ...]
        self.X = []      # the X matrix
        self.resultTable = [] # save the results table ... maybe not?
        self.minFamilies = 10 # the number of informative familes
        self.famidx = []      # index into pedigrees found in tfam
        self.childidx = []      # index into children
        self.paridx = []        # index of the parents
        self.markers = []     # a list of the markers
        self.pedPhenotype = []  # phenotypes just like the ped
        self.phenotype = []   # the phenotypes
        self.offset = 0       # T_ij = Y_ij - offset
        
    def load (self, tfam_file, tped_file):
        """ load the data here """
        self.tpedfile=tped_file
        self.tfamfile=tfam_file
        self.tped = open(tped_file,'r')
        tfam = open(tfam_file,'r').read().strip().split("\n")
        self.tfam = map(lambda x: x.split("\t"), tfam)
        if len(self.tfam[1]) != 6:
            print("Error: tfam has wrong number of columns")
            sys.exit(2)
        self.pedPhenotypes = map(lambda x: (x[5]), self.tfam)
        self.phenotypes = self.adjPhenotypes(self.pedPhenotypes)
        self.famidx, self.childidx, self.paridx = self.familyIndex(self.tfam)
        self.checkForTrios()
        self.printFBAT()

    def loadIndex (self, index_file):
        

    def checkForTrios(self):
        gap = []
        for i in range(len(self.famidx)-1):
            gap.append(self.famidx[i+1] - self.famidx[i])
        if any([x != 3 for x in gap]):
            print("\nWARNING: non-trio pedigrees detected!\n")
            print("check around family: ",
                  str([i for i in range(len(gap)) if gap[i] != 3]))
            sys.exit(0)

    def setFreqCutoff(self, f):
        self.freqcutoff = float(f)

    def setOffset(self, o):
        self.offset = o
        self.applyOffset()

    def adjPhenotypes(self,p):
        # fbat takes ped format: 2 == affected, 1 == unaffected
        # and codes them as 1 and 0 respectively.
        ps = []
        for pi in p:
            if pi == '2':
                ps.append('1')
            elif pi == '1':
                ps.append('0')
            else:
                ps.append('0')
        return(ps)
        
    def familyIndex(self, tfam):
            """ a list of families and indices into tfam """
            famid = [] # index into families
            chid = []  # index into children within families
            parid = []
            old = -1
            c = []
            p = []
            for i in range(len(self.phenotypes)):
                  if old != tfam[i][0]: # if we are starting a new family
                        if len(c) > 0:
                              chid.append(c) # if we were already making an index save it.
                              parid.append(p)
                        j = 0   # new child index starting for new family
                        c = []  # index of children
                        p = []  # index of parents
                        famid.append(i) # record the index
                        old = tfam[i][0]         # remember we're in family
                  if tfam[i][2] == '0' and tfam[i][3] == '0': # then its a parent
                        p.append(j)
                  else:  # it's a kid
                        c.append(j)
                  j += 1
                  if i == (len(self.phenotypes)-1): # end of the file
                        chid.append(c)
                        parid.append(p)
            return([famid,chid,parid])                  
    
    def applyOffset(self):
        ps = []
        for i in range(len(self.pedPhenotypes)):
            if self.pedPhenotypes[i] != '0':
                ps.append(float(self.phenotypes[i])-float(self.offset))
            else:
                ps.append(float(0))
        self.phenotypes = ps

    def single(self):
        """ perform the single marker test """
        print("Marker\tAllele\tafreq\tNfams\tS-E(S)\tVar(S)\tZ\tP")
        for thisg in self.tped:
            gs = thisg.strip().split("\t")
            marker = gs[1]
            s = singleTest.SingleTest(marker, gs[4:], self.phenotypes,
                                      self.famidx, self.childidx, self.paridx,
                                      False, self.freqcutoff)
            s.test(True)

    def rare(self, regionfile, freqfile, weighted):
        """ reads the region file and performs rare variant tests """
        print("Region\tW\tVarW\tZ\tpvalue\tallele_freqs\tweights")
        regions = open(regionfile,'r').read().strip().split("\n")
        regions = map(lambda x: x.strip(), regions)
        if freqfile != "none":
            freqs = open(freqfile,'r').read().strip().split("\n")
        for r in regions:
            regionname = r.split("\t")[0]
            theseMarkers = r.split("\t")[1:]
            # find them in the index ...
            # then make the list of data using the file offsets .. call it thisData
            # get the allele freqs from the freqfile ... if requested.
            #    the freq file needs to be 
            t = rareTest.RareTest(regionname, theseMarkers, freqs, thisData,
                                  self.phenotypes, self.famidx, self.childidx,
                                  self.paridx, weighted)
            t.test()
            
    def printFBAT(self):
        """ print a little statement about the object """
        print("#FBAT")
        print("#files: " + self.tfamfile + "  " + self.tpedfile)
        print("#Number of families: " + str(len(self.famidx)))
        print("#Offset: " + str(self.offset))
        
    def printResults(self):
        """ print the results to a tab separated table """


