#! /usr/bin/python
#
# Trajectory.py 	
#
# NOTES: 9/15/2005 - fixed a bug in rep vs. realrep trajectory writing...

from Config import *
from Chain import *
from Monty import *
from Replica import *

import random
import string
import math
import sys
import os


#################################################

class Trajectory:
    """A set of functions for creating, reading, writing, and organizing trajcetory files"""

    def __init__(self,replicas,config):
        """Initialize the trajectory object"""    
    
	self.trjqueue = []	# A queue of SWAPEVERY coordinates
        self.enequeue = []	# A queue of SWAPEVERY ene fields
	
        self.BUFFERSIZE = 500	# number of points to keep in queues
				# before dumping to file
	
        for i in range(0,len(replicas)):
	    self.trjqueue.append([])
	    self.enequeue.append([])
	    	    
	
        self.repfiles_trj = []		# file handles for coord trajs for each replica
        self.repfiles_ene = []		# file handles for ene trajs for each replica
	
        self.bytemp_trj = []		# file handles for coord trajs arranged by temp
        self.bytemp_ene = []		# file handles for ene trajs arranged by temp
        self.bytemp_replica = []	# file handle - keeps track of replica number
	
        self.byreplica_trj = []		# file handles for coord trajs arranged by temp
        self.byreplica_ene = []		# file handles for ene trajs arranged by temp
        self.byreplica_temp = []	# file handle - keeps track of temp number

        self.mkdir(config.EXPDIR)
	    
        # make /data /anal and /setup if they don't exist
	self.mkdir(config.DATADIR)
	self.mkdir(config.ANALDIR)
	self.mkdir(config.SETUPDIR)
	
        # make /workspace, /by_temp, /by_replica
        # Open trajectory files for writing in the worksapce
	self.mkdir(os.path.join(config.DATADIR, 'workspace') )
	self.mkdir(os.path.join(config.DATADIR, 'by_temp') )
	self.mkdir(os.path.join(config.DATADIR, 'by_replica') )
	
	
        # Open trajectory files for writing in the worksapce
	for i in range(0,len(replicas)):
	
	    workspacedir = os.path.join(config.DATADIR, 'workspace', str(replicas[i].repnum))
	    self.mkdir(workspacedir)
	    
	    trjname = os.path.join(workspacedir, 'mc.trj')
	    enename = os.path.join(workspacedir, 'mc.ene')

	    # fout = open(trjname,'w')
	    self.repfiles_trj.append(open(trjname,'w'))
	    # fout = open(enename,'w')
	    self.repfiles_ene.append(open(enename,'w'))

            ### write a header file explaining what the ene fields are
	    eneheader = os.path.join(workspacedir, 'header.ene')
	    self.write_eneheader(eneheader,replicas[i])
	    
	
        ### Open trajectory files for writing in the by_temp directory
	for i in range(0,len(replicas)):
	
	    bytempdir = os.path.join(config.DATADIR, 'by_temp' )
	    self.mkdir(bytempdir)
	    
	    trjname = bytempdir + '/' + str(replicas[i].repnum) + '.trj'
	    enename = bytempdir + '/' + str(replicas[i].repnum) + '.ene'
	    replicaname = bytempdir + '/' + str(replicas[i].repnum) + '.replica'
	    
	    self.bytemp_trj.append(open(trjname,'w'))
	    self.bytemp_ene.append(open(enename,'w'))
	    self.bytemp_replica.append(open(replicaname,'w'))
	    
	    ### write a header file explaining what the ene fields are
	    eneheader = bytempdir + '/header.ene'
	    self.write_eneheader(eneheader,replicas[i])

	    
	
        ### Open trajectory files for writing in the by_replica directory
	for i in range(0,len(replicas)):
	
	    byreplicadir = os.path.join(config.DATADIR, 'by_replica')
	    self.mkdir(byreplicadir)
	    
	    trjname = byreplicadir + '/' + str(replicas[i].repnum) + '.trj'
	    enename = byreplicadir + '/' + str(replicas[i].repnum) + '.ene'
	    tempname = byreplicadir + '/' + str(replicas[i].repnum) + '.temp'
	    
	    self.byreplica_trj.append(open(trjname,'w'))
	    self.byreplica_ene.append(open(enename,'w'))
	    self.byreplica_temp.append(open(tempname,'w'))

	    ### write a header file explaining what the ene fields are
	    eneheader = byreplicadir + '/header.ene'
	    self.write_eneheader(eneheader,replicas[i])
	    	
        ### print 'REpfiles:', self.repfiles_trj



    def mkdir(self,pathname):
	"""Automatically create directory if it doesn't exist"""
        if os.path.exists(pathname)==False:		
            os.mkdir(pathname)

    def write_eneheader(self,filename,replica):
        """Write column headers for the energy file."""
 
        fheader = open(filename,'w')
        fheader.write('E_pot\tE_rest(D)\tD\tcontact_state\ttemp\n')
        fheader.write('# Energy units: Joules/mol\n')
        fheader.write('# Restrained contact state: ' + repr(replica.mc.restraint.contacts) + '\n')
        fheader.write('# kspring: '+str(replica.mc.restraint.kspring) + '\n')
	fheader.close()

    def queue_trj(self,replica):
        """Queue a trajectory point to the buffer for writing to file"""
	point = []
	for i in range(0, len(replica.chain.coords) ):
	    point.append((replica.chain.coords[i][0],replica.chain.coords[i][1]))
	self.trjqueue[replica.repnum].append(point)
	if len(self.trjqueue[replica.repnum]) >= self.BUFFERSIZE:
	    self.dump_trjqueue(replica)
	    


    def queue_ene(self,replica):
        """Queue an energy value to the buffer, for writing to file."""
	self.enequeue[replica.repnum].append([ replica.mc.energy(replica.chain), replica.mc.restraint.energy(replica.chain), replica.mc.restraint.D(replica.chain),replica.chain.contactstate()])

	if len(self.enequeue[replica.repnum]) == self.BUFFERSIZE:
	    self.dump_enequeue(replica)

    def dump_trjqueue(self,replica):
        """Dump the queue to the the respective files and clear them for future use."""
	
	# write coords and enes to the workspace, by_temp and by_replica
        rep = replica.repnum
	
	### WORKSPACE FILES ###
	for pt in range(0,len(self.trjqueue[rep])):

	    self.repfiles_trj[rep].write(repr(self.trjqueue[rep][pt]))
	    self.repfiles_trj[rep].write('\n')

	### BY_TEMP and BY_REPLICA FILES ###
        realrep = replica.mc.tempfromrep
	
	self.byreplica_temp[rep].write(str(rep))
	self.byreplica_temp[rep].write('\n')

	self.bytemp_replica[rep].write(str(realrep))
	self.bytemp_replica[rep].write('\n')

	for pt in range(0,len(self.trjqueue[realrep])):
	    self.bytemp_trj[rep].write(repr(self.trjqueue[realrep][pt]))
	    self.bytemp_trj[rep].write('\n')

	for pt in range(0,len(self.trjqueue[rep])):
	    self.byreplica_trj[rep].write(repr(self.trjqueue[rep][pt]))
	    self.byreplica_trj[rep].write('\n')

        ### clear the trj and ene queues
        self.trjqueue[rep] = []

	
    def dump_enequeue(self,replica):
        """Dumps the queued energy values to the respective files and clears them for further use."""

	# write coords and enes to the workspace, by_temp and by_replica
        rep = replica.repnum

	### WORKSPACE FILES ###
	for pt in range(0,len(self.enequeue[rep])):
            ## tab-delimit the ene fields
	    for field in self.enequeue[rep][pt]:
	    	self.repfiles_ene[rep].write(str(field) + '\t')
	    self.repfiles_ene[rep].write(str(replica.mc.temp) + '\n')

	### BY_TEMP and BY_REPLICA FILES ###
        realrep = replica.mc.tempfromrep
	
	for pt in range(0,len(self.enequeue[rep])):
            ## tab-delimit the ene fields
	    for field in self.enequeue[rep][pt]:
	    	self.bytemp_ene[realrep].write(str(field) + '\t')
	    self.bytemp_ene[realrep].write(str(replica.mc.temp) + '\n')

            ## tab-delimit the ene fields
	    for field in self.enequeue[rep][pt]:
	    	self.byreplica_ene[rep].write(str(field) + '\t')
	    self.byreplica_ene[rep].write(str(replica.mc.temp) + '\n')

        ### clear the trj and ene queues
        self.enequeue[rep] = []


    def cleanup(self,replicas):
        """Write any remaining points in the trajectory and energy buffers to file,
        and close any open file handles."""
	
	# write the last of the trj and ene buffers 
        for rep in range(0,len(replicas)):	
            self.dump_trjqueue(replicas[rep])
            self.dump_enequeue(replicas[rep])

	# close any open file handles
        for i in range(0,len(self.repfiles_trj)):
	    self.repfiles_trj[i].close() 	# file handles for coord trajs for each replica
            self.repfiles_ene[i].close()	# file handles for ene trajs for each replica
	
            self.bytemp_trj[i].close()		# file handles for coord trajs arranged by temp
            self.bytemp_ene[i].close()		# file handles for ene trajs arranged by temp
            self.bytemp_replica[i].close()	# file handle - keeps track of replica number
	
            self.byreplica_trj[i].close()	# file handles for coord trajs arranged by temp
            self.byreplica_ene[i].close()	# file handles for ene trajs arranged by temp
            self.byreplica_temp[i].close()	# file handle - keeps track of temp number

	
