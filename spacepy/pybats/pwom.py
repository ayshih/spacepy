#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the Polar Wind Outflow Model
(PWOM), one of the choices for the PW module in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray

class Line(PbData):
    '''
    Class for loading a single field line output file.
    At instantiation time, user may wish to set the start date and time of 
    the simulation using the starttime kwarg.  If not given, start time
    will default to Jan. 1st, 2000, 00:00UT.
    '''
    def __init__(self, filename, starttime=None, *args, **kwargs):
        super(Line, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(starttime)

    def __repr__(self):
        return 'PWOM single field line output file %s' % (self.attrs['file'])

    def _read(self, starttime):
        '''
        Read ascii line file; should only be called upon instantiation.
        '''
        
        # Slurp whole file.
        f=open(self.attrs['file'], 'r')
        lines=f.readlines()
        f.close()

        # Determine size of file.
        nTimes=lines.count(lines[0])
        nAlts =int(lines[2].strip())
        self.attrs['nAlt']=nAlts; self.attrs['nTime']=nTimes

        # Start building time array.
        self['time']=np.zeros(nTimes, dtype=object)

        # Get variable names; pop radius (altitude).
        var=(lines[4].split())[1:-1]
        self._rawvar=var

        # Get altitude, which is constant at all times.
        self['r']=dmarray(np.zeros(nAlts),{'units':'km'})
        for i, l in enumerate(lines[5:nAlts+5]):
            self['r'][i]=float(l.split()[0])

        # Create 2D arrays for data that is time and alt. dependent.
        # Set as many units as possible.
        for v in var:
            self[v]=dmarray(np.zeros((nTimes, nAlts)))
            if v=='Lat' or v=='Lon':
                self[v].attrs['units']='deg'
            elif v[0]=='u':
                self[v].attrs['units']='km/s'
            elif v[0:3]=='lgn':
                self[v].attrs['units']='log(cm-3)'
            elif v[0]=='T':
                self[v].attrs['units']='K'
            else:
                self[v].attrs['units']=None

        # Loop through rest of data to fill arrays.
        for i in range(nTimes):
            t=float((lines[i*(nAlts+5)+1].split())[1])
            self['time'][i]=starttime+dt.timedelta(seconds=t)
            for j, l in enumerate(lines[i*(nAlts+5)+5:(i+1)*(nAlts+5)]):
                parts=l.split()
                for k,v in enumerate(var):
                    self[v][i,j]=float(parts[k+1])

class Lines(PbData):
    '''
    A class for reading and plotting a complete set of PWOM line output files.
    Open using a glob string that encompasses all of the lines that are
    intended to be read, e.g. 'north*.out'.  Note that unix wildcards are
    accepted.
    '''

    def __init__(self, lines, starttime=None, *args, **kwargs):
        super(Lines, self).__init__(*args, **kwargs) # Init as PbData.
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(lines)

    def __repr__(self):
        return 'PWOM field line group of %i separate line files.' %\
            (self.attrs('nFile'))

    def _read(self, lines):
        '''
        Read all ascii line files; should only be called upon instantiation.
        '''

        from glob import glob
        
        self['files']=glob(lines)
        
        # Read first file; use it to set up all arrays.
        l1=Line(self['files'][0])
        nA=l1.attrs['nAlt']
        nT=l1.attrs['nTime']
        nF=len(self['files'])
        
        self.attrs['nAlt']=nA; self.attrs['nTime']=nT; self.attrs['nFile']=nF
        self['time']=l1['time']

        # Create all arrays.
        for v in l1._rawvar:
            self[v]=dmarray(np.zeros((nF, nT, nA)))
            self[v][0,:,:]=l1[v]
            if v=='Lat' or v=='Lon':
                self[v].attrs['units']='deg'
            elif v[0]=='u':
                self[v].attrs['units']='km/s'
            elif v[0:3]=='lgn':
                self[v].attrs['units']='log(cm-3)'
            elif v[0]=='T':
                self[v].attrs['units']='K'
            else:
                self[v].attrs['units']=None

        # Place data into arrays.
        for i, f in enumerate(self['files'][1:]):
            l=Line(f)
            for v in l1._rawvar:
                self[v][i+1,:,:]=l[v]
