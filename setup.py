#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# setup.py to install spacepy

__version__ = "$Revision: 1.2 $, $Date: 2010/05/20 17:36:48 $"
__log__ = """
$Log: setup.py,v $
Revision 1.2  2010/05/20 17:36:48  smorley
Updated setup.py for install of oneralib

Revision 1.12  2010/05/14 21:02:27  jkoller
removed conflict on package_data line

Revision 1.11  2010/04/28 20:55:26  smorley
Updates to radbelt documentation, spacetime, __init__ and setup

Revision 1.10  2010/04/14 18:56:24  jkoller
fixed libtool and f2py compiler flag required for MacOS

Revision 1.9  2010/03/11 22:40:12  smorley
Updates to Toolbox, __init__ and updated documentation.
Dan Welling's first function in SpacePy -- Now Available!
(While stocks last)

Revision 1.8  2010/03/11 20:21:14  jkoller
fixed entry

"""
__author__ = 'Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)'



# -------------------------------------
def compile_oneralib():
	
	# compile oneralib
	import os, sys
	
	
	os.chdir('spacepy/oneralib')
	
	F90files = ['source/onera_desp_lib.f', 'source/CoordTrans.f']
	functions = ['make_lstar1', 'make_lstar_shell_splitting1', \
		'coord_trans1 find_magequator1', 'find_mirror_point1', 
		'get_field1']
	
	# call f2py
	os.system('f2py --overwrite-signature -m onerapylib -h onerapylib.pyf '+' '.join(F90files) \
			+' only: ' + ' '.join(functions) + ' :')
	
	# intent(out) substitute list
	outlist = ['lm', 'lstar', 'blocal', 'bmin', 'xj', 'mlt', 'xout', 'bmin', 'posit', \
			'xgeo', 'bmir', 'bl', 'bxgeo']
	
	inlist = ['sysaxesin', 'sysaxesout', 'iyr', 'idoy', 'secs', 'xin']
	
	fln = 'onerapylib.pyf'
	
	print 'Substituting fortran intent(in/out) statements'
	f = open(fln, 'r')
	filestr = f.read()
	f.close()
        
	for item in inlist:
		filestr = subst( ':: '+item, ', intent(in) :: '+item, filestr)
		
	for item in outlist:
		filestr = subst( ':: '+item, ', intent(out) :: '+item, filestr)
	
	f = open(fln, 'w')
	f.write(filestr)
	f.close()

	
	# compile (platform dependent)
	os.chdir('source')
	if sys.platform == 'darwin': # then mac OS
   		os.system('gfortran -c -w -O2 -fPIC *.f')
   		os.system('libtool -static -o libBL2.a *.o')
   		os.chdir('..')
   		os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
   		
	elif sys.platform == 'linux2': # then linux   	
		os.system('gfortran -c -w -O2 -fPIC *.f')
		os.system('ar -r libBL2.a *.o')
		os.system('ranlib libBL2.a')
		os.chdir('..')
		os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
	else:
		print sys.platform, ' not supported at this time'
		sys.exit(1)
	
	os.system('mv -f onerapylib.so ../')
	os.chdir('../..')
	
	return

# -------------------------------------
def subst(pattern, replacement, filestr,
          pattern_matching_modifiers=None):
          
    """
    replace pattern by replacement in file
    pattern_matching_modifiers: re.DOTALL, re.MULTILINE, etc.
    """
    
    import os, re, sys, shutil
    
    
    if pattern_matching_modifiers is not None:
    	cp = re.compile(pattern, pattern_matching_modifiers)
    else:
    	cp = re.compile(pattern)
	
	if cp.search(filestr):  # any occurence of pattern?
		filestr = cp.sub(replacement, filestr)
		
    return filestr
	
# -------------------------------------
from distutils.core import setup
import os, sys, glob, re

# run compile for onera_desp_lib first
compile_oneralib()
pkg_files = ['oneralib/onerapylib.so','data/omnidata.pbin', 'data/tai-utc.dat',
         'data/PSDdb.pbin', 'doc/*.*']

# run setup from distutil
setup(name='spacepy',
      version='0.1',
      description='SpacePy: Tools for Space Science Applications',
      author='Steve Morley, Josef Koller, Dan Welling, Brian Larsen',
      author_email='spacepy_mail@lanl.gov',
      url='http://spacepy.lanl.gov',
      packages=['spacepy','spacepy.sandbox'],
      package_data={'spacepy': pkg_files}
      ) 
