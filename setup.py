#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# setup.py to install spacepy

__version__ = "$Revision: 1.38 $, $Date: 2010/11/16 23:58:48 $"
__author__ = 'The SpacePy Team, Los Alamos National Lab (spacepy@lanl.gov)'

import os, sys, shutil, getopt, warnings
from distutils.core import setup
from os import environ as ENVIRON

# -------------------------------------
def compile_pybats():
    os.chdir('spacepy/pybats')
    os.system('f2py -c ctrace2d.pyf trace2d.c')
    os.chdir('../..')
    
# -------------------------------------
def compile_irbempy(fcompiler):

    # compile irbemlib
    os.chdir('spacepy/irbempy/irbem-lib-2010-09-14-rev189')

    F90files = ['source/onera_desp_lib.f', 'source/CoordTrans.f', 'source/AE8_AP8.f']
    functions = ['make_lstar1', 'make_lstar_shell_splitting1', \
        'coord_trans1 find_magequator1', 'find_mirror_point1', 
        'get_field1', 'get_ae8_ap8_flux', 'fly_in_nasa_aeap1']

    # call f2py
    os.system('f2py --overwrite-signature -m irbempylib -h irbempylib.pyf '+' '.join(F90files) \
            +' only: ' + ' '.join(functions) + ' :')

    # intent(out) substitute list
    outlist = ['lm', 'lstar', 'blocal', 'bmin', 'xj', 'mlt', 'xout', 'bmin', 'posit', \
            'xgeo', 'bmir', 'bl', 'bxgeo', 'flux']

    inlist = ['sysaxesin', 'sysaxesout', 'iyr', 'idoy', 'secs', 'xin']

    fln = 'irbempylib.pyf'

    print('Substituting fortran intent(in/out) statements')
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
        if fcompiler == 'pg':
            raise NotImplementedError('Portland Group compiler option "pg" is not supported on Mac OS with f2py')
            #os.system('pgf77 -c -Mnosecond_underscore -w -fastsse *.f')
            #os.system('libtool -static -o libBL2.a *.o')
            #os.chdir('..')
            #os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=pg')
        elif fcompiler == 'gnu':
            os.system('g77 -c -w -O2 -fPIC -fno-second-underscore *.f')
            os.system('libtool -static -o libBL2.a *.o')
            os.chdir('..')
            os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu --f77flags=-fno-second-underscore ')
        else:
            os.system('gfortran -c -w -O2 -fPIC -ffixed-line-length-none *.f')
            os.system('libtool -static -o libBL2.a *.o')
            os.chdir('..')
            os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
        
    elif sys.platform == 'linux2': # then linux
        if fcompiler == 'pg':
            os.system('pgf77 -c -Mnosecond_underscore -w -fastsse -fPIC *.f')
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')
            os.chdir('..')
            os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=pg')
        elif fcompiler == 'gnu':           
            os.system('g77 -c -w -O2 -fPIC -fno-second-underscore *.f')
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')
            os.chdir('..')
            os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu --f77flags=-fno-second-underscore')
        else:
            os.system('gfortran -c -w -O2 -fPIC -ffixed-line-length-none *.f')
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')
            os.chdir('..')
            os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
    else:
        print('%s not supported at this time' % sys.platform)
        sys.exit(1)
        
    err = os.system('mv -f irbempylib.so ../')
    if err:
        print '------------------------------------------------------'
        print 'WARNING: Something went wrong with compiling irbemlib.'
        print '------------------------------------------------------'
    
    os.chdir('../../..')

    return


def compile_oneralib(fcompiler):

    # compile oneralib
    os.chdir('spacepy/onerapy/onera_lib_V4.1')

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

    print('Substituting fortran intent(in/out) statements')
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
        if fcompiler == 'pg':
            os.system('pgf77 -c -Mnosecond_underscore -w -fastsse *.f')
            os.system('libtool -static -o libBL2.a *.o')
            os.chdir('..')
            os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=pg')
        elif fcompiler == 'gnu':
            os.system('g77 -c -w -O2 -fPIC -fno-second-underscore *.f')
            os.system('libtool -static -o libBL2.a *.o')
            os.chdir('..')
            os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu --f77flags=-fno-second-underscore ')
        else:
            os.system('gfortran -c -w -O2 -fPIC *.f')
            os.system('libtool -static -o libBL2.a *.o')
            os.chdir('..')
            os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
        
    elif sys.platform == 'linux2': # then linux
        if fcompiler == 'pg':
            os.system('pgf77 -c -Mnosecond_underscore -w -fastsse -fPIC *.f')
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')
            os.chdir('..')
            os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=pg')
        elif fcompiler == 'gnu':
            os.system('g77 -c -w -O2 -fPIC -fno-second-underscore *.f')
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')
            os.chdir('..')
            os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu --f77flags=-fno-second-underscore')
        else:
            os.system('gfortran -c -w -O2 -fPIC *.f')
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')
            os.chdir('..')
            os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
    else:
        print('%s not supported at this time' % sys.platform)
        sys.exit(1)

    os.system('mv -f onerapylib.so ../')
    os.chdir('../../..')

    return

# -------------------------------------
def subst(pattern, replacement, filestr,
          pattern_matching_modifiers=None):
          
    """
    replace pattern by replacement in file
    pattern_matching_modifiers: re.DOTALL, re.MULTILINE, etc.
    """
    
    import re, shutil
    
    
    if pattern_matching_modifiers is not None:
        cp = re.compile(pattern, pattern_matching_modifiers)
    else:
        cp = re.compile(pattern)

    if cp.search(filestr):  # any occurence of pattern?
        filestr = cp.sub(replacement, filestr)
        
    return filestr

# -------------------------------------

# check for compiler flags
fcompiler = 'gnu95' # standard compiler flag
if len(sys.argv) > 2:
    for i in range(len(sys.argv)):
        if sys.argv[i] == '--fcompiler=pg':
            fcompiler = 'pg'
        if sys.argv[i] == '--fcompiler=gnu':
            fcompiler = 'gnu'
            
for ff in ['pg', 'gnu95', 'gnu']:
    try:
        sys.argv.remove('--fcompiler='+ff)
    except:
        pass

#import tooblox by reading file from repository
# this will provide mostly the query_yes_no function and update()
exec(compile(open('spacepy/toolbox.py').read(), 'spacepy/toolbox.py', 'exec'))

#test for python version 2.x where x>=5
try:
    dum = sys.version_info
    if dum[0]==2:
        assert dum[1]>=5
    import numpy
except:
    raise Exception("""SpacePy requires Python 2.X, where X>=5.\n
    Numpy, Scipy and Matplotlib(>=0.99) are also required\n
    Please install suitable versions.""")
try:
    import scipy
    import matplotlib
except:
    warnings.warn('''Missing Packages: SciPy and MatPlotLib are 
    required for large parts of this library.''')

# run compile for irbem-lib first
if os.path.exists('spacepy/irbempy/irbempylib.so'):
    ans = query_yes_no('\nDo you want to recompile the IRBEM library?', default="no")
    if ans=='yes':
        compile_irbempy(fcompiler)
else:
    compile_irbempy(fcompiler)
    
# run compile for onera_desp_lib first
#if os.path.exists('spacepy/onerapy/onerapylib.so'):
#    ans = query_yes_no('\nDo you want to recompile the ONERA-DESP library?', default="no")
#    if ans=='yes':
#        compile_oneralib(fcompiler)
#else:
#    compile_oneralib(fcompiler)


# Compile PyBats
compile_pybats()

# create .spacepy in $HOME and move data
# read-in .rc file first
exec(compile(open('spacepy/data/spacepy.rc').read(), 'spacepy/data/spacepy.rc', 'exec'))
if 'SPACEPY' in ENVIRON:
    DOT_FLN = ENVIRON['SPACEPY']+'/.spacepy'
else:
    DOT_FLN = ENVIRON['HOME']+'/.spacepy'

if os.path.exists(DOT_FLN):
    ans = query_yes_no('\n'+DOT_FLN+' already exists. Do you want to start fresh?', default="no")
    if ans=='no':
        fresh_install = False
    else:
        fresh_install = True
        i = 0
        while 1:
            if os.path.exists(DOT_FLN+'.bak.'+str(i)):
                i = i+1
            else:
                shutil.move(DOT_FLN, DOT_FLN+'.bak.'+str(i))
                break
else:
    fresh_install = True

if fresh_install:
    os.mkdir(DOT_FLN)
    os.chmod(DOT_FLN, 0o777)
    os.mkdir(DOT_FLN+'/data')
    os.chmod(DOT_FLN+'/data', 0o777)
    shutil.copy('spacepy/data/spacepy.rc', DOT_FLN+'/')
    shutil.copy('spacepy/data/tai-utc.dat', DOT_FLN+'/data')

pkg_files = ['irbempy/irbempylib.so', 'irbempy/*.py', 'doc/*.*', 'pybats/*.py', 'pybats/*.so'
    'pybats/*.out', 'pycdf/*.py']
#pkg_files = ['onerapy/onerapylib.so','onerapy/*.py', 'doc/*.*', 'pybats/*.py', 'pybats/*.so', 'pycdf/*.py']

# run setup from distutil
setup(name='spacepy',
      version='0.1',
      description='SpacePy: Tools for Space Science Applications',
      author='Steve Morley, Josef Koller, Dan Welling, Brian Larsen, Mike Henderson',
      author_email='spacepy@lanl.gov',
      url='http://www.spacepy.lanl.gov',
      requires=['numpy','scipy','matplotlib (>=0.99)'],
      packages=['spacepy','spacepy.sandbox'],
      package_data={'spacepy': pkg_files},
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Software Development :: Libraries :: Python Modules'
          ]
      )

# update/download packages
if sys.version_info[0]<3:
    if fresh_install:
        dir = update()
        print("Data installed to " + dir)
    else:
        ans = query_yes_no("\nDo you want to update OMNI database and leap seconds table? (Internet connection required)", default = "no")
        if ans=='yes':
            dir = update()
        else:
            print("\nRemember to update OMNI and leap seconds table occasionally by running spacepy.toolbox.update()")
else:
    print('''Updating OMNI and leap seconds on install is not currently supported for Python 3.X.''')

print("\nThanks for installing SpacePy.")
