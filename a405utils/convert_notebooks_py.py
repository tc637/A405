#!/usr/bin/env python
"""
   convert all ipython notebooks to python scripts (first time)
   convert all new or modified notebooks (subsequently)

   usage:

   ./convert_notebooks.py  notebook_folder

   result:

      notebook_folder/python  contains the converted notebooks
"""

import os
import glob, stat
import datetime
import tzlocal  #pip install tzlocal
import pytz
import subprocess
import argparse
import textwrap
import sys, traceback


def find_modtime(the_file):
    """
       remove the .py or .ipynb extenstion from the file name
       to get the head, and return that name, plus the modification
       date in UTC.  
    """
    head,ext=os.path.splitext(the_file)
    print('finding modtime for {}'.format(head))
    #
    #  see os.stat docs for the format of the stat function.  It returns
    #  multiple fields (owner, date created, size, etc.) that are indexed by the stat object
    #
    the_date=datetime.datetime.fromtimestamp(os.stat(the_file)[stat.ST_MTIME])
    #
    # finding the local timezone is suprisingly hard -- need to install a
    # special module called tzlocal using pip install tzlocal
    #
    local_tz = tzlocal.get_localzone()
    the_date=local_tz.localize(the_date)
    #
    # convert every date to UTC
    #
    the_date = the_date.astimezone(pytz.utc)
    #
    # remove everything but the root filename
    #
    head = head.split('/')[-1] 
    return head,the_date

if __name__ == "__main__":
    
    linebreaks=argparse.RawTextHelpFormatter
    descrip=textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,description=descrip)
    parser.add_argument('--folder','-f',type=str,help='folder containing ipynb files (default current directory)',default='.')
    args=parser.parse_args()
    try:
        currdir=os.getcwd()
        os.chdir(args.folder)
        #
        # make the python folder if it doesn't exist
        #
        python_files='python'
        if not os.path.exists(python_files):
            print('creating python folder')
            os.makedirs(python_files)

        #
        # get the ipynb and py files
        #
        notebooklist=glob.glob('*.ipynb')
        pythonlist=glob.glob('./python/*.py')
        #
        #  build two dictionaries containing the root names and the modtimes
        #
        py_dict={}
        nb_dict={}
        for the_file in notebooklist:
            head, the_date = find_modtime(the_file)
            nb_dict[head] = the_date

        for the_file in pythonlist:
            head, the_date = find_modtime(the_file)
            py_dict[head] = the_date

        # use python sets to
        #find notebooks not in pythonlist
        #

        py_files=set(py_dict.keys())
        nb_files=set(nb_dict.keys())
        print('{py_files!r}'.format_map(vars()))
        #
        # find all notebooks that don't have py files
        #
        make_py=nb_files - py_files

        #
        #  run ipython nbconvert on these notebook and put output in python folder
        #
        print('rebuilding {}'.format(make_py))
        cmdstring='ipython nbconvert --stdout --to python {0:s}.ipynb > python/{0:s}.py'
        for the_file in make_py:
            command=cmdstring.format(the_file)
            print('executing: ',command)
            out=subprocess.getstatusoutput(command)
            print(out)
        #
        # now go through notebooks once more, converting all that are newer than their python files
        #
        pythonlist=glob.glob('./python/*.py')
        for the_file in pythonlist:
            head, the_date = find_modtime(the_file)
            py_dict[head] = the_date

        for the_root in nb_files:
            nb='{}.ipynb'.format(the_root)
            py='python/{}.py'.format(the_root)
            if nb_dict[the_root] > py_dict[the_root]:
                print('rebuilding {}'.format(the_root))
                command=cmdstring.format(the_root)
                out=subprocess.getstatusoutput(command)
                print(out)
    #
    # trap all exceptions and print a traceback
    #
    except Exception as e:
         ex_type, ex_val, tb = sys.exc_info()
         print('bummer: ',ex_val)
         print('\nhere is the traceback:\n')
         traceback.print_tb(tb)
    #
    # make sure we get back to the original folder
    #
    finally:
        os.chdir(currdir)
