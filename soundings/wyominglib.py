"""
  module to read the wyoming upperair site and parse it into
  dataframes  -- see newsoundings.ipynb for an example
"""
import requests
from bs4 import BeautifulSoup
import re
import datetime
from  datetime import timezone as tz
import pandas as pd
from collections import OrderedDict
from a405thermo.thermlib import find_esat
from a405thermo.constants import constants as con
import numpy as np
import sys
import h5py


# We need to parse a set of lines that look like this:

#                              Station number: 82965
#                            Observation time: 110201/0000
#                            Station latitude: -9.86
#                           Station longitude: -56.10
#                           Station elevation: 288.0
#                             Showalter index: -0.37
#                                Lifted index: -0.68

# Here's a regular expresion that does that:


re_text="""
           .*Station\snumber\:\s(.+?)\n
           \s+Observation\stime\:\s(.+?)\n
           \s+Station\slatitude\:\s(.+?)\n
           \s+Station\slongitude\:\s(.+?)\n
           \s+Station\selevation\:\s(.+?)\n
           .*
        """

#
# DOTALL says: make . include \n
# VERBOSE says: ignore whitespace and comments within the regular expression
#         to help make the expression more readable
#
header_re=re.compile(re_text,re.DOTALL|re.VERBOSE)

def parse_header(header_text):
    """
    Read a header returned by make_frames and
    return station information

    Parameters
    ----------

    header_text : str
                  string containing wyoming header read with make_frames

    Returns
    -------

    theDate : datetime
              date of sounding
                                  
    the_id  : int
              5 digit souding id

    lat     : float
              latitude (degrees North)

    lon     : float
              longitude (degrees east)
    
    elev    : float
              station elevation (meters)
    """
    header_text=header_text.strip()
    the_match = header_re.match(header_text)
    try:
        the_id,string_time,lat,lon,elev=the_match.groups()
    except AttributeError:
        print('parse failure with: \n',header_text)
    the_id,string_time,lat,lon,elev=the_match.groups()
    elev=elev.split('\n')[0]  #some soundings follow elev with Shoalwater, not Lifted
    lat=float(lat)
    lon=float(lon)
    elev=float(elev)
    day,hour = string_time.strip().split('/')
    print('here is the day: ',day)
    year=int(day[:2]) + 2000
    month=int(day[2:4])
    day=int(day[4:6])
    minute=int(hour[2:])
    hour=int(hour[:2])
    theDate=datetime.datetime(year,month,day,hour,minute,0,0,tz.utc)
    return theDate,the_id,lat,lon,elev

def parse_data(data_text):
    """
    Read a single sounding into a dataframe

    Parameters
    ----------

    data_text : str
                sounding text

    Returns
    -------

    df_out : dataframe
             11 column data frame with sounding values

    unit_name : list
                list of strings with name of units of each column
    
    """

    """
      read lines with 11 numbers and convert to dataframe

      data_text looks like:

        -----------------------------------------------------------------------------
           PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
            hPa     m      C      C      %    g/kg    deg   knot     K      K      K 
        -----------------------------------------------------------------------------
         1000.0    100                                                               
          979.0    288   24.0   23.0     94  18.45      0      0  299.0  353.0  302.2
          974.0    333   25.2   21.1     78  16.46    348      0  300.6  349.1  303.6
          932.0    719   24.0   16.0     61  12.42    243      3  303.2  340.3  305.4
          925.0    785   23.4   15.4     61  12.03    225      3  303.2  339.2  305.4
    """
    all_lines=data_text.strip().split('\n')
    count=0
    theLine=all_lines[count]
    try:
        while theLine.find('PRES   HGHT   TEMP   DWPT') < 0:
            count += 1
            theLine = all_lines[count]
        header_names=all_lines[count].lower().split()
    except IndexError:
        print("no column header line found in sounding")
        sys.exit(1)
    count += 1  #go to unit names
    unit_names=all_lines[count].split()
    count+=2  #skip a row of ------
    data_list=[]
    while True:
        try:
            the_line=all_lines[count]
            dataFields = the_line.split()
            if len(dataFields) == 11:
                try:
                    dataFields = [float(number) for number in dataFields]
                    es = find_esat(dataFields[3] + 273.15)*0.01  #get vapor pressure from dewpoint in hPa
                    dataFields[5] = (con.eps*es/(dataFields[0] - es))*1.e3   #g/kg
                except ValueError:
                    print('trouble converting dataFields to float')
                    print(dataFields)
                    sys.exit(1)
                data_list.append(dataFields)
                #
                # get the next line
                #
            count += 1
            theLine = all_lines[count]
        except IndexError:
            break
    df_out=pd.DataFrame.from_records(data_list,columns=header_names)
    return df_out,unit_names

def make_frames(html_doc):
    """
    turn an html page retrieved from the wyoming site into a dataframe

    Parameters
    ----------

    html_doc : string
               web page from wyoming upperair sounding site
               http://weather.uwyo.edu/cgi-bin/sounding retrieved by
               the requests module

    Returns
    -------

    attr_dict : dict
               attr_dict dictionary with ['header', 'site_id','longitude','latitude', 'elevation', 'units']
              
    sound_dict : dict  
                 sounding dictionary with sounding times as keys and sounding as dataframes
    """  
    soup=BeautifulSoup(html_doc,'html.parser')
    keep=list()
    header= (soup.find_all('h2'))[0].text
    print('header is: ',header)
    for item in soup.find_all('pre'):
        keep.append(item.text)

    sounding_dict=OrderedDict()
    evens=np.arange(0,len(keep),2)
    for count in evens:
        df,units=parse_data(keep[count])
        theDate,the_id,lat,lon,elev=parse_header(keep[count+1])    
        sounding_dict[theDate]=df
        
    attr_dict=dict(units=';'.join(units),site_id=the_id,
                   latitude=lat,longitude=lon,elevation=elev,
                   header = header)
    
    return attr_dict,sounding_dict
    
def test_wyoming():
    """
    function to test downloading a sounding
    from http://weather.uwyo.edu/cgi-bin/sounding and
    creating an hdf file with soundings and attributes

    see the notebook newsoundings.ipynb for use
    """

    url_template=("http://weather.uwyo.edu/cgi-bin/sounding?"
              "region={region:s}"
              "&TYPE=TEXT%3ALIST"
              "&YEAR={year:s}"
              "&MONTH={month:s}"
              "&FROM={start:s}"
              "&TO={stop:s}"
              "&STNM={station:s}")

    values=dict(region='samer',year='2013',month='2',start='0100',stop='2800',station='82965')
    #values=dict(region='nz',year='2013',month='2',start='0100',stop='2800',station='93417')
    #values=dict(region='naconf',year='2013',month='2',start='0100',stop='2800',station='71802')
    #values=dict(region='ant',year='2013',month='07',start='0100',stop='2800',station='89009')
    values=dict(region='naconf',year='2015',month='7',start='1600',stop='1800',station='72451')
    
#naconf, samer, pac, nz, ant, np, europe,africa, seasia, mideast
    url=url_template.format_map(values)

    do_web = True
    backup_file='backup.txt'
    if do_web:
        html_doc = requests.get(url).text
        print(len(html_doc))
        with open(backup_file,'w') as f:
            f.write(html_doc)
        if len(html_doc) < 2000:
            print('debug: short html_doc, something went wrong:',html_doc)
            sys.exit(1)
    else:
        with open(backup_file,'r') as f:
            html_doc=f.read()

    attr_dict,sounding_dict = make_frames(html_doc)
    attr_dict['history']="written by test_requests.py"
    key_list=['header', 'site_id','longitude','latitude', 'elevation', 'units','history']

    name = 'out.h5'    
    with pd.HDFStore(name,'w') as store:
        for key,value in sounding_dict.items():
            thetime=key.strftime("Y%Y_%b_%d_%HZ")
            store.put(thetime,value,format='table')

    with h5py.File(name,'a') as f:
        for key in key_list:
            print('writing key, value: ',key,attr_dict[key])
            f.attrs[key]=attr_dict[key]
        f.close()

    print('hdf file {} written'.format(name))
    print('reading attributes: ')
    with h5py.File(name,'r') as f:
        keys=f.attrs.keys()
        for key in keys:
            try:
                print(key,f.attrs[key])
            except OSError:
                pass


if __name__ == "__main__":
    test_wyoming()
