import re

testfull="""
Station number: 89009
                           Observation time: 130701/0000
                           Station latitude: -90.00
                          Station longitude: 0.00
                          Station elevation: 2835.0
                               Lifted index: 22.62
    LIFT computed using virtual temperature: 22.63
      Convective Available Potential Energy: 0.00
             CAPE using virtual temperature: 0.00
                      Convective Inhibition: 0.00
             CINS using virtual temperature: 0.00
                     Bulk Richardson Number: 0.00
          Bulk Richardson Number using CAPV: 0.00
  Temp [K] of the Lifted Condensation Level: 222.03
Pres [hPa] of the Lifted Condensation Level: 606.36
     Mean mixed layer potential temperature: 256.21
              Mean mixed layer mixing ratio: 0.08
              1000 hPa to 500 hPa thickness: 5107.00
Precipitable water [mm] for entire sounding: 0.33
"""
testfull=testfull.strip()[:400]
line = 'Station Baker 1 2009-11-17       1223.0'
re_text='Station\snumber\:\s(.+)'
the_re=re.compile(re_text,re.DOTALL|re.VERBOSE)
print(the_re.findall(testfull))




      



