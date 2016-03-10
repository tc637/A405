To reinstall conda from scratch with the A405 packages:

Go to http://conda.pydata.org/miniconda.html and download the Python 3.5 version

Install into a new, non-existing directory (conda will create it), for example
mini35_test

Make sure this version of miniconda is at the front of your path.  In a windows
cmd shell, type

where conda

and make sure that you are getting mini35_test/bin/conda

For osx/linux bash shells, the command is

which conda

Install all our packages with this command:

make sure you are in the folder containing this readme, then do 

conda install --file class_specs.txt

to get the packages

