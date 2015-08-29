pyWAD -- plugin wrapper for the WAD framework
=============================================

Documentation can be found in the Documentation folder:
e.g.: firefox Documentation/_build/pywad2.html

To create a plugin version of the package that can be used with the WAD Software framework:

python Tools/pywad2zip.py zip

and upload the resulting pywadplugin.py.zip file

To test if the package 'works', the -hello world- example is:
python pywadplugin.py Testing/PackageTesting/input.xml

