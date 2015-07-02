import sys
from dicom import tag
import lxml.html
import xml.etree.ElementTree as ET
import lxml.etree as etree

def print_xml(xmlroot):
   for child in xmlroot:
        print '=='*20
        print child.tag, child.attrib, child.text
        
        for subchild in child:
            print '1 \t', subchild.tag, subchild.attrib, subchild.text

            for value in subchild:
                print '2 \t\t', value.tag, value.attrib, value.text
                for subvalue in value:
                    print '3 \t\t\t', subvalue.tag, subvalue.attrib, subvalue.text

                    for subsubvalue in subvalue:
                        print '4 \t\t\t\t', subsubvalue.tag, subsubvalue.attrib, subsubvalue.text

#flood field uniformity
#spatial linearity
#slice profile
#spatial resolution

def section_SpatialResulution(html):
   #name
   #date
   #list type
   #applied verif files
   

   #Table
   #patient
   #scan name
   #scan date 
   # ...

   #'Field Strength','Gradient Chain','Gradient Coil','Magnet','Tested by'

    return


def split_piqt_report(html):
   return ''.join(open(html).readlines()).split('\r\n<br><b>')

def parse_html(html):
   output = []
   tmplist = []


   resultdict = {}

   html2etree = lxml.html.fromstring(html)
   for child in html2etree.getchildren():
      if child.tag in resultdict.keys():
         resultdict[child.tag].append(child)
      else:
         resultdict[child.tag] = [child]

   if 'font' in resultdict.keys():
      name = resultdict['font'][0].text_content()
      #print 'naam:', naam
   
   tmpout = []
   if 'table' in resultdict.keys():
      for table in resultdict['table']:
         tmpdict = {}
         if name in ['Flood Field Uniformity','Spatial Linearity','Slice Profile']:
            print 'name',name

            if name == 'Flood Field Uniformity':
               prename = 'FFU'
            if name == 'Spatial Linearity':
               prename = 'SL'
            if name == 'Slice Profile':
               prename = 'SP'


            rows = table.findall("tr")
            for row in rows:
               output.append(tmpdict)
               tmprow = [prename+'_'+c.text_content() for c in row.getchildren()]
               tmpkey = tmprow[0]
               #print tmpkey

               if tmpkey == 'Scan_Name':
                  tmprow = row.getchildren()
                  scanlayout = [c.attrib['colspan'] for c in tmprow]
                  print 'scanlayout',scanlayout


               if tmpkey == 'Echo_No':
                  tmprow = row.getchildren()
                  echolayout = [c.attrib['colspan'] for c in tmprow]
                  print 'echolayout',echolayout



               if tmpkey in ['Field Strength','Gradient Chain','Gradient Coil','Magnet','Tested by']:
                  tmpdict[tmprow[0]] = tmprow[1]
                  tmpdict[tmprow[2]] = tmprow[3]
               else:
                  tmpdict[tmpkey]=tmprow[1:]
            tmpout.append(tmpdict)

   print tmpout
   #for key in tmpdict.keys():
   #   print key,tmpdict[key]

   headerpars = ['']

   outdict = {}

   #for key in tmpdict.keys():
   #   print key

   return tmpout



def parseqcreport(data,results,**kwargs):
    params = kwargs.get('params', None)
    p = {}
    for param in params:
        p[param.tag] = (param.text,param.attrib)

    print p
    relevantfile = data.getAllInstances()[0]
    xmltext = relevantfile[tag.Tag(p.get('use_private_tag')[0].split(','))]

    root = etree.fromstring(xmltext.value)
    print_xml(root)

    #Sections:
    #Title

    #Scandate


    #Phantomparameters
    #phantompars = root.find('cPhantomParameters')
    #Isotope = phantompars.find('aIsotope').text
    #results.addChar('Isotope',Isotope,level=1)

    #results.addFloat('BlockTimingWidth',BlockTimingWidth,level=2)
    #TimeAlignmentResidual =  detres.find('lTAResidual').find('cBlkValue').find('aValue').text
    #results.addFloat('Time alignment residual',TimeAlignmentResidual,level=2)


def main():
   filename = '/data/SoftwareTools/WAD_Software/WADplugins/pywad2/Tools/PIQT/testdata/SPT report2.htm'
   
   splitfile = split_piqt_report(filename)



   for elem in splitfile:
      parse_html(elem)


if __name__ == "__main__":
   sys.exit(main())
