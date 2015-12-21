# PyWAD is open-source software and consists of a set of plugins written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
# Description:
# This plugin parses the daily QC report (XML format) generated by the Siemens biograph TOF PET-CT.
# To send the data from the scanner to dcm4chee a separate tool that has to be installed on the scanner has been developed by Rob van Rooij and Dennis Dickerscheid.

__version__ = '20151221'
__author__ = 'TdW'



from dicom import tag
import xml.etree.ElementTree as ET
import lxml.etree as etree

def print_xml(xmlroot):
   for child in xmlroot:
        print '=='*20
        print child.tag, child.attrib, child.text
        
        for subchild in child:
            print '\t', subchild.tag, subchild.attrib, subchild.text

            for value in subchild:
                print '\t\t', value.tag, value.attrib, value.text
                for subvalue in value:
                    print '\t\t\t', subvalue.tag, subvalue.attrib, subvalue.text

                    for subsubvalue in subvalue:
                        print '\t\t\t\t', subsubvalue.tag, subsubvalue.attrib, subsubvalue.text

def parseqcreport(data,results,**kwargs):
    params = kwargs.get('params', None)
    p = {}
    for param in params:
        p[param.tag] = (param.text,param.attrib)

    print p
    
    status_ok=True

    relevantfile = data.getAllInstances()[0]
    xmltext = relevantfile[tag.Tag(p.get('use_private_tag')[0].split(','))]

    root = etree.fromstring(xmltext.value)
    print_xml(root)
    
    datagroup = root.find('DataGroup')
    scandate = datagroup.get('DateTime')
    
    datapackages = root.findall(".//DataPackage[@Package_Name='ID_MODE']")
    for data in datapackages:
        status_ok=True
        kV=data.find(".//*/DataEntry[@Variable_Name='ID_VOLT']/LongVariable").get('Long_Value')
        mA=data.find(".//*/DataEntry[@Variable_Name='ID_CURRENT']/LongVariable").get('Long_Value')
        region=data.find(".//*/DataEntry[@Variable_Name='ID_REGION']/StringVariable").get('String_Value').replace('ID_','')

        center_row_n1_ok=data.find(".//*/DataSeries[@Series_Name='ID_CENTER']/SeriesEntry[@Series_ID='ID_ROW_N1']/DataEntry").get('Variable_In_Spec')
        center_row_n2_ok=data.find(".//*/DataSeries[@Series_Name='ID_CENTER']/SeriesEntry[@Series_ID='ID_ROW_N2']/DataEntry").get('Variable_In_Spec')
        center_row_n3_ok=data.find(".//*/DataSeries[@Series_Name='ID_CENTER']/SeriesEntry[@Series_ID='ID_ROW_N3']/DataEntry").get('Variable_In_Spec')
        status_ok = status_ok and center_row_n1_ok=='INTOL' and center_row_n2_ok=='INTOL' and center_row_n3_ok=='INTOL'

        center_row_n1=data.find(".//*/DataSeries[@Series_Name='ID_CENTER']/SeriesEntry[@Series_ID='ID_ROW_N1']/*/DoubleVariable").get('Double_Value')
        center_row_n2=data.find(".//*/DataSeries[@Series_Name='ID_CENTER']/SeriesEntry[@Series_ID='ID_ROW_N2']/*/DoubleVariable").get('Double_Value')
        center_row_n3=data.find(".//*/DataSeries[@Series_Name='ID_CENTER']/SeriesEntry[@Series_ID='ID_ROW_N3']/*/DoubleVariable").get('Double_Value')

        perif_3h_row_n1_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF3']/SeriesEntry[@Series_ID='ID_ROW_N1']/DataEntry").get('Variable_In_Spec')
        perif_3h_row_n2_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF3']/SeriesEntry[@Series_ID='ID_ROW_N2']/DataEntry").get('Variable_In_Spec')
        perif_3h_row_n3_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF3']/SeriesEntry[@Series_ID='ID_ROW_N3']/DataEntry").get('Variable_In_Spec')
        status_ok = status_ok and perif_3h_row_n1_ok=='INTOL' and perif_3h_row_n2_ok=='INTOL' and perif_3h_row_n3_ok=='INTOL'

        perif_3h_row_n1=data.find(".//*/DataSeries[@Series_Name='ID_DIFF3']/SeriesEntry[@Series_ID='ID_ROW_N1']/*/DoubleVariable").get('Double_Value')
        perif_3h_row_n2=data.find(".//*/DataSeries[@Series_Name='ID_DIFF3']/SeriesEntry[@Series_ID='ID_ROW_N2']/*/DoubleVariable").get('Double_Value')
        perif_3h_row_n3=data.find(".//*/DataSeries[@Series_Name='ID_DIFF3']/SeriesEntry[@Series_ID='ID_ROW_N3']/*/DoubleVariable").get('Double_Value')

        perif_6h_row_n1_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF6']/SeriesEntry[@Series_ID='ID_ROW_N1']/DataEntry").get('Variable_In_Spec')
        perif_6h_row_n2_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF6']/SeriesEntry[@Series_ID='ID_ROW_N2']/DataEntry").get('Variable_In_Spec')
        perif_6h_row_n3_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF6']/SeriesEntry[@Series_ID='ID_ROW_N3']/DataEntry").get('Variable_In_Spec')
        status_ok = status_ok and perif_6h_row_n1_ok=='INTOL' and perif_6h_row_n2_ok=='INTOL' and perif_6h_row_n3_ok=='INTOL'

        perif_6h_row_n1=data.find(".//*/DataSeries[@Series_Name='ID_DIFF6']/SeriesEntry[@Series_ID='ID_ROW_N1']/*/DoubleVariable").get('Double_Value')
        perif_6h_row_n2=data.find(".//*/DataSeries[@Series_Name='ID_DIFF6']/SeriesEntry[@Series_ID='ID_ROW_N2']/*/DoubleVariable").get('Double_Value')
        perif_6h_row_n3=data.find(".//*/DataSeries[@Series_Name='ID_DIFF6']/SeriesEntry[@Series_ID='ID_ROW_N3']/*/DoubleVariable").get('Double_Value')

        perif_9h_row_n1_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF9']/SeriesEntry[@Series_ID='ID_ROW_N1']/DataEntry").get('Variable_In_Spec')
        perif_9h_row_n2_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF9']/SeriesEntry[@Series_ID='ID_ROW_N2']/DataEntry").get('Variable_In_Spec')
        perif_9h_row_n3_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF9']/SeriesEntry[@Series_ID='ID_ROW_N3']/DataEntry").get('Variable_In_Spec')
        status_ok = status_ok and perif_9h_row_n1_ok=='INTOL' and perif_9h_row_n2_ok=='INTOL' and perif_9h_row_n3_ok=='INTOL'

        perif_9h_row_n1=data.find(".//*/DataSeries[@Series_Name='ID_DIFF9']/SeriesEntry[@Series_ID='ID_ROW_N1']/*/DoubleVariable").get('Double_Value')
        perif_9h_row_n2=data.find(".//*/DataSeries[@Series_Name='ID_DIFF9']/SeriesEntry[@Series_ID='ID_ROW_N2']/*/DoubleVariable").get('Double_Value')
        perif_9h_row_n3=data.find(".//*/DataSeries[@Series_Name='ID_DIFF9']/SeriesEntry[@Series_ID='ID_ROW_N3']/*/DoubleVariable").get('Double_Value')

        perif_12h_row_n1_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF12']/SeriesEntry[@Series_ID='ID_ROW_N1']/DataEntry").get('Variable_In_Spec')
        perif_12h_row_n2_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF12']/SeriesEntry[@Series_ID='ID_ROW_N2']/DataEntry").get('Variable_In_Spec')
        perif_12h_row_n3_ok=data.find(".//*/DataSeries[@Series_Name='ID_DIFF12']/SeriesEntry[@Series_ID='ID_ROW_N3']/DataEntry").get('Variable_In_Spec')
        status_ok = status_ok and perif_12h_row_n1_ok=='INTOL' and perif_12h_row_n2_ok=='INTOL' and perif_12h_row_n3_ok=='INTOL'

        perif_12h_row_n1=data.find(".//*/DataSeries[@Series_Name='ID_DIFF12']/SeriesEntry[@Series_ID='ID_ROW_N1']/*/DoubleVariable").get('Double_Value')
        perif_12h_row_n2=data.find(".//*/DataSeries[@Series_Name='ID_DIFF12']/SeriesEntry[@Series_ID='ID_ROW_N2']/*/DoubleVariable").get('Double_Value')
        perif_12h_row_n3=data.find(".//*/DataSeries[@Series_Name='ID_DIFF12']/SeriesEntry[@Series_ID='ID_ROW_N3']/*/DoubleVariable").get('Double_Value')


        results.addChar('{0} {1}kV {2}mA'.format(region,kV,mA),'Passed' if status_ok else 'Failed',level=1)

        results.addFloat("{0} {1}kV {2}mA ROW_N1 center".format(region,kV,mA),center_row_n1,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N1 perif_3h".format(region,kV,mA),perif_3h_row_n1,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N1 perif_6h".format(region,kV,mA),perif_6h_row_n1,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N1 perif_9h".format(region,kV,mA),perif_9h_row_n1,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N1 perif_12h".format(region,kV,mA),perif_12h_row_n1,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N2 center".format(region,kV,mA),center_row_n2,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N2 perif_3h".format(region,kV,mA),perif_3h_row_n2,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N2 perif_6h".format(region,kV,mA),perif_6h_row_n2,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N2 perif_9h".format(region,kV,mA),perif_9h_row_n2,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N2 perif_12h".format(region,kV,mA),perif_12h_row_n2,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N3 center".format(region,kV,mA),center_row_n3,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N3 perif_3h".format(region,kV,mA),perif_3h_row_n3,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N3 perif_6h".format(region,kV,mA),perif_6h_row_n3,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N3 perif_9h".format(region,kV,mA),perif_9h_row_n3,level=2,units='HU')
        results.addFloat("{0} {1}kV {2}mA ROW_N3 perif_12h".format(region,kV,mA),perif_12h_row_n3,level=2,units='HU')
