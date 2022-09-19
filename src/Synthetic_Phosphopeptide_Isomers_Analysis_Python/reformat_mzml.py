import xml.etree.ElementTree as ET
import xml.dom as dom
import xml
import click
from lxml import etree

@click.command()
@click.option('--mzml-file')
@click.option('--output-mzml-file')
def main(mzml_file, output_mzml_file):
    """A script to reformat the Bruker converted mzML file to an mzML file that is compatible with OpenSwathWorkflow. This should be done after the ms1 spectra are removed from the converted file."""
    mzml = ET.parse(mzml_file)
    root = mzml.getroot()
    ns = {'hupo' : 'http://psi.hupo.org/ms/mzml'}
    
    #spectrum_list = root.find('hupo:run/hupo:spectrumList',ns)

    spectrum_list = root.find('hupo:mzML/hupo:run/hupo:spectrumList',ns)   # the pyopenms filtered mzml file has a slightly different structure 

    for elem in spectrum_list.findall('hupo:spectrum',ns):
        precursor = elem.find('hupo:precursorList/hupo:precursor',ns)
        selectedIonList = precursor.find('hupo:selectedIonList', ns)
        cvparam = elem.find('hupo:precursorList/hupo:precursor/hupo:selectedIonList/hupo:selectedIon/hupo:cvParam',ns)
        isolation_window = ET.SubElement(precursor, '{http://psi.hupo.org/ms/mzml}isolationWindow')
        mz_value = cvparam.get('value')
        if float(mz_value) == 412.5:
                cv_param_a = ET.SubElement(isolation_window, '{http://psi.hupo.org/ms/mzml}cvParam', attrib = {"cvRef":"MS", "accession":"MS:1000827", "name":"isolation window target m/z", "value": mz_value, "unitAccession":"MS:1000040", "unitName":"m/z", "unitCvRef":"MS"})
                cv_param_b = ET.SubElement(isolation_window, '{http://psi.hupo.org/ms/mzml}cvParam', attrib = {"cvRef":"MS", "accession":"MS:1000828", "name":"isolation window lower offset", "value":"12.5", "unitAccession":"MS:1000040", "unitName":"m/z", "unitCvRef":"MS"})
                cv_param_c =  ET.SubElement(isolation_window, '{http://psi.hupo.org/ms/mzml}cvParam', attrib = {"cvRef":"MS", "accession":"MS:1000829", "name":"isolation window upper offset", "value":"12.5", "unitAccession":"MS:1000040", "unitName":"m/z", "unitCvRef":"MS"})
        else:
                cv_param_a = ET.SubElement(isolation_window, '{http://psi.hupo.org/ms/mzml}cvParam', attrib = {"cvRef":"MS", "accession":"MS:1000827", "name":"isolation window target m/z", "value": mz_value, "unitAccession":"MS:1000040", "unitName":"m/z", "unitCvRef":"MS"})
                cv_param_b = ET.SubElement(isolation_window, '{http://psi.hupo.org/ms/mzml}cvParam', attrib = {"cvRef":"MS", "accession":"MS:1000828", "name":"isolation window lower offset", "value":"13.0", "unitAccession":"MS:1000040", "unitName":"m/z", "unitCvRef":"MS"})
                cv_param_c =  ET.SubElement(isolation_window, '{http://psi.hupo.org/ms/mzml}cvParam', attrib = {"cvRef":"MS", "accession":"MS:1000829", "name":"isolation window upper offset", "value":"13.0", "unitAccession":"MS:1000040", "unitName":"m/z", "unitCvRef":"MS"})
        
        
        precursor.append(isolation_window) 
        isolation_window.append(cv_param_a)
        isolation_window.append(cv_param_b)
        isolation_window.append(cv_param_c)
        precursor.remove(selectedIonList)

    def strip_namespace(xml_file): 
    root = xml_file.getroot()
    for elem in root.getiterator():
        elem.tag = etree.QName(elem.tag).localname

    strip_namespace(mzml)

    mzml.write(output_mzml_file)

if __name__=="__main__":
    main()