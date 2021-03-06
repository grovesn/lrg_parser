'''lrg_parser.py
authors: Kelly E, Natalie G, Kirsty R

This script is designed to extract the intron sequences from an LRG file and putput this as either an XML of Fasta file

To run the script, the locations of the "files_to_be_analysed" and "analysis_results" directories must be set in the script.

By default it is assumed that these directories are within the same directory as the script is run from e.g. "./files_to_be_analysed"

The script also requires an arguement that defines the output required. This must either be "fasta" or "xml" 

Fasta output will generate a .fa file with the ID of the LRG file containing a header and the corresponding sequence for each intron.

The header will be in the following format: ">id=LRG_7|transcript_name=t1|intron_number=12" and the following line will contain the sequence

XML files will contain the information within the header of the .fa file but structured as follows: 

<id name="LRG_7">
  <transcript name="t1">
    <intron length="2034bp">
      <sequence>
'''






#LIBRARIES IMPORTED

import xml.etree.ElementTree as tree
#imports ElementTree for xml parsing
import sys
#imports sys to access command line arguments
import glob
#imports glob to read multiple files
import re
#imports regular expression for pattern matching
import string


#ASSERTION FUNCTIONS

def test_version(root):
  assert root.attrib['schema_version'] == '1.9', 'Schema version must be 1.9'
#assert function to ensure correct file version is used


def check_sequence(sequence):
  assert sequence is not None, 'sequence tag empty'
  pattern = re.match('^[ATCG]+$', sequence)
  assert pattern is not None, 'sequence not valid'
#assert function to ensure only ACGTs are found in the seq


def test_for_none(list, name):
  message = 'No '+ name + ' found'
  assert len(list) > 0, message
#assert function to check that the transcript and exon lists are not empty


def test_dics(dics):
  count = len(dics[0])/2 
  i = range(1,count + 1)
  exon_dic = dics[0]
  intron_dic = dics[1]
  trans_dic = dics[2]
  assert intron_dic['intron1end'] == 5000, 'Length of intron 1 is not 5000 as expected'
  for exon in i:
    exon_start = 'exon' + str(exon) + 'start'
    exon_end = 'exon' + str(exon) + 'end'
    cdna_start = 'cdna' + str(exon) + 'start'
    cdna_end = 'cdna' + str(exon) + 'end'
    intron_start = 'intron' + str(exon + 1) + 'start'
    intron_end = 'intron' + str(exon) + 'end'
    assert exon_dic[exon_end] - exon_dic[exon_start] == trans_dic[cdna_end] - trans_dic[cdna_start], 'Exon coordinates for exon ' + str(exon) + ' do not match with cdna coordinates'
    assert exon_dic[exon_start] == intron_dic[intron_end]+1, 'Exon and intron coordinates do not match for exon ' + exon 
    assert exon_dic[exon_end] == intron_dic[intron_start]-1, 'Exon and intron coordinates do not match for exon ' + exon
#ensures the dictionary values have been assigned to the correct keys

def get_exon_coords(exons, id, trans):
  intron_dic = {}
  exon_dic = {}
  trans_dic = {}
  numbers = []
  dics = [exon_dic, intron_dic, trans_dic, numbers]
#create an array for dictionaries to pass back to main code

  for exon in exons:  
    label = exon.attrib['label']
    number = 0
    if re.search('[a-zA-Z]', label):
      letter = label[1:]
      point = string.lowercase.index(letter)
      label = label[0:len(label)-1] + '.' + str(point)
      number = float(label)
    else:
      number = int(label)
    numbers.append(number)
    coord = exon.getchildren()
    name = 'exon' + label
#sets variables for keys in the dictionary

    for points in coord:
      if points.attrib['coord_system'] == id:
        key = name + 'start'
        start = int(points.attrib['start'])
        exon_dic[key] = start
        intron_key = 'intron' + str(number) + 'end'
        intron_dic[intron_key] = start - 1
#gets coords for exon start and intron end        

        key = name + 'end'
        end = int(points.attrib['end'])
        exon_dic[key] = end
        intron_number = number + 1
        intron_key = 'intron' + str(number) + 'start'
        intron_dic[intron_key] = end + 1
#gets coords for exon end and intron start

      elif points.attrib['coord_system'] == id + trans:
        trans_dic['cdna'+ label + 'start'] = points.attrib['start']
        trans_dic['cdna'+ label + 'end'] = points.attrib['end']
#adds coords for cdna values for testing of exon variables        

  return dics
      

#EXPORT FUNCTIONS

def export_xml(header):
  export_root = tree.Element(id)
    
  trans_attrib = {'name' : trans}
  trans_tag = tree.SubElement(export_root, 'transcript', trans_attrib)
  
  for intron in ultimate_dic:
    intron_attrib = {'number' : str(intron[0]), 'start' : str(intron[1]), 'end' : str(intron[2])}
    intron_tag = tree.SubElement(trans_tag, 'intron', intron_attrib)
    intron_tag.text = intron[3]
  
  xml_file = filepath + 'analysis_results/' + id + '_results.xml'
  file = open(xml_file, 'w')  
  file.write(tree.tostring(export_root))

def export_fasta():
  fasta_file = filepath + 'analysis_results/' + id + '_results.fa'
  file = open(fasta_file, 'w')
  for intron in ultimate_dic:
    intron_header = '>id=' + id + '|transcript=' + trans + '|intron_number=' + str(intron[0]) + '|start=' + str(intron[1]) + '|end=' + str(intron[2]) + '\n'
    intron_sequence = intron[3] + '\n'

    file.write(intron_header)
    file.write(intron_sequence)



#CODE START
assert len(sys.argv) == 2, 'Please specify either xml or fasta as export type'

args = sys.argv
export_type = None
if args[1] == 'xml':
  export_type = 'xml'
elif args[1] == 'fasta':
  export_type = 'fasta'
else:
  print 'Please specify either xml or fasta as export type'
  exit()

filepath = '/home/swc/Desktop/lrg/'
filenames = glob.glob(filepath + 'files_to_be_analysed/LRG*.xml')
#imports all files in the specified folder

for file in filenames:
  print file
  lrg_tree = tree.parse(file)
  root = lrg_tree.getroot()
  test_version(root)
  fixed_annotation = root[0]
  id = fixed_annotation[0].text
  print 'id = ', id
  fa_children = fixed_annotation.getchildren()
#parses xml file and creates variables for important elements

  transcripts = []
#creates array for transcripts

  for tags in fa_children:
    if tags.tag == 'sequence':
      sequence = tags.text
      check_sequence(sequence)
      print sequence[0:10]
      print len(sequence)
    elif tags.tag== 'transcript':
      transcripts.append(tags)
#extracts sequence from xml and adds transcript elements to array

  test_for_none(transcripts, 'transcripts')
#ensures array is not empty
      
  for transcript in transcripts:
   trans = transcript.attrib['name']
   exons = transcript.findall('exon')
   test_for_none(exons, 'exons')
#ensures there are exons in the xml   

   dics = get_exon_coords(exons, id, trans)
   exon_dic = dics[0]
   intron_dic = dics[1]
   print intron_dic
   numbers = dics[3]
   print numbers
#assign dictionaries in  array to the correct variable names

   number_introns = len(intron_dic)/2 +1
#get the number of introns - this is one less than the actual number as the start of intron 1 and end of the last intron are not added to the dictionary
  
#create an array of values to loop through each intron sequentially - must be one more than the number of introns so that the last intron number is included in the for loop
   ultimate_dic = []
   for intron in numbers:
     print intron
     start_key = 'intron' + str(intron) + 'start' #define the keys for the start and end of the intron
     end_key = 'intron' + str(intron) + 'end'
     if intron == 1: #the start of intron one is not included in the dictionary as it is the beginning of the sequence
       start = 0
       end = intron_dic[end_key]
       intron_header = 'intron 1 sequence:'
 #extract the range of sequence that corresponds to the intron
     elif intron == number_introns: #the end of the last intron is not included in the dictionary as it is the end of the sequence
       start = intron_dic[start_key]
       end = len(sequence)
       intron_header = 'intron ' + str(intron) + ' sequence:'
     else:
       start = intron_dic[start_key]
       end = intron_dic[end_key]
       intron_header = 'intron ' + str(intron) + ' sequence:'
     intron_sequence = sequence[start:end + 1]
     value = [intron, start, end, intron_sequence]
     ultimate_dic.append(value)
     print intron_header
     print intron_sequence
   if export_type == 'xml':
     export_xml()
   elif export_type == 'fasta':
     export_fasta()
