#LIBRARIES IMPORTED

import xml.etree.ElementTree as tree
#imports ElementTree for xml parsing
import sys
#imports sys to access command line arguments
import glob
#imports glob to read multiple files
import re
#imports regular expression for pattern matching



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
  message = 'No'+ name + 'found'
  assert len(list) > 0, message
#assert function to check that the transcript and exon lists are not empty


def test_dics(dics):
  count = sum(dics[0].values())/2 
  i = range(1,count)
  exon_dic = dics[0]
  intron_dic = dics[1]
  trans_dic = dics[2]
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
  dics = [exon_dic, intron_dic, trans_dic]
#create an array for dictionaries to pass back to main code
  
  for exon in exons:  
    label = exon.attrib['label']
    name = 'exon' + label
    print name
    number = int(label)
    coord = exon.getchildren()
#sets variables for keys in the dictionary

    for points in coord:
      if points.attrib['coord_system'] == id:
        key = name + 'start'
        start = int(points.attrib['start'])
        exon_dic[key] = start
        intron_key = 'intron' + label + 'end'
        intron_dic[intron_key] = start - 1
#gets coords for exon start and intron end        

        key = name + 'end'
        end = int(points.attrib['end'])
        exon_dic[key] = end
        intron_number = int(label) + 1
        intron_key = 'intron' + str(intron_number) + 'start'
        intron_dic[intron_key] = end + 1
#gets coords for exon end and intron start

      elif points.attrib['coord_system'] == id + trans:
        trans_dic['cdna'+ label + 'start'] = points.attrib['start']
        trans_dic['cdna'+ label + 'end'] = points.attrib['end']
#adds coords for cdna values for testing of exon variables        

  return dics
      

filenames = glob.glob('/home/swc/Desktop/lrg_parser/files_to_be_analysed/LRG*.xml')
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
    elif tags.tag== 'transcript':
      transcripts.append(tags)
      print type(tags)
#extracts sequence from xml and adds transcript elements to array

  test_for_none(transcripts, 'transcripts')
#ensures array is not empty
      
  for transcript in transcripts:
   trans = transcript.attrib['name']
   print type(transcript)
   exons = transcript.findall('exon')
   test_for_none(exons, 'exons')
   print exons
#ensures there are exons in the xml   

   dics = get_exon_coords(exons, id, trans)
   print dics[0]
   print dics[1]
#prints the coords for all introns and exons
    
