#LIBRARIES IMPORTED

import xml.etree.ElementTree as tree
import sys
import glob
import re

#ASSERTION FUNCTIONS

def test_version(root):
  assert root.attrib['schema_version'] == '1.9', 'Schema version must be 1.9'

def check_sequence(sequence):
  assert sequence is not None, 'sequence tag empty'
  pattern = re.match('^[ATCG]+$', sequence)
  assert pattern is not None, 'sequence not valid'

def test_for_none(list, name):
  message = 'No'+ name + 'found'
  assert len(list) > 0, message

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

def get_exon_coords(exons, id, trans):
  intron_dic = {}
  exon_dic = {}
  trans_dic = {}
  dics = [exon_dic, intron_dic, trans_dic]
  for exon in exons:  
    label = exon.attrib['label']
    name = 'exon' + label
    print name
    number = int(label)
    coord = exon.getchildren()
    for points in coord:
      if points.attrib['coord_system'] == id:
        key = name + 'start'
        start = int(points.attrib['start'])
        exon_dic[key] = start
        intron_key = 'intron' + label + 'end'
        intron_dic[intron_key] = start - 1
        
        key = name + 'end'
        end = int(points.attrib['end'])
        exon_dic[key] = end
        intron_number = int(label) + 1
        intron_key = 'intron' + str(intron_number) + 'start'
        intron_dic[intron_key] = end + 1
     
      elif points.attrib['coord_system'] == id + trans:
        trans_dic['cdna'+ label + 'start'] = points.attrib['start']
        trans_dic['cdna'+ label + 'end'] = points.attrib['end']
        
  return dics
      

filenames = glob.glob('/home/swc/Desktop/lrg_parser/files_to_be_analysed/LRG*.xml')

for file in filenames:
  print file
  lrg_tree = tree.parse(file)
  root = lrg_tree.getroot()
  test_version(root)
  fixed_annotation = root[0]
  id = fixed_annotation[0].text
  print 'id = ', id
  fa_children = fixed_annotation.getchildren()

  transcripts = []

  for tags in fa_children:
    if tags.tag == 'sequence':
      sequence = tags.text
      check_sequence(sequence)
      print sequence[0:10]
    elif tags.tag== 'transcript':
      transcripts.append(tags)
      print type(tags)
 
  test_for_none(transcripts, 'transcripts')
      
  for transcript in transcripts:
   trans = transcript.attrib['name']
   print type(transcript)
   exons = transcript.findall('exon')
   test_for_none(exons, 'exons')
   print exons
   
   dics = get_exon_coords(exons, id, trans)
   print dics[0]
   print dics[1]
        
    
