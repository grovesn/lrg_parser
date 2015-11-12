import xml.etree.ElementTree as tree
import sys
import glob
import re

def test_version(root):
  assert root.attrib['schema_version'] == '1.9', 'Schema version must be 1.9'
def check_sequence(sequence):
  assert sequence is not None, 'sequence tag empty'
  pattern = re.match('^[ATCG]+$', sequence)
  assert pattern is not None, 'sequence not valid' 


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
 
  print transcripts[0].attrib['name']
      
      
