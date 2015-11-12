import xml.etree.ElementTree as tree
import sys
import glob

def test_version(root):
  assert root.attrib['schema_version'] == '1.9', 'Schema version must be 1.9'

filenames = glob.glob('/home/swc/Desktop/lrg_parser/files_to_be_analysed/LRG*.xml')

for file in filenames:
  lrg_tree = tree.parse(file)
  root = lrg_tree.getroot()
  test_version(root)
  fixed_annotation = root[0]
  id = fixed_annotation[0].text
  print 'id = ', id
  for child in fixed_annotation.iter('sequence'):
    sequence = child.text
    print 'sequence found'
    
