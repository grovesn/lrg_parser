LRG Parser

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
