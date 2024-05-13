#!/usr/bin/python3

import subprocess
import sys
import csv
import os
from datetime import datetime
import pickle
import itertools
import shutil
import gzip
import glob

log_dir = 'log/'
metadata_dir = '00.Metadata/'

def run_single_job_and_wait(job_id, job_general_name, cmd_line, should_capture_stdout = True):
    print_log_style('run_single_job_and_wait - About to launch procedure "' + job_general_name  + '" with job id "' + job_id + '"')
    print_log_style('run_single_job_and_wait - Executing the command: "' + cmd_line  + '"')
    if should_capture_stdout:
      out_fh = open('log/'+ job_general_name + '_' + job_id + '.out', 'w')
    else:
      out_fh = None
    err_fh = open('log/'+ job_general_name + '_' + job_id + '.err', 'w')

    exit_code = subprocess.call(cmd_line, bufsize=- 1, executable = '/bin/bash', stdin=None, stdout=out_fh, stderr = err_fh, preexec_fn=None, close_fds=True, shell=True, cwd=None, env=None, universal_newlines=None, startupinfo=None, creationflags=0, restore_signals=True, start_new_session=False)

    if should_capture_stdout:
      out_fh.close()
    err_fh.close()

    if exit_code != 0:
        print_log_style('run_single_job_and_wait - Error procedure "' + job_general_name + '" with job-id "' + job_id + '" returned a non-zero return code: ' + str(exit_code))
        sys.exit(1)
    else:
        print_log_style('run_single_job_and_wait - Procedure "' + job_general_name + '" with job-id "' + job_id + '" returned 0 (normal completetion)')

def print_log_style(s):
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(str(dt_string) + '|' + s)


def produceFastaTableFrom_ampAbundTable(fn, min_count, ref_amp_seq_dict):
  print_log_style('Extracting sequences from file "' + fn + '" an producing fasta file..')
  
  ifh = open(fn, 'r')
  csvreader = csv.reader(ifh, delimiter = "\t")
  next(csvreader)
  
  ofh_dict = {}
  i = 0
  
  
  
  for line in csvreader:
    count = int(line[0])
    seq = line[1]
    seqname = line[2]
    
    if seqname == 'Unidentified':
      continue
      
    if count < min_count:
      continue
    i += 1
    curr_ofn = 'amp_seq_' + seqname + '.fa'
    
    if not curr_ofn in ofh_dict:
      ofh_dict[curr_ofn] = open(curr_ofn, 'w')
      my_ref_seq = ref_amp_seq_dict[seqname]
      ofh_dict[curr_ofn].write('>ref_seq\t0\n' + my_ref_seq + '\n')
    
    ofh_dict[curr_ofn].write('>seq_index_' + str(i) + '\t' + str(count) + '\n' + seq + '\n')
    
  for fh in ofh_dict.values():
    fh.close()
  ifh.close()

def peform_msa():
  my_fns = glob.glob('amp_seq_*.fa')
  
  for ifn in my_fns:
    print_log_style('Performing MSA for sequences in file "' + ifn + '"..')
    ofn = ifn[:-2] + 'align'
    cmd = 'mafft ' + ifn + ' > ' + ofn
    
    run_single_job_and_wait(job_id = ifn, job_general_name = 'align_whole_amp', cmd_line = cmd, should_capture_stdout = True)    
 
def read_ref_ampseq():
  my_dict = {}
  fh = open('../00.Metadata/amplicon_sequences.tsv', 'r')
  
  csvreader = csv.reader(fh, delimiter = "\t")
  
  next(csvreader)    
  
  for line in csvreader:   
    amplicon_id = line[0]
    seq = line[1]
    my_dict[amplicon_id] = seq
  
  
  fh.close()
  
  return my_dict

def getCoordinatesForRefBases(aligned_ref_seq):
  print(aligned_ref_seq)
  i = -1
  aligned_red_seq_existing_pos = []
  for b in aligned_ref_seq:
    i += 1
    if b != '-':
      aligned_red_seq_existing_pos.append(i)
  return aligned_red_seq_existing_pos

def addInfoToBaseCountDict(myBaseCount_dict, curr_seq, my_ref_seq_base_coordinates, seq_count):
  #print(curr_seq)
  for coord in myBaseCount_dict.keys():
    b = curr_seq[coord].upper()
    if b == '-':
      b = 'Del'
    myBaseCount_dict[coord][b] += seq_count
    
  return myBaseCount_dict

def writeBaseCountDictToMatrixFile(myBaseCount_dict, my_ref_seq_base_coordinates, out_fn):
  out_fh = open(out_fn, 'w')
  csvwriter = csv.writer(out_fh, delimiter = "\t")
  csvwriter.writerow(['Pos', 'Ref', 'Major', 'A', 'C', 'T', 'G', 'N', 'Del'])
  Pos = 0
  for coord in my_ref_seq_base_coordinates:
    Pos += 1
    csvwriter.writerow([str(Pos), myBaseCount_dict[coord]['Ref'], myBaseCount_dict[coord]['Major'], myBaseCount_dict[coord]['A'], myBaseCount_dict[coord]['C'], myBaseCount_dict[coord]['T'], myBaseCount_dict[coord]['G'], myBaseCount_dict[coord]['N'], myBaseCount_dict[coord]['Del']])
    
  out_fh.close()
    
def identifyMajorAllelesInCountDict(myBaseCount_dict):
  for coord in myBaseCount_dict.keys():
    ACTGND = [myBaseCount_dict[coord]['A'], myBaseCount_dict[coord]['C'], myBaseCount_dict[coord]['T'], myBaseCount_dict[coord]['G'], myBaseCount_dict[coord]['N'], myBaseCount_dict[coord]['Del']]
    max_value = max(ACTGND)
    peaks_count = ACTGND.count(max_value)
    if peaks_count == 1:
      myBaseCount_dict[coord]['Major'] = 'ACTGND'[ACTGND.index(max_value)]
  
  return myBaseCount_dict
  
def convertCountDict2Percents(myBaseCount_dict):
  for coord in myBaseCount_dict.keys():
    ACTGND_Sum = myBaseCount_dict[coord]['A'] + myBaseCount_dict[coord]['C'] + myBaseCount_dict[coord]['T'] + myBaseCount_dict[coord]['G'] + myBaseCount_dict[coord]['N'] + myBaseCount_dict[coord]['Del']
    
    for b in ['A', 'C', 'T', 'G' ,'N' ,'Del']:
      myBaseCount_dict[coord][b] = float(myBaseCount_dict[coord][b]) / ACTGND_Sum * 100
  
  return myBaseCount_dict


def quantifyBasesAlongAmpliconSeq():
  align_fns = glob.glob('amp_seq_*.align')
  
  for align_fn in align_fns:
    print_log_style('Now quantifying bases in amplicon sequences for "' + align_fn + '"')
    ofn = align_fn[ : -6] + '_baseCountAlongSeq.tsv'
    ifh = open(align_fn, 'r')
    curr_seq = ''
    seqname = '*UNDEF*'
    for line in ifh:
      if len(line) > 0 and line[0] == '>':
        
        line = line [1:]
        line_split = line.split("\t")
        seq_count = int(line_split[1])

        if seqname == 'ref_seq':
          print('"' + seqname + '"')
          ref_seq = curr_seq
          my_ref_seq_base_coordinates = getCoordinatesForRefBases(ref_seq)
          myBaseCount_dict = {}
          for coord in my_ref_seq_base_coordinates:
            myBaseCount_dict[coord] = {'A' : 0, 'C' : 0, 'T' : 0, 'G' : 0, 'N' : 0, 'Del' : 0, 'Ref' : ref_seq[coord].upper(),'Major' : ''}
             
            
            
        elif seqname != '*UNDEF*':
          myBaseCount_dict = addInfoToBaseCountDict(myBaseCount_dict, curr_seq, my_ref_seq_base_coordinates, seq_count)
          
        curr_seq = ''
        seqname = line_split[0]
        
      else:
        curr_seq = curr_seq + line.rstrip()

    myBaseCount_dict = identifyMajorAllelesInCountDict(myBaseCount_dict)
    writeBaseCountDictToMatrixFile(myBaseCount_dict = myBaseCount_dict, my_ref_seq_base_coordinates = my_ref_seq_base_coordinates, out_fn = ofn)
    myBaseCount_dict = convertCountDict2Percents(myBaseCount_dict)
    ofn = ofn[ : -4] + '_percents.tsv'
    writeBaseCountDictToMatrixFile(myBaseCount_dict = myBaseCount_dict, my_ref_seq_base_coordinates = my_ref_seq_base_coordinates, out_fn = ofn)
    
    
    
    
    ifh.close()
    

min_count = 10
ref_amp_seq_dict = read_ref_ampseq()
print(ref_amp_seq_dict)

#produceFastaTableFrom_ampAbundTable(fn = 'amplicon_abundance_table_w_primer_ids.tsv', min_count = min_count, ref_amp_seq_dict = ref_amp_seq_dict)
#peform_msa() 
quantifyBasesAlongAmpliconSeq()