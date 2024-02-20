#!/usr/bin/python3
import subprocess
import sys
import csv
import os
from datetime import datetime
import pickle
import itertools
import shutil

log_dir = 'log/'
#cut_off
#conda_init_script_path =  str(pipeline_settings_dict['analysis_settings_dict']['conda_init_string_path']['Value'])
conda_init_script_path = '/home/synthesizer/miniconda3/etc/profile.d/conda.sh'
def merge_r1_r2_reads(sample_id):
  cmd = 'source /home/anaconda/miniconda3/bin/activate ngmerge_env && NGmerge -1 ' + sample_id + '_header_corrected_r1.fq -2 ' + sample_id + '_header_corrected_r2.fq -o ' + sample_id + '_merged.fq.gz -l ' + sample_id + '_stich.log -f ' + sample_id + '_failed_stich.fq.gz -j ' + sample_id + '_stitch_alignment.txt'
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'merge_r1_r2', cmd_line = cmd)
  

def is_a_gz_file(fn):
  if fn[-3:] == '.gz' or fn[-5:] == '.gzip' or fn[-3:] == '.GZ' and fn[-5:] == '.GZIP':
    return True
  else:
    return False 

def fix_mgi_fastq_headers(sample_id, fq1fn, fq2fn):
  if is_a_gz_file(fq1fn):
    cat_tool = 'zcat'
  else:
    cat_tool = 'cat'
  cmd = cat_tool + " " + fq1fn + '|awk \'{if (NR%4==1) {gsub(/\//, " ");print } else print}\'> ' + sample_id + "_header_corrected_r1.fq"
  run_single_job_and_wait(job_id = 's1', job_general_name = 'fix_fastq_headers_r1', cmd_line = cmd)

  if is_a_gz_file(fq2fn):
    cat_tool = 'zcat'
  else:
    cat_tool = 'cat'
  cmd = cat_tool + " " + fq2fn + '|awk \'{if (NR%4==1) {gsub(/\//, " ");print } else print}\'> ' + sample_id + "_header_corrected_r2.fq"
  run_single_job_and_wait(job_id = 's1', job_general_name = 'fix_fastq_headers_r2', cmd_line = cmd)

def  read_primers_info():
  primers_dict = {}
  inv_primers_dict = {}
  
  fh = open('r1_primers.tsv', 'r')
  
  csvreader = csv.reader(fh, delimiter = "\t")
  next(csvreader)
  
  for line in csvreader:
    primer_id = line[0]
    primer_seq = line[1]
    primers_dict[primer_id] = primer_seq 
    inv_primers_dict[primer_seq] = primer_id
  
  fh.close()

  return primers_dict, inv_primers_dict
  
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


def produce_initial_amp_abundance_table(sample_id):
  cmd = "cat " + sample_id + "_merged.fq.gz|awk \'{if (NR %  4 == 2) print}'|sort|uniq -c| tr -s ' '|tr ' ' $'\t'|cut -c2-|sort -rnk1 >" + sample_id + "_amplicon_abundance_table.tsv"
  run_single_job_and_wait(job_id = 's1', job_general_name = 'produce_amp_abundance_table', cmd_line = cmd)

def choose_matching_prefix(line, prefixes_dict):
  ret_val = 'Unidentified'
  for k in prefixes_dict.keys():
    v = prefixes_dict[k]
    if line[0:len(k)] == k:
      ret_val = v
      break
  return ret_val

def assign_seqs2primers_ids(inv_primers_dict, sample_id):
  ifh = open(sample_id + '_amplicon_abundance_table.tsv', 'r')
  ofh = open(sample_id + '_amplicon_abundance_table_w_primer_ids.tsv', 'w')
  
  csvreader = csv.reader(ifh, delimiter = "\t")
  csvwriter = csv.writer(ofh, delimiter = "\t")
  
  for line in csvreader:
    #print(line)
    primer_id = choose_matching_prefix(line[1], inv_primers_dict)
    line.append(primer_id)
    #print(line)
    csvwriter.writerow(line)
  
  ifh.close()
  ofh.close()

def produce_primers_dist_table_and_pie_plot(sample_id):
  cmd = "Rscript " + scripts_dir + '/visualize_primer_dist.R ' + sample_id
  run_single_job_and_wait(job_id = 's1', job_general_name = 'produce_primers_dist_table_and_pie_plot', cmd_line = cmd)

scripts_dir = os.path.dirname(os.path.realpath(__file__))
fix_mgi_fastq_headers(sample_id = 's1', fq1fn = '../raw_data/Carmi_Nir_484/S250051509_L01_NGS484_s1_amp_M1_1.fq.gz', fq2fn = '../raw_data/Carmi_Nir_484/S250051509_L01_NGS484_s1_amp_M1_2.fq.gz')
merge_r1_r2_reads(sample_id = 's1')
produce_initial_amp_abundance_table(sample_id = 's1')
primers_dict, inv_primers_dict = read_primers_info()
assign_seqs2primers_ids(inv_primers_dict, sample_id = 's1')
produce_primers_dist_table_and_pie_plot(sample_id = 's1')