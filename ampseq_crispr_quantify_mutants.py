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

log_dir = 'log/'
metadata_dir = '00.Metadata/'
#cut_off
#conda_init_script_path =  str(pipeline_settings_dict['analysis_settings_dict']['conda_init_string_path']['Value'])
conda_init_script_path = '/home/synthesizer/miniconda3/etc/profile.d/conda.sh'

def unify_individual_samples_tables():
  print_log_style('About to unify tables from different samples..')
  cmd = 'Rscript ' + scripts_dir + '/unify_individual_samples_tables.R'
  run_single_job_and_wait(job_id = 'all', job_general_name = 'unify_individual_samples_tables', cmd_line = cmd)

def visualize_mutation_rates():
  print_log_style('Visualizing results..')
  cmd = 'Rscript ' + scripts_dir + '/visualize_mutation_rates.R'
  run_single_job_and_wait(job_id = 'all', job_general_name = 'visualize_mutation_rates', cmd_line = cmd)



def copy_primer_dist_image_files(my_pipeline_settings_dict):
  for sample_id in my_pipeline_settings_dict['sample_seq_files_dict'].keys():
    shutil.copyfile(sample_id + '/' + 'primers_dist.png', 'plots/' + 'primers_dist_' + sample_id + '.png')
    

def countLinesInGZIP(fn):
  with gzip.open(fn, 'rb') as f:
    for i, l in enumerate(f):
      pass
  return i+1

def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

def countLinesInTextFile(fn):
  with open(fn, "r",encoding="utf-8",errors='ignore') as f:
    return sum(bl.count("\n") for bl in blocks(f))

def countLinesInFile(fn):
  if is_a_gz_file(fn):
    return countLinesInGZIP(fn)
  else:
    return countLinesInTextFile(fn)
    
def producePercentStitchedReport(my_pipeline_settings_dict):
  print_log_style('Counting reads in fastq files..')
  my_report_table = [['sample_id', 'reads_count', 'reads_stitched', 'percent_reads_stitched']]

  for sample_id in my_pipeline_settings_dict['sample_seq_files_dict'].keys():
    print_log_style('Processing sample "' + sample_id + '"..')
    fq1fn = pipeline_settings_dict['sample_seq_files_dict'][sample_id]['r1_filename']
    reads_count_original = countLinesInFile(fq1fn) / 4
    reads_count_stitched = countLinesInTextFile(sample_id + '/' + sample_id + '_merged.fq.gz') / 4
    percent_stitched = reads_count_stitched / reads_count_original * 100
    my_report_table.append([sample_id, reads_count_original, reads_count_stitched, percent_stitched])
  
  fh = open('reads_merge_statistics.tsv', 'w')
  csvwriter = csv.writer(fh, delimiter = "\t")
  csvwriter.writerows(my_report_table)
  fh.close()
    

def make_req_dirs():
  os.makedirs(log_dir, exist_ok = True) 
  os.makedirs(metadata_dir, exist_ok = True)
  os.makedirs('plots', exist_ok = True)

def merge_r1_r2_reads(sample_id):
  print_log_style('About to merge r1 and r2 read files for sample "' + sample_id + '"..')
  cmd = 'source /home/anaconda/miniconda3/bin/activate ngmerge_env && NGmerge -1 ' + sample_id + '_header_corrected_r1.fq -2 ' + sample_id + '_header_corrected_r2.fq -o ' + sample_id + '_merged.fq.gz -l ' + sample_id + '_stich.log -f ' + sample_id + '_failed_stich.fq.gz -j ' + sample_id + '_stitch_alignment.txt'
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'merge_r1_r2', cmd_line = cmd)

def convertExcelFile2Pickle(scripts_dir):
    stdout_fn = log_dir + '/create_pickle_settings_file.out'
    stderr_fn = log_dir + '/create_pickle_settings_file.err'
    stderr_fh = open(stderr_fn, 'w')
    cmd = scripts_dir + '/' + 'analysis_settings_file_reader.py'
    print_log_style('Running analysis_settings_file_reader.py..')
    print_log_style('Command: "' + cmd + '"')
    exit_code = subprocess.run([cmd], shell = True, stderr = stderr_fh)
    stderr_fh.close()
    print_log_style('analysis_settings_file_reader.py finished with exit code: ' + str(exit_code.returncode))

    if (exit_code.returncode != 0):
        print_log_style('Error: analysis_settings_file_reader.py program returned a non-zero exit code!')
        sys.exit(1)


def PipeLineSettingsPickleFileReader():
    fh = open('pipeline_run_settings.pickle', 'rb')
    my_settings_dict = pickle.load(fh)
    fh.close()

    return my_settings_dict

def is_a_gz_file(fn):
  if fn[-3:] == '.gz' or fn[-5:] == '.gzip' or fn[-3:] == '.GZ' and fn[-5:] == '.GZIP':
    return True
  else:
    return False 

def fix_mgi_fastq_headers(sample_id, fq1fn, fq2fn):
  print_log_style('About to fix headers in fastq files for sample "' + sample_id + '"..')
  if is_a_gz_file(fq1fn):
    cat_tool = 'zcat'
  else:
    cat_tool = 'cat'
  cmd = cat_tool + " " + fq1fn + '|awk \'{if (NR%4==1) {gsub(/\//, " ");print } else print}\'> ' + sample_id + "_header_corrected_r1.fq"
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'fix_fastq_headers_r1', cmd_line = cmd)

  if is_a_gz_file(fq2fn):
    cat_tool = 'zcat'
  else:
    cat_tool = 'cat'
  cmd = cat_tool + " " + fq2fn + '|awk \'{if (NR%4==1) {gsub(/\//, " ");print } else print}\'> ' + sample_id + "_header_corrected_r2.fq"
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'fix_fastq_headers_r2', cmd_line = cmd)

def  read_primers_info():
  print_log_style('reading primers info..')
  primers_dict = {}
  inv_primers_dict = {}
  
  fh = open(metadata_dir + '/primer_sequences.tsv', 'r')
  
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
  print_log_style('Producing initial abundance table for sample "' + sample_id + '"..')
  cmd = "cat " + sample_id + "_merged.fq.gz|awk \'{if (NR %  4 == 2) print}'|sort|uniq -c| tr -s ' '|tr ' ' $'\t'|cut -c2-|sort -rnk1 >" + "amplicon_abundance_table.tsv"
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'produce_amp_abundance_table', cmd_line = cmd)

def choose_matching_prefix(line, prefixes_dict):
  ret_val = 'Unidentified'
  for k in prefixes_dict.keys():
    v = prefixes_dict[k]
    if line[0:len(k)] == k:
      ret_val = v
      break
  return ret_val

def assign_seqs2primers_ids(inv_primers_dict, sample_id):
  print_log_style('Assigning amplicon-ids according to primers in sample "' + sample_id + '"..')
  ifh = open('amplicon_abundance_table.tsv', 'r')
  ofh = open('amplicon_abundance_table_w_primer_ids.tsv', 'w')
  
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
  print_log_style('About to produce primers distribution table for sample "' + sample_id + '"..')
  cmd = "Rscript " + scripts_dir + '/visualize_primer_dist.R ' + sample_id
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'produce_primers_dist_table_and_pie_plot', cmd_line = cmd)

def cal_percent_mutated_per_grna(sample_id):
  print_log_style('About to calculate percent mutants for sample "' + sample_id + '"..')
  cmd = 'Rscript ' + scripts_dir + '/cal_percent_mutants_for_grnas_regions.R'
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'cal_percent_mutated_per_grna', cmd_line = cmd)

def extract_grna_regions_from_amplicons(sample_id, minseqdepth, flank_len):
  print_log_style('About to extract gRNA region sequences for sample "' + sample_id + '"..')
  cmd = 'Rscript ' + scripts_dir + '/extract_gRNA_region_from_seqs.R ' + sample_id + ' ' + str(minseqdepth) + ' ' + str(flank_len)
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'extract_grna_regions_from_amplicons', cmd_line = cmd)

def count_mutation_rates_per_sample_and_exon(sample_id, has_wt_sample):
  print_log_style('About to count mutation rates for sample "' + sample_id + '"..')
  cmd = 'Rscript ' + scripts_dir + '/count_mutation_rates_per_sample_and_exon_by_intact_gRNA_seq_presence.R' + ' ' + sample_id + ' ' + str(has_wt_sample)
  run_single_job_and_wait(job_id = sample_id, job_general_name = 'count_mutation_rates_per_sample_and_exon', cmd_line = cmd)

def produce_ref_grna_region_seqs_from_wt_sample():
  print_log_style('About to produce refence gRRNA region sequences from wild-type sample..')
  cmd = 'Rscript ' + scripts_dir + '/produce_ref_grna_region_seqs_from_wt_sample.R'
  run_single_job_and_wait(job_id = 'wt', job_general_name = 'produce_ref_grna_region_seqs_from_wt_sample', cmd_line = cmd)

def produce_bases_dist_per_position(sid):
  print_log_style('About to produce gRNA bases distribution per position..')
  cmd = 'Rscript ' + scripts_dir + '/producePerBaseNuclDist.R' + ' ' + sid
  run_single_job_and_wait(job_id = sid, job_general_name = 'produce_gRNA_bases_distribution_per_position', cmd_line = cmd)


def analyze_individual_sample(sample_id, fq1, fq2, minseqdepth, flank_len, has_wt_sample):
  print_log_style('Now analyzing sample "' + sample_id + '"..')
  os.chdir(project_dir)
  os.makedirs(sample_id, exist_ok = True)
  os.chdir(sample_id)
  os.makedirs('log', exist_ok = True)
  #fix_mgi_fastq_headers(sample_id = sample_id, fq1fn = fq1, fq2fn = fq2)
  #merge_r1_r2_reads(sample_id = sample_id)
  #produce_initial_amp_abundance_table(sample_id = sample_id)
  #assign_seqs2primers_ids(inv_primers_dict, sample_id = sample_id)
  #produce_primers_dist_table_and_pie_plot(sample_id = sample_id)
  #extract_grna_regions_from_amplicons(sample_id = sample_id, minseqdepth = int(minseqdepth), flank_len = int(flank_len))
  #count_mutation_rates_per_sample_and_exon(sample_id = sample_id, has_wt_sample = has_wt_sample)
  #cal_percent_mutated_per_grna(sample_id = sample_id)
  if (sample_id == 'wt'):
    produce_ref_grna_region_seqs_from_wt_sample()
  else:
    produce_bases_dist_per_position(sample_id)

def produce_amp_seqs_dist_table_all_samples():
  print_log_style('Producing amplicon sequences distribution table for all samples..')
  cmd = ' cat ampliconDistribution_w_gRNA_region_and_mutation_status_all.tsv|cut -f3,4,5,6,7 > amp_aonly_bundance_table_all_samples_unsorted.tsv && head -n1 amp_aonly_bundance_table_all_samples_unsorted.tsv > amp_aonly_bundance_table_all_samples_header.tsv && tail -n+2 amp_aonly_bundance_table_all_samples_unsorted.tsv | sort -rnk5 | uniq > amp_aonly_bundance_table_all_samples_sorted.tsv && cat amp_aonly_bundance_table_all_samples_header.tsv amp_aonly_bundance_table_all_samples_sorted.tsv > amp_aonly_bundance_table_all_samples.tsv && rm amp_aonly_bundance_table_all_samples_sorted.tsv amp_aonly_bundance_table_all_samples_header.tsv amp_aonly_bundance_table_all_samples_unsorted.tsv'
  run_single_job_and_wait(job_id = 'all', job_general_name = 'produce_amp_only_abundance_table.tsv', cmd_line = cmd)

project_dir = os. getcwd()
make_req_dirs()
scripts_dir = os.path.dirname(os.path.realpath(__file__))
convertExcelFile2Pickle(scripts_dir)
pipeline_settings_dict = PipeLineSettingsPickleFileReader()
primers_dict, inv_primers_dict = read_primers_info()
print(pipeline_settings_dict)
minseqdepth = pipeline_settings_dict['analysis_settings_dict']['min_depth_per_sequence']['Value']
flank_len = pipeline_settings_dict['analysis_settings_dict']['grna_flanking_bases_count']['Value']
if 'wt' in  pipeline_settings_dict['sample_seq_files_dict'].keys():
    has_wt_sample = True
    analyze_individual_sample(sample_id = 'wt', fq1 = pipeline_settings_dict['sample_seq_files_dict']['wt']['r1_filename'], fq2 = pipeline_settings_dict['sample_seq_files_dict']['wt']['r2_filename'], minseqdepth = minseqdepth, flank_len = flank_len, has_wt_sample = has_wt_sample)
else:
    has_wt_sample = False

for sample_id in pipeline_settings_dict['sample_seq_files_dict'].keys():
  if (sample_id == 'wt'):
    continue
  print(sample_id)
  fq1fn = pipeline_settings_dict['sample_seq_files_dict'][sample_id]['r1_filename']
  fq2fn = pipeline_settings_dict['sample_seq_files_dict'][sample_id]['r2_filename']
  analyze_individual_sample(sample_id = sample_id, fq1 = fq1fn, fq2 = fq2fn, minseqdepth = minseqdepth, flank_len = flank_len, has_wt_sample = has_wt_sample)

os.chdir(project_dir)
#producePercentStitchedReport(my_pipeline_settings_dict = pipeline_settings_dict)
#copy_primer_dist_image_files(my_pipeline_settings_dict = pipeline_settings_dict)
#unify_individual_samples_tables()
visualize_mutation_rates()
produce_amp_seqs_dist_table_all_samples()
