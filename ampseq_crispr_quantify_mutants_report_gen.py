# -*- coding: utf-8 -*-

"""
Spyder Editor

This is a temporary script file.
"""

#%% import
import os
import shutil
import pandas as pd
import numpy as np
import glob
from pathlib import Path 
from detect_delimiter import detect 
import uuid
import re
import sys
import argparse
from datetime import datetime
log_dir = os.path.abspath('log/')
import subprocess
import gzip
import unicodedata
from datetime import date

def copyFileHead(sfn, tfn, nLines):
  ifh = open(sfn, 'r')
  ofh = open(tfn, 'w')
  
  lineCount = 0
  
  for line in ifh:
    lineCount += 1
    if lineCount <= nLines:
      ofh.write(line)
    else:
      break
  
  ifh.close()
  ofh.close()

def cutHeaderLinesFromFile(sfn, tfn, lineSkipCount):
  ifh = open(sfn, 'r')
  ofh = open(tfn, 'w')
  
  lineCount = 0
  
  for line in ifh:
    lineCount += 1
    if lineCount > lineSkipCount:
      ofh.write(line)
  
  ifh.close()
  ofh.close()

def get_used_sw_vers():
    run_single_job_and_wait(job_id = '_', job_general_name = 'get_used_sw_vers', cmd_line = 'python3 ' + scripts_dir + '/identify_versions/identify_versions.py')

def replaceHTMLInjections(line):
    injectionID = line[line.find('INJECTABLE_NAME'):line.rfind('INJECTABLE_NAME')-1].replace('INJECTABLE_NAME[', '')
    if injectionID == 'place_date':
        cur_date = str(date.today().strftime("%d-%b-%Y"))
        new_line = line[:line.find('INJECTABLE_NAME')] + cur_date + line[line.rfind('INJECTABLE_NAME')+15:]
        return new_line
    elif injectionID == 'place_project_name':
        new_line = line[:line.find('INJECTABLE_NAME')] + cmdline_args.name + line[line.rfind('INJECTABLE_NAME')+15:]
        return new_line
    elif injectionID == 'place_host_genome_ref':
        new_line = line[:line.find('INJECTABLE_NAME')] + cmdline_args.Host_genome_ref + line[line.rfind('INJECTABLE_NAME')+15:]
    elif injectionID == 'place_host_viral_ref':
        new_line = line[:line.find('INJECTABLE_NAME')] + cmdline_args.ref_viral + line[line.rfind('INJECTABLE_NAME')+15:]
    else:
        print("Error: unknown injectionID parameter value: \"" + injectionID  + "\"")
        sys.exit(1)

def readPerformedAnalyses():
    my_dict = {}
    log_dir = 'log/'
    performed_analyses_fn = log_dir + 'performed_analyses.txt'
    fh = open(performed_analyses_fn, 'r')
    for line in fh:
        my_dict[line.strip()] = 1
    fh.close()

    return my_dict

def ext_analysis_id(line):
    analysis_id = line[line.rfind('[') + 1:line.rfind(']')]
    return analysis_id


def ext_reference_text_and_replace_with_num(line, num):
	ref_text = line[line.find('[ref:') + 5:line.find(':ref]')]
	text2replace = line[line.find('[ref:') + 1:line.find(':ref]')+4]
	line = line.replace(text2replace, str(num))
	return line, ref_text

def print_log_style(s):
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(str(dt_string) + '|' + s)


def run_single_job_and_wait(job_id, job_general_name, cmd_line):
    print_log_style('run_single_job_and_wait - About to launch procedure "' + job_general_name  + '" with job id "' + job_id + '"')
    print_log_style('run_single_job_and_wait - Executing the command: "' + cmd_line  + '"')
    out_fh = open(log_dir + '/' + job_general_name + '_' + job_id + '.out', 'w')
    err_fh = open(log_dir + '/' + job_general_name + '_' + job_id + '.err', 'w')

    exit_code = subprocess.call(cmd_line, bufsize=- 1, executable = '/bin/bash', stdin=None, stdout=out_fh, stderr = err_fh, preexec_fn=None, close_fds=True, shell=True, cwd=None, env=None, universal_newlines=None, startupinfo=None, creationflags=0, restore_signals=True, start_new_session=False)

    out_fh.close()
    err_fh.close()

    if exit_code != 0:
        print_log_style('run_single_job_and_wait - Error procedure "' + job_general_name + '" with job-id "' + job_id + '" returned a non-zero return code: ' + str(exit_code))
        sys.exit(1)
    else:
        print_log_style('run_single_job_and_wait - Procedure "' + job_general_name + '" with job-id "' + job_id + '" returned 0 (normal completetion)')


def processCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, help="Name for the report", required = True)
    args = parser.parse_args()

    return args

def copyFileHeader(sfn, dfn, nlines):
    if sfn[-2 : ] == 'gz':
        my_tool = 'zcat'
    else:
        my_tool = 'cat'
	#cmd = 'set -e;set -o pipefail;' + my_tool + ' ' + sfn + ' | head -n' + str(nlines) + ' > ' + dfn
	#run_single_job_and_wait(job_id = '0', job_general_name = 'Copy file header', cmd_line = cmd)

    if my_tool == 'zcat':
        fh = gzip.open(sfn, 'rt')
        print('!')
    else:
        fh = open(sfn, 'r')

    ofh = open(dfn, 'w')
    
    lc = 0


    for line in fh:
        #if my_tool == 'zcat':
            #line = unicodedata.normalize('', line).encode('ascii', 'ignore')
        ofh.write(line)
        lc += 1
        if lc >= nlines:
            break



    fh.close()
    ofh.close()

def table_section_write(table, all_or_top = 'all', pageSize = 15):    
    
    
    # if not 'all' keep only top # lines
    if not all_or_top == 'all' and not type(all_or_top) == int:
        raise Exception("Sorry, 'all_or_top' can be 'all' or a number of top lines to be printed")
        
    elif type(all_or_top) == int:
        table = table[:all_or_top]
    
    
    ### write html section
    
    # for table of more then pageSize (num of lines), add page select nemu (pagination)
    if table.shape[0] > pageSize:
        
        # unique id
        id = str(uuid.uuid4()).split("-", 1)[0]
        
        # Javascript for table pagination 
        with open(report_temp_dir+'/table_pagination.html', "r", encoding='utf-8') as file:
            html_table = file.read()
        
        # set id for uniqueid variables
        html_table = html_table.replace('uniqueid', id)
        
        # reset page size (numbe of lines)
        html_table = html_table.replace('var pageSize_'+id+' = 15;', 'var pageSize_'+id+' = '+str(pageSize)+';')
        
        # table section start 
        html_table += '\n\n<p class="center">\n<div style="overflow-x:auto">\n<table class="gy" id="mytable_'+id+'">\n<thead><tr>\n'
        
        html_table += '\n		<a id="btn0_'+id+'"></a>'
        html_table += '\n		<a id="sjzl_'+id+'"></a>'
        html_table += '\n		<a  href="javascript:void(0)" id="btn1_'+id+'"> | First | </a>'
        html_table += '\n		<a  href="javascript:void(0)" id="btn2_'+id+'"> Previous | </a>'
        html_table += '\n		<a  href="javascript:void(0)" id="btn3_'+id+'"> Next | </a>'
        html_table += '\n		<a  href="javascript:void(0)" id="btn4_'+id+'"> Last | </a>'
        html_table += '\n		<a>Go to </a>'
        html_table += '\n		<input id="changePage_'+id+'" type="text" size="1" maxlength="4"/>'
        html_table += '\n		<a> page </a>'
        html_table += '\n		<a  href="javascript:void(0)" id="btn5_'+id+'"> Jump </a>\n'
        
    else:
        # table section start without pagination
        html_table = '<p class="center">\n<div style="overflow-x:auto">\n<table class="gy">\n<thead><tr>\n'
    
    
    # add headers
    colnames = np.insert(table.columns.values, 0, table.index.name)
    
    for column_name in colnames:
        html_table += '<th style="text-align:center">'+str(column_name)+'</th>'
    
    html_table += "\n</tr></thead>"
    
    
    # add rows
    for idx, row in table.iterrows():
        
        # row opening
        html_table += '\n<tr><td style="text-align:left">'+str(row.name)+'</td>'
        
        # insert row values
        for value in row:
            html_table += '<td>'+str(value)+'</td>'
        
        # row closing: replase last value of the row to last row value format 
        html_table = html_table[:html_table.rfind('<td>'+str(value)+'</td>')] + '<td rowspan=1>'+str(value)+'</td></tr>'
    
    html_table += "\n</table>\n</div>\n</p>"
    
    
    # write to html file
    report_file.writelines(html_table)


#%% constract and write boxSelect section

def boxSelect_section_write(line):
    
    # get location and name formats
    subdir, prefix, suffix = line[line.rfind('[')+1:line.rfind(']')].split(',')
    
    # boxSelect section ids (each section needs to have it's owen unique names)
    id = str(uuid.uuid4()) # unique id
    boxSelect_id = 'boxSelect-'+id
    imageToSwap_id = 'imageToSwap-'+id
    
    
    ## write html section
    
    # boxSelect select section start 
    boxSelect = '<html>\n<select id="'+boxSelect_id+'">'
    
    # add plots to select
    paths = glob.glob(subdir+"/"+prefix+"*"+suffix)
    for path in paths:
        name = os.path.basename(path).replace(suffix, "", 1).replace(prefix, "", 1)
        
        # plots to select from
        boxSelect += '\n<option value="'+path+'">'+name+'</option>'
    
    # boxSelect select section end, front img and functions    
    print(line)
    boxSelect += '\n</select>\n<center><img name="'+imageToSwap_id+'" src="' + paths[0] + '"/></center>'
    boxSelect += '\n<script>$(document).ready(function(){'
    boxSelect += '\n                 $("#'+boxSelect_id+'").change(function(){'
    boxSelect += '\n                 $("img[name='+imageToSwap_id+']").attr("src",$(this).val());'
    boxSelect += '\n                 });'
    boxSelect += '\n  });</script>\n</html>'
    
    # write boxSelect section to report file  
    report_file.writelines(boxSelect)

#%% replace in line, "VAR_NAME['some_path/*/and_some_more*.suffix']VAR_NAME" with the variable name in place of the first *, and write to html



def var_name_write(line):
    
    # pattern for files search
    pattern = report_dir+line[line.find('VAR_NAME'):line.rfind('VAR_NAME')-1].replace('VAR_NAME[', '')
    
    # all paths matching
    paths = glob.glob(pattern)
    
    # first file name in list as base for var_name
    var_name = paths[0]
    
    # prefix of all files to remove from var_name
    prefix = os.path.commonprefix(paths+[pattern])
    var_name = var_name.replace(prefix, '')
    
    # suffix
    suffix = pattern[pattern.find('*')+1:]
    
    # prefix of suffix to remove
    pre_suffix = suffix[:suffix.find('*')]
    
    # remove suffix (from its prefix)
    var_name = var_name[:var_name.find(pre_suffix)]    
    
    # replace in line
    new_line = line[:line.find('VAR_NAME')] + var_name + line[line.rfind('VAR_NAME')+8:]
    
    # write boxSelect section to report file  
    report_file.writelines(new_line)
    

#%% replace in line, "EXMPL_PATH['some_path/*/and_some_more*.suffix']EXMPL_PATH" with the first path that match, and write to html

def exmpl_path_write(line):
    # pattern for files search
    pattern = report_dir+line[line.find('EXMPL_PATH'):line.rfind('EXMPL_PATH')-1].replace('EXMPL_PATH[', '')
    
    # all paths matching
    paths = glob.glob(pattern)
    
    # first file name in list as base for var_name
    print(line)
    print(paths)
    first_path = paths[0].replace(report_dir, '')
    
    # replace in line
    new_line = line[:line.find('EXMPL_PATH')] + first_path + line[line.rfind('EXMPL_PATH')+10:]
    
    # write boxSelect section to report file  
    report_file.writelines(new_line)

def writeMethods(performedAnalysesList):
    performedAnalysesList = list(performedAnalysesList.keys())
    print(performedAnalysesList)
    inMethodsFN = '/home/goors/Scripts/mRNA_medical_v2/Report/methods.html'
    outMethodsFN  = project_dir + '/01.Metadata/methods.html'
    outMethodsFNInReport = project_dir + '/04.Report/src/data_for_report/methods.html'

    ifh = open(inMethodsFN, 'r')
    ofh = open(outMethodsFN, 'w')

    shouldWriteChapter = True
    chapterIndex = 0
    referenceIndex = 0

    ref_list = []

    for line in ifh:
        if 'ADDCHAPTER[[' in line:
            newChapterID = line[line.find('ADDCHAPTER[[') + 12:line.find(']') ]
            print('"' + newChapterID + '"')
            if newChapterID in performedAnalysesList:
                shouldWriteChapter = True
                chapterIndex += 1
            else:
                shouldWriteChapter = False
            print(shouldWriteChapter)
            continue
        if not shouldWriteChapter:
            continue
        if 'CHAPTER_ID[[]]' in line:
            chapter_id_indices = line.replace('CHAPTER_ID[[]]', str(chapterIndex))
        if 'REF[[' in line:
            ref_coords = [(m.start(), m.end()) for m in re.finditer('REF\[\[(.*?)\]\]', line)]
            #for ref_coord in ref_coords:
            while 'REF[[' in line:
                ref_coords = [(m.start(), m.end()) for m in re.finditer('REF\[\[(.*?)\]\]', line)]
                ref_coord = ref_coords[0]
                referenceIndex += 1
                text_start = ref_coord[0] + 5
                text_end = ref_coord[1] - 2
                ref_text = line[text_start : text_end]
                line = line[0 : ref_coord[0]] + '[' + str(referenceIndex) + ']' + line [ref_coord[1]  :]
                ref_list.append(ref_text)
        ofh.write(line)
    i = 0
    for ref in ref_list:
        i += 1
        line = '<p>' + str(i) + '. ' + ref + '</p>\n'
        ofh.write(line)
    ofh.write('</ol>\n')
    ofh.write('</body>\n</html>\n')
    ifh.close()
    ofh.close()

    shutil.copy(outMethodsFN, outMethodsFNInReport)

def write_file_content_as_is(line):
    content_fn = line[line.find('[[') + 2:line.rfind(']]')].replace('EXMPL_PATH[', '')
    
    content_fh = open(content_fn, 'r')
    
    for content_file_line in content_fh:
      report_file.writelines(content_file_line)
    
    content_fh.close()


def read_tsv_file_as_dict(fn):
  #metadata_dir + '/ref_dbs_in_use.tsv'
  fh = open(fn, 'r')
  my_dict = {}
  for line in fh:
    line = line.strip().split("\t")
    my_dict[line[0]] = line[1]
  fh.close()
  return my_dict

def copy_results_files():
  shutil.copy('../00.Metadata/primer_sequences.tsv', '.')
  shutil.copy('../00.Metadata/grna_sequences.tsv', '.')
  shutil.copy('../00.Metadata/analysis_settings.tsv', '.')
  shutil.copy('../00.Metadata/amplicon_sequences.tsv', '.')
  shutil.copy('../reads_merge_statistics.tsv', '.')
  shutil.copy('../grna_mutation_rates_and_counts_per_amplicon_all_samples_all.tsv', '.')
  shutil.copy('../ampliconDistribution_w_gRNA_region_and_mutation_status_all.tsv', '.')
  shutil.copy('../mutation_rate_and_reads_count_per_sample_and_exon_pie_lattice_based_on_intact_gRNA_region_only_mutated_samples.png', '.')
  shutil.copytree('../plots', './plots', dirs_exist_ok = True)
  copyFileHead(sfn = 'ampliconDistribution_w_gRNA_region_and_mutation_status_all.tsv', tfn = 'ampliconDistribution_w_gRNA_region_and_mutation_status_all_head.tsv', nLines = 25)
  shutil.copy('../amp_aonly_bundance_table_all_samples.tsv', '.')
  copyFileHead(sfn = 'amp_aonly_bundance_table_all_samples.tsv', tfn = 'amp_aonly_bundance_table_all_samples_head.tsv', nLines = 25)
  
  for file in glob.glob(r'../*/*gRNAbaseCountPerPosition*.tsv'):
    shutil.copy(file, '.')
  
  

scripts_dir = os.path.dirname(os.path.realpath(__file__))
metadata_dir = '00.Metadata/'

#%% set directories


cmdline_args = processCommandLine()

#ref_dbs_used_dict = read_tsv_file_as_dict(metadata_dir + '/ref_dbs_in_use.tsv')

project_dir = os.path.abspath(os.getcwd()) # project_dir = "/home/stav/Projects/Syntezza/Goor"
#report_temp_dir = os.path.expanduser("~/Scripts/mRNA_medical_v2/Report/")
report_temp_dir = scripts_dir

# create report directory
report_dir = os.path.abspath(project_dir+"/02.Report/")
if not os.path.exists(report_dir): os.mkdir(report_dir)




# copy src directory
shutil.copytree(report_temp_dir+"/src", report_dir+"/src", symlinks=False, ignore=None, ignore_dangling_symlinks=False, dirs_exist_ok=True)




    
#%%  write html report

# change to report directory
#performed_analyses_dict = readPerformedAnalyses()
performed_analyses_dict = {}
performed_analyses_dict['Overview'] = 1
performed_analyses_dict['run_parameters'] = 1
performed_analyses_dict['insert_seqs_statistics'] = 1
performed_analyses_dict['insert_seqs_count'] = 1
performed_analyses_dict['insert_seqs_dist_plot'] = 1
performed_analyses_dict['amp_seqs_count'] = 1
performed_analyses_dict['bases_dist_per_pos'] = 1
os.chdir(report_dir)
copy_results_files()



template_file = open(report_temp_dir+'/report_temp_r2.html',"r", encoding='utf-8')
report_file = open(report_dir+'/Analysis_report.html',"w", encoding='utf-8')

expectedHeader = 'Syntezza RNASeq virus detection pipeline'


shouldWrite = True
chapter_num = 0
reference_num = 0
writeCurrentChapter = True



for line in template_file.readlines():

    if 'Add chapter here' in line:
        chapter_id = ext_analysis_id(line)
        print(chapter_id)
        if chapter_id in performed_analyses_dict:
            if not 'no_increase' in line:
                chapter_num += 1
            writeCurrentChapter = True
            print('Writing this chapter')
        else:
            writeCurrentChapter = False
            print('Not Writing this chapter')



    if writeCurrentChapter == False:
        continue


    if 'INJECTABLE_NAME' in line:
        line = replaceHTMLInjections(line)

    if '[ChapterNumber]' in line:
        line = line.replace('[ChapterNumber]', str(chapter_num))


    if 'Add gene count table here' in line:
        ## open gene_counts
        gene_counts = pd.read_csv(report_dir+"/src/data_for_report/fpkm_table_with_gene_info.tsv", sep="\t", header=0, index_col=0)
        # sorted by row sum
        gene_counts = gene_counts.loc[gene_counts.drop(['gene_name', 'description'], axis=1).sum(axis=1).sort_values(ascending=False).index]
        # write table
        table_section_write(gene_counts, all_or_top = 50)
        
        
    elif 'Add as-is table here' in line:
        
        ## in case of multiple paths matching, the first will be written
        print(line)
        filename = glob.glob(report_dir+line[line.rfind('[')+1:line.rfind(']')])[0]
        
        # detect delimiter
        with open(filename) as f:
            firstline = f.readline()
            delimiter = detect(firstline)
        print('Delimiter: "' + delimiter + '"')    
        if filename[-3 : ] == 'txt':
            delimiter = '\t'
        # write table
        print(filename)
        if delimiter == "\t":
            table = pd.read_csv(filename, sep = delimiter, header=0, index_col=0, engine = 'python', quoting = 3, quotechar = '"')
        else:
            table = pd.read_csv(filename, sep = delimiter, header=0, index_col=0, engine = 'python', quotechar = '"')
        if 'JCEC' in filename:
            all_or_top_val = 25
        else:
            all_or_top_val = 'all' 
        table_section_write(table, all_or_top = all_or_top_val)
        
        
    elif 'Add boxSelect section here' in line:
        print(line)
        boxSelect_section_write(line)
    
    elif 'VAR_NAME' in line:
        var_name_write(line)
        
    elif 'EXMPL_PATH' in line:
        exmpl_path_write(line)
        
    elif 'Add as-is content from file here' in line:
        write_file_content_as_is(line)

    
    else: report_file.writelines(line)

template_file.close()
report_file.close()
