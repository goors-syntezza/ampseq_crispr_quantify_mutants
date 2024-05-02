#!/usr/bin/env python3
import xlrd
import sys
import re
import os
import pickle
import csv

log_dir = 'log/'
metadata_dir = '00.Metadata'
expected_sheet_headers = {'sample_seq_files' : ['sample_id', 'r1_filename', 'r2_filename'], 'analysis_settings' : ['Parameter', 'Value'], 'amplicon_sequences' : ['amplicon_id', 'seq'], 'primers' : ['amplicon_id', 'seq'], 'grna' : ['amplicon_id', 'grna_id', 'grna_seq']}
expected_analysis_settings_params = ["grna_flanking_bases_count", "min_depth_per_sequence"]


def convertTable2Dict(wb, ws_name, ncol):
	sheet = wb.sheet_by_name(ws_name)
	if ncol == -1:
		ncol = sheet.ncols
	
	header = [sheet.cell_value(0, i) for i in range(0, ncol)]

	my_dict = {}
	nrows = sheet.nrows
	
	for row_idx in range(1, nrows):
		colCount = 0
		for col_idx in range(0, ncol):
			cell_val = sheet.cell_value(row_idx, col_idx)
			if col_idx == 0:
				my_main_key = str(cell_val)
				my_dict[my_main_key] = {}
			else:
				my_key = str(header[col_idx])
				if not my_key in my_dict[my_main_key].keys():
					my_dict[my_main_key][my_key] = {}
				
				my_dict[my_main_key][my_key] = cell_val
			
	return my_dict

def validate_files_in_range_exist(wb, ws_name, cindex_1, cindex_2, error_msg):

    sheet = wb.sheet_by_name(ws_name)
    vals = []
    for ws_cindex in range(cindex_1, cindex_2 + 1):
        vals = vals + sheet.col_values(ws_cindex)[1 : ]

    for fn in vals:
        fn = str(fn)
        fn_split = fn.split(",")

        for curr_fn in fn_split:
            curr_fn = curr_fn.strip()
            if not os.path.isfile(curr_fn):
                print(error_msg + " (\"" + curr_fn + "\")")
                sys.exit(1)


def validate_unique_values_in_range(wb, ws_name, cindex_1, cindex_2, error_msg):
	sheet = wb.sheet_by_name(ws_name)
	vals = []
	for ws_cindex in range(cindex_1, cindex_2 + 1):
		vals = vals + sheet.col_values(ws_cindex)[1 : ]
	if len(set(vals)) != len(vals):
		print(error_msg)
		sys.exit(1)

def comp_values_in_two_columns(wb, ws1name, ws2name, ws1cindex, ws2cindex, header):
	sheet1 = wb.sheet_by_name(ws1name)
	sheet2 = wb.sheet_by_name(ws2name)
	
	
	val1 = sheet1.col_values(ws1cindex)[1 : ]
	val2 = sheet2.col_values(ws2cindex)[1 : ]

	if len(val1) != len(val2):
		print("Error: " + header + " has different number of items in workbook \"" + ws1name + "\" and \"" + ws2name + "\"!")
		sys.exit(1)
		
	if not set(val1) == set(val2):
		print("Error: " + header + " has different items in workbook \"" + ws1name + "\" and \"" + ws2name + "\"!")
		sys.exit(1)
		
def comp_distinct_values_in_two_columns(wb, ws1name, ws2name, ws1cindex, ws2cindex, header):
    sheet1 = wb.sheet_by_name(ws1name)
    sheet2 = wb.sheet_by_name(ws2name)


    val1 = sheet1.col_values(ws1cindex)[1 : ]
    val2 = sheet2.col_values(ws2cindex)[1 : ]

    if not set(val1) == set(val2):
        print("Error: " + header + " has different items in workbook \"" + ws1name + "\" and \"" + ws2name + "\"!")
        sys.exit(1)


def verify_rows_width(wb, wsn, exp_col_count):
	sheet = wb.sheet_by_name(wsn)
	nrows = sheet.nrows
	if exp_col_count == -1:
		exp_col_count = sheet.ncols
		
	for row_idx in range(1, nrows):
		#print(row_idx)
		colCount = 0
		for col_idx in range(0, exp_col_count):
			cell_val = sheet.cell_value(row_idx, col_idx)
			#print(cell_val)
			if cell_val != '':
				colCount += 1
		#print(colCount)
		if colCount < exp_col_count:
			print("Error: sheet \"" + wsn + "\", row " + str(row_idx) + " has less columns than expected!")
			sys.exit(1)
def verify_sheets_headers(wb, expected_sheet_headers_param):
	for wsn in expected_sheet_headers_param.keys():
		worksheet = workbook.sheet_by_name(wsn)
		actual_ws_header_len = len(worksheet.row(0))
		
		expected_ws_header_len = len(expected_sheet_headers_param[wsn])
		if expected_ws_header_len > actual_ws_header_len:
			print('Error: too short header in sheet "' + wsn + '": header should be ', str(expected_ws_header_len), ' items long, but it\'s', str(actual_ws_header_len))
			sys.exit(1)
		
		for i in range(0, len(expected_sheet_headers_param[wsn])):
			#print(i)
			cell_val = worksheet.cell_value(0, i)
			#print(cell_val)
			if cell_val != expected_sheet_headers_param[wsn][i]:
				print(expected_sheet_headers_param[wsn])
				print('Error: invalid header for column ' + str(i) + ' in sheet "' + wsn + '": "', cell_val, '". Should be "', expected_sheet_headers_param[wsn][i], '"')

def get_sheet_names(wb):
	sns = []
	for sheet in workbook.sheets():
		sns.append(sheet.name)
	return(sns)

def verify_items_in_list(my_list, ref_list):
	miss_items = []
	for item in ref_list:
		if not item in my_list:
			miss_items.append(item)
	return(miss_items)

def verify_legal_name(s, err_message_addition):
	s = str(s)
	is_legal_char_content = bool(re.match("^[A-Za-z0-9_]*$", s))
	if not is_legal_char_content:
		print("Error: invalid name used: \"" + s + "\"")
		print("Valid names must only contain alphanumeric characters and underscore ('_')")
		sys.exit(1)
	if s[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
		print("Error: invalid first character in name: \"" + s[0] + "\"")
		print("Valid name must begin with a characters or underscore ('_') (" + err_message_addition + ")")
		sys.exit(1)


def convertTables2DictsAndWriteToFile(wb):
    sample_seq_files_dict = convertTable2Dict(wb = workbook, ws_name = 'sample_seq_files', ncol = 3)
    analysis_settings_dict = convertTable2Dict(wb = workbook, ws_name = 'analysis_settings', ncol = 2)
    grna_dict = convertTable2Dict(wb = workbook, ws_name = 'grna', ncol = 3)
    amplicon_sequences_dict = convertTable2Dict(wb = workbook, ws_name = 'amplicon_sequences', ncol = 2)
    primers_dict = convertTable2Dict(wb = workbook, ws_name = 'primers', ncol = 2)


    excel_params_dict = {'sample_seq_files_dict' : sample_seq_files_dict, 'analysis_settings_dict' : analysis_settings_dict, 'primers_dict' : primers_dict, 'grna_dict' : grna_dict, 'amplicon_seqs_dict' : amplicon_sequences_dict}

    fh = open('pipeline_run_settings.pickle', 'wb')
    pickle.dump(excel_params_dict, fh)
    fh.close()

    return excel_params_dict
	
def verifyLegalNamesInSheetColumns(wb, sn, col_indices):
	sheet = wb.sheet_by_name(sn)
	header = sheet.row(0)

	nrows = sheet.nrows
	if col_indices == -1:
		col_indices = range(0, sheet.ncols)
	for col_idx in col_indices:
		for row_idx in range(1, nrows):
			cell_val = sheet.cell_value(row_idx, col_idx)
			verify_legal_name(cell_val, 'sheet ' + sn)

def verifySpecificValueExistsInSheetColumns(wb, sn, col_indices, what):
    sheet = wb.sheet_by_name(sn)
    header = sheet.row(0)
    found = False

    nrows = sheet.nrows
    if col_indices == -1:
        col_indices = range(0, sheet.ncols)
    for col_idx in col_indices:
        for row_idx in range(1, nrows):
            cell_val = sheet.cell_value(row_idx, col_idx)
            if cell_val == what:
                found = True
                break
    if not found:
        print("Error: value \"" + what + "\" not found in sheet \"" + sn + "\"!")
        sys.exit(1)

def verifyLegalNamesInSheets(wb):
	verifyLegalNamesInSheetColumns(wb = wb, sn = 'sample_seq_files', col_indices = [0])


def get_cell_value_by_col0_val_and_col_idx(wb, sn, param_name, value_col_idx, required = False):
	sheet = wb.sheet_by_name(sn)
	nrows = sheet.nrows
	my_param_value = None
	for row_idx in range(1, nrows):
		param_Value = sheet.cell_value(row_idx, 0)
		if  param_Value == param_name:
			my_param_value = sheet.cell_value(row_idx, value_col_idx)
	if required and my_param_value is None:
		print("Error: sheet \"" + sn + "\" has a missing value: \"" + param_name + "\"")
		sys.exit(1)
	return my_param_value

def verify_params_in_col0(wb, sn, param_names_as_list, value_col_idx):
	sheet = wb.sheet_by_name(sn)
	nrows = sheet.nrows
	
	for param_name in param_names_as_list:
		my_param_value = None
		for row_idx in range(1, nrows):
			param_Value = sheet.cell_value(row_idx, 0)
			if  param_Value == param_name:
				my_param_value = sheet.cell_value(row_idx, value_col_idx)
		if my_param_value is None:
			print("Error: sheet \"" + sn + "\" has a missing value: \"" + param_name + "\"")
			sys.exit(1)

def validate_analysis_settings_sheet_content(wb):
    verify_params_in_col0(wb = wb, sn = 'analysis_settings', param_names_as_list = expected_analysis_settings_params, value_col_idx = 1)
    
	
    grna_flanking_bases_count = get_cell_value_by_col0_val_and_col_idx(wb, 'analysis_settings', 'grna_flanking_bases_count', 1, required = True)
    grna_flanking_bases_count = int(grna_flanking_bases_count)
    if grna_flanking_bases_count < 0:
        print("Error:  expected 0 < grna_flanking_bases_count")
        sys.exit(1)

    min_depth_per_sequence = get_cell_value_by_col0_val_and_col_idx(wb, 'analysis_settings', 'min_depth_per_sequence', 1, required = True)
    min_depth_per_sequence = int(min_depth_per_sequence)
    if min_depth_per_sequence < 1:
        print("Error:  expected 0 < min_depth_per_sequence")
        sys.exit(1)

    sheet = wb.sheet_by_name('sample_seq_files')

    has_wt_sample = 'No'
    nrows = sheet.nrows
    for row_idx in range(1, nrows):
        sample_id = sheet.cell_value(row_idx, 0)
        if  sample_id == 'wt':
            has_wt_sample = 'Yes'
            break

    return has_wt_sample


def save_table_as_csv(wb, ws_name, out_fn):
	sheet = wb.sheet_by_name(ws_name)

	ncols = sheet.ncols
	nrows = sheet.nrows
	
	my_list = [[''] * ncols for x in range(0, nrows)]
	
	for row_idx in range(0, nrows):
		for col_idx in range(0, ncols):
			cell_val = sheet.cell_value(row_idx, col_idx)
			my_list[row_idx][col_idx] = cell_val
	
	fh = open(out_fn, 'w')
	csvwriter = csv.writer(fh, delimiter = "\t")
	csvwriter.writerows(my_list)
	fh.close()

def write_groups2comp(wb):
	sheet = wb.sheet_by_name('samples_metadata')
	ncols = sheet.ncols
	
	my_list = ['']* (ncols -1)
	
	for col_idx in range(1, ncols):
		cell_val = sheet.cell_value(0, col_idx)
		my_list[col_idx - 1] = cell_val
	
	out_fn = metadata_dir + 'groups_to_compare.tsv'
	fh = open(out_fn, 'w')
	for line in my_list:
		fh.write(line + "\n")
	fh.close()
	

def write_genome_name2file(excel_params_dict):
	genome_name = excel_params_dict['genome_info_dict']['genome_name']['Value']
	
	out_fn = metadata_dir + 'genome_name.txt'
	fh = open(out_fn, 'w')
	fh.write(genome_name + "\n")
	fh.close()

def write_qualimap_info_filename2file(excel_params_dict):
    qualimap_info_file_name = excel_params_dict['analysis_settings_dict']['qualimap_info_file']['Value']

    out_fn = metadata_dir + 'qualimap_info_filename.txt'
    fh = open(out_fn, 'w')
    fh.write(qualimap_info_file_name + "\n")
    fh.close()
	

workbook = xlrd.open_workbook("pipeline_settings.xlsx")
expected_sheet_names = expected_sheet_headers.keys()
sheet_names = get_sheet_names(workbook)
miss_sheets = verify_items_in_list(sheet_names, expected_sheet_names)
if (len(miss_sheets) > 0):
	print("Error: missing sheets. Missing sheet names: ", ",".join(miss_sheets))
	sys.exit(1)

verify_sheets_headers(workbook, expected_sheet_headers)
#worksheet = workbook.sheet_by_index(0)
verify_rows_width(workbook, 'sample_seq_files', 3)
verify_rows_width(workbook, 'analysis_settings', 2)
verify_rows_width(workbook, 'amplicon_sequences', 2)
verify_rows_width(workbook, 'primers', 2)
verify_rows_width(workbook, 'grna', 3)
verifyLegalNamesInSheets(wb = workbook)
has_wt_sample = validate_analysis_settings_sheet_content(wb = workbook)
comp_values_in_two_columns(workbook, 'amplicon_sequences', 'primers', 0, 0, 'amplicon_id')
comp_distinct_values_in_two_columns(workbook, 'primers', 'amplicon_sequences', 0, 0, 'amplicon_id')

validate_unique_values_in_range(wb = workbook, ws_name = 'sample_seq_files', cindex_1 = 1, cindex_2 = 2, error_msg = 'Error: duplicate fastq files!')
#verifySpecificValueExistsInSheetColumns(wb = workbook, sn = 'sample_seq_files', col_indices = [0], what = 'wt')
excel_params_dict = convertTables2DictsAndWriteToFile(wb = workbook)
#print(excel_params_dict)
excel_params_dict['analysis_settings_dict']['has_wt_sample'] = {'Value' : has_wt_sample}
save_table_as_csv(wb = workbook, ws_name = 'amplicon_sequences', out_fn = metadata_dir + '/amplicon_sequences.tsv')
save_table_as_csv(wb = workbook, ws_name = 'grna', out_fn = metadata_dir + '/grna_sequences.tsv')
save_table_as_csv(wb = workbook, ws_name = 'primers', out_fn = metadata_dir + '/primer_sequences.tsv')
save_table_as_csv(wb = workbook, ws_name = 'analysis_settings', out_fn = metadata_dir + '/analysis_settings.tsv')

