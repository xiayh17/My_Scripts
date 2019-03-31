import os,sys
import glob
import pandas as pd
import io


class AutoVivification(dict):             # my_dict[1][2] = 3   ==> {1: {2: 3}}
	def __getitem__(self, item):
		try:
			return dict.__getitem__(self, item)
		except KeyError:
			value = self[item] = type(self)()
			return value

def get_fqStat(fqStat):
#(total_base_1,total_base_2,Q20_1,Q20_2,Q30_1,Q30_2,GC_1,GC_2,err_1,err_2) = 0
	(total_base,Q20,Q30,GC,err,AT_div,GC_div) = (0,0,0,0,0,0,0)
	(base_A,base_C,base_G,base_T,base_N) = (0,0,0,0,0)
	(a_minus_t,c_minus_g,icnt) = (0,0,0)
	with open(fqStat,'r') as fh:
		for line in fh:
			line = line.rstrip()
			lines = line.split('\t')
			if line.startswith('#'):
				if line.startswith('#BaseNum'):
					total_base = lines[1]
				elif line.startswith('#GC%'):
					GC = lines[1]
				elif line.startswith('#Q20%'):
					Q20 = lines[1]
				elif line.startswith('#Q30%'):
					Q30 = lines[1]
				elif line.startswith('#EstErr%'):
					err = float(lines[1])
					err = "%.2f" % err
			elif float(lines[0]) >=11:
				(base_A,base_C,base_G,base_T,base_N) = lines[1:6]
				base_total = int(base_A) + int(base_T) + int(base_C) + int(base_G) + int(base_N)
				a_minus_t += float(abs(int(base_A) - int(base_T))/base_total)
				c_minus_g += float(abs(int(base_C) - int(base_G))/base_total)
				icnt+=1
	
	AT_div = float(a_minus_t/icnt*100)
	GC_div = float(c_minus_g/icnt*100)
	AT_div = "%.2f" % AT_div
	GC_div = "%.2f" % GC_div
	return(total_base,Q20,Q30,GC,err,AT_div,GC_div)


def get_lib(samplelist):
	all_dicts = AutoVivification()

	with open(samplelist, 'r') as fh:
		for line in fh:
			lines = line.rstrip().split('\t')
			sample_name,p_lib = lines[:2]
			paths = lines[2:]

			lib_name = p_lib.split('-')[0]
			barcodes = p_lib.split('-')[1:]
			lib_barcodes = ['-'.join([lib_name,i]) for i in barcodes]
			if not all_dicts[sample_name]['lib_barcodes']:
				all_dicts[sample_name]['lib_barcodes'] = []
			all_dicts[sample_name]['lib_barcodes'].extend(lib_barcodes)
			all_dicts[sample_name]['lib_barcodes'] = list(set(all_dicts[sample_name]['lib_barcodes']))

			if not all_dicts[sample_name]['paths']:
				all_dicts[sample_name]['paths'] = []
			all_dicts[sample_name]['paths'].extend(paths)
			all_dicts[sample_name]['paths'] = list(set(all_dicts[sample_name]['paths']))

	return(all_dicts)

def get_all_AT(all_dicts):
	all_AT = []
	for sample_name in all_dicts:
		for lib_barcode in all_dicts[sample_name]['lib_barcodes']:
			lib = lib_barcode.split('-')[0]
			for path in all_dicts[sample_name]['paths']:
				if os.path.exists(path):
					fqStat1,fqStat2 = (0,0)
					dirs = os.listdir(path)
					barcode_dirs = [os.path.join(path, d) for d in dirs if lib_barcode in d]
					#print(barcode_dirs)
					for barcode_dir in barcode_dirs:
						for stat_file in os.listdir(barcode_dir):
							if stat_file.endswith('_1.fq.fqStat.txt'):
								fqStat1 = os.path.abspath(os.path.join(barcode_dir,stat_file))
							elif stat_file.endswith('_2.fq.fqStat.txt'):
								fqStat2 = os.path.abspath(os.path.join(barcode_dir,stat_file))
						#print(fqStat1,fqStat2)
						qc1 = get_fqStat(fqStat1)
						qc2 = get_fqStat(fqStat2)
						qc_1_2 = list(zip(qc1,qc2))
						qc_str = [i+";"+j for i,j in qc_1_2]
						format_str = "%s\t%s\t%s" + len(qc_str)*"\t%s"
						#print(format_str %(sample_name,lib,lib_barcode,*qc_str))
						all_AT.append(format_str %(sample_name,lib,lib_barcode,*qc_str))
	return all_AT

def split_col(dat,sep=';'):
	colnames = list(dat.columns)
	#print(colnames)
	for col in colnames[3:]:
		col_1 = "%s_1" % col
		col_2 = "%s_2" % col
		dat[col_1],dat[col_2] = dat[col].str.split(sep,1).str				#split one column to two
		dat[[col_1,col_2]] = dat[[col_1,col_2]].apply(pd.to_numeric)                    #transform str to numeric
		dat.drop(col, axis=1, inplace=True)                                             #remove original column
	return(dat)


def ave_all(groups, columns_ave, df):
	columns_ave_res = []
	data_size = df.groupby(groups).apply(lambda x:(x['data_size_1']*2).sum())
	data_size.name = "data_size"
	columns_ave_res.append(data_size)
	for col in columns_ave:                                                                 # avrage of columns one by one 
		func = lambda x:(x['data_size_1']/x['data_size_1'].sum()*x[col]).sum()          #lambda fuction
		df_tmp = df.groupby(groups).apply(func)                                         #groupby then apply function, get Series
		df_tmp.name = col                                                               #name the Series by df column name
		columns_ave_res.append(df_tmp)                                                  #reserve it 
	df_ave = pd.concat(columns_ave_res,axis=1)                                              #merge all Series as dataframe columns
	return(df_ave)


if __name__ == '__main__':
	all_dicts = get_lib(sys.argv[1])
	outfile = sys.argv[2]

	all_AT = get_all_AT(all_dicts)
	col = ['sample','lib','barcode','data_size','Q20','Q30','GC','error_rate','AT_div','GC_div']
	col = '\t'.join(col)
	all_AT = [col] + all_AT
	#print(all_AT)
	df = pd.read_csv(io.StringIO('\n'.join(all_AT)), delim_whitespace=True)
	print(df)
	df2 = split_col(df)
	print(df2)

	df2.drop('data_size_2',axis=1,inplace=True)
	groups = list(df.columns)[:2]              #'sample','lib'
	columns_ave = list(df.columns)[4:]         #'Q20','Q30','GC','error_rate','AT_div''GC_div'
	df3 = ave_all(groups,columns_ave,df2)
	df3.to_csv(outfile,sep='\t',float_format='%.2f')
