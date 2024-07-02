#!/usr/bin/python3
import sys,os,subprocess,time,argparse

##########################################################################

parser = argparse.ArgumentParser(description="Tool for mammalian U2 branch point prediction.")

parser.add_argument("-i", "--input",
                    dest="fname",
                    action="store",
                    required=True,
                    help="Input file name. "
                         "FASTA format only.")

parser.add_argument("-s", "--species",
                    dest="species",
                    action="store",
                    type=str,
                    choices=['Hsap', 'Ptro', 'Mmul', 'Mmus', 'Rnor', 'Cfam', 'Btau'],
                    required=True,
                    help="Species under study. Case sensitive. "
                         "Options: "
                         "Hsap (Homo sapiens) "
                         "Ptro (Pan troglodytes) "
                         "Mmul (Macaca mulatta) "
                         "Mmus (Mus musculus) "
                         "Rnor (Rattus norvegicus) "
                         "Cfam (Canis familiaris) "
                         "Btau (Bos taurus).")

parser.add_argument("-l", "--max-len",
                    dest="maxslen",
                    action="store",
                    type=int,
                    default=100,
                    help="Number of bases at the 3' end of the input sequences that will be scanned. "
		    "SVM-BPfinder assumes they correspond to the 3' end of introns. "
		    "For sequences of length smaller than this value, the entire sequence will be scanned. Default value is 100.")

parser.add_argument("-d", "--min-dist",
                    dest="mindist3ss",
                    action="store",
                    type=int,
                    default=15,
                    help="Distance in nucleotides allowed between the branch point A and the 3' splice-site. Default value is 15.")

options = parser.parse_args()

##########################################################################

base_dir = os.path.dirname(os.path.realpath(__file__))

model=base_dir+'/MODELS/'+options.species+'.svm.model'
infile="%d"%(time.time()*100)+'.in'
outfile="%d"%(time.time()*100)+'.out'
gfscript=base_dir+'/SCRIPTS/svm_getfeat.py'
svmlclass=base_dir+'/SCRIPTS/svm_classify'


##########################################################################

def print_merge(ffn,sfn):
	print('\t'.join(['seq_id','agez','ss_dist','bp_seq','bp_scr','y_cont','ppt_off','ppt_len','ppt_scr','svm_scr']))
	for l in subprocess.getoutput('paste '+ffn+' '+sfn).split('\n'):
		if l!='':
			la=l.split(' #') 
			feats=[x[2:] for x in la[0].split(' ')[1:]]
			extras=la[1].split('\t')
			pa=extras[:-2]+feats[:-1]+[extras[-2]]+[feats[-1]]+[extras[-1]]
			print('\t'.join(pa))		

##########################################################################

os.system(gfscript+' '+options.fname+' '+options.species+' '+str(options.maxslen)+' '+str(options.mindist3ss)+' > '+infile)
os.system(svmlclass+' '+infile+' '+model+' '+outfile+' > /dev/null') 
print_merge(infile,outfile)
os.remove(infile)
os.remove(outfile)

