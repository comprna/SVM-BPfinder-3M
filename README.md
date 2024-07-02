# SVM-BPfinder  - A tool for mammalian BP prediction. 

(code also available from https://bitbucket.org/regulatorygenomicsupf/svm-bpfinder)

Andre Corvelo and Eduardo Eyras | Regulatory Genomics @ Universitat Pompeu Fabra, Barcelona, Spain | 2010 

This tool is free for purpose of academic, non-commercial research. The software must not be further distributed without prior permission of the authors.

### CONTACTS:
``` 
acorvelo[at]gmail.com
eduardo.eyras[at]anu.edu.au
```
   
### REFERENCE:
```
A. Corvelo, M. Hallegger, C.W.J. Smith, E. Eyras. (2010). Genome-wide Association between Branch Point Properties and Alternative Splicing
```
------------------------------------------


### IMPORTANT:

SVM-BPfinder requires SVMlight, which can be downloaded at:
```
http://osmot.cs.cornell.edu/svm_light/current/
```

After downloading SVMlight, copy the executable 'svm_classify' to the ./SCRIPTS/ folder.

Make sure you have permission to execute these three scripts:
1) 'svm_bpfinder.py'
2) 'SCRIPTS/svm_getfeat.py'
3) 'SCRIPTS/svm_classify'  

To know more about SVMlight, please visit http://svmlight.joachims.org/

------------------------------------------

### USAGE:
```
svm_bpfinder.py [-h] -i FNAME -s {Hsap,Ptro,Mmul,Mmus,Rnor,Cfam,Btau}
                [-l MAXSLEN] [-d MINDIST3SS]
```

### ARGUMENTS:
```
  -h, --help            show this help message and exit

  -i FNAME, --input FNAME
                        Input file name. (Multi)FASTA format only.

  -s {Hsap,Ptro,Mmul,Mmus,Rnor,Cfam,Btau}, --species {Hsap,Ptro,Mmul,Mmus,Rnor,Cfam,Btau}
                        Species under study. Case sensitive.
                        Options: Hsap (Homo sapiens) Ptro (Pan troglodytes)
                        Mmul (Macaca mulatta) Mmus (Mus musculus) Rnor (Rattus
                        norvegicus) Cfam (Canis familiaris) Btau (Bos taurus).

  -l MAXSLEN, --max-len MAXSLEN
                        Number of bases at the 3' end of the input sequences that
                        will be scanned. SVM-BPfinder assumes they correspond
                        to the 3' end of introns. For sequences of length
                        smaller than this value, the entire sequence will be
                        scanned. Default value is 100.

  -d MINDIST3SS, --min-dist MINDIST3SS
                        Distance in nucleotides allowed between the branch
                        point A and the 3' splice-site. Default value is 15.

```

### OUTPUT:

Results are printed to STDOUT, tab delimited, one line per BP candidate. Header included.

Output fields:
```
seq_id  - Sequence Identifier
agez    - AG dinucleotide Exclusion Zone length
ss_dist - Distance to 3' splice-site
bp_seq  - BP sequence (9-mer; from -5 to +3 relative to the BP adenine)
bp_scr  - BP sequence score using a variable order Markov model
y_cont  - Pyrimidine content between the BP adenine and the 3' splice-site
ppt_off - Polypyrimidine tract offset relative to the BP adenine
ppt_len - Polypyrimidine tract length
ppt_scr - Polypyrimidine tract score
svm_scr - Final BP score using the SVM classifier
```
    
### EXAMPLE:

The command

``` 
./svm_bpfinder.py -i test_A1CF_intron.fa -s Hsap -l 100 -d 10
```

scans the 3'-most 100nts of every sequence contained in the FASTA file named 'test_A1CF_intron.fa', 
using a human-specific model ('Hsap'), and predicts BPs as close as 10nt from the 3' splice-site.  

The output should be

```
seq_id    agez      ss_dist   bp_seq    bp_scr    y_cont    ppt_off   ppt_len   ppt_scr   svm_scr
A1CF;chr10:52623841-52645340:-     63   89   atatgattc  -1.18388285163      0.761904761905      1         11        25        0.35390249
A1CF;chr10:52623841-52645340:-     63   76   ctgtaatac  0.931808096058      0.760563380282      13        57        99        1.1116549
A1CF;chr10:52623841-52645340:-     63   70   tactaaccg  3.97557240521       0.784615384615      7         57        99        2.6909905
A1CF;chr10:52623841-52645340:-     63   31   ttttgaccc  1.77449397266       0.807692307692      1         24        44        1.7040359
A1CF;chr10:52623841-52645340:-     63   16   ctctcaccc  2.75935904557       0.727272727273      1         9         16        1.8028404
```
 
# NOTE:

SVM_BPfinder works by calling two scripts/programs: 
  1) 'SCRIPTS/svm_getfeat.py' - collects features 
  2) 'SCRIPTS/svm_classify'   - scores candidate BPs

This creates two files in the working directory, which are removed once the final results are displayed.
In case the process gets interrupted, these two files might be left on your system.  		 
	
------------------------------------------
