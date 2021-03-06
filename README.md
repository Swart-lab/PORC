# PORC
![Prediction Of Reassigned Codons](example/bsto_atcc.weblogo.png)

This code predicts genetic codes with an alternative normalization scheme to that implemented in the program [FACIL](http://facil.cmbi.umcn.nl/facil/cgi-bin/display.pl?disp=home). 

PORC was used to predict a few new genetic codes in cilates using [MMETSP transcriptomes](https://pubmed.ncbi.nlm.nih.gov/27426948/), recorded as genetic codes 27-29 in [NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi), and also in analyses of [karyorelict genetic codes](https://doi.org/10.1101/2022.04.12.488043).

## Installation

`git clone https://github.com/Swart-lab/PORC/`

`cd PORC`

## Running

PORC can be run with the following command:

`./PORC.sh $PATH_TO_PORC_SCRIPTS HMMER_PROCESSES GENOME_FILE_PREFIX GENOME_FILE_SUFFIX PFAM_PATH`

See comments in `PORC.sh` for more details on how to run. 

