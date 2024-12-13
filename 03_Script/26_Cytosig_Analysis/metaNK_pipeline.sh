### 1. Cytosig run to correlate different NK subtypes (new version of Cytosig to include more cytokines in beta phase)
conda activate py39

DATA=/Users/dan.sun/Documents/Research/Niklas_projects/5.metaNK/data
cd $DATA
for i in {1..11}
do
echo $i
CytoSig_run.py -i $DATA/seu_cytosig_mat_part$i.txt -o cytosig_part$i -s 2
done


# 2023/06/28 new clustering from Lucas
DATA=/Users/dan.sun/Documents/Research/Niklas_projects/5.metaNK/data
cd $DATA
for i in {1..4}
do
echo $i
CytoSig_run.py -i $DATA/new_seu_cytosig_mat_part$i.txt -o new_cytosig_part$i -s 2
done
