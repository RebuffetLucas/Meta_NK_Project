BARCODES=/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/MetaNK_consortium/GSM5584154/03_Output/01_cellranger/GSM5584154/outs/filtered_feature_bc_matrixbarcodes.tsv
#BARCODES='/mnt/DOSI/JPGLAB/BIOINFO/Project/LysoDC_precursors/211217_VH00228_78_AAATHTGM5/05_Output/00_CellRanger/mm10/outs/filtered_feature_bc_matrix/barcodes.tsv'
OUT=/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/04_RNA_Velocity/Velocyto_FirstStep/$1
#MASK=/mnt/DOSI/JPGLAB/BIOINFO/Project/LysoDC_precursors/211217_VH00228_78_AAATHTGM5/01_Reference/mm10_rmsk.gtf
POSSORTED_GENOME=/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/MetaNK_consortium/CMVpos4_donorE/03_Output/01_cellranger/CMVpos4_donorE/outs/possorted_genome_bam.bam
ANNOTATIONS=/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/01_REFERENCE/01_GENOME/CellRanger/human/2020-A/refdata-gex-GRCh38-2020-A/genes/genes.gtf


echo BARCODES = $BARCODES
echo MASK = $MASK
echo GENOME = $POSSORTED_GENOME
echo ANNOTATIONS = $ANNOTATIONS

echo -------------------

echo OUT = $OUT

mkdir $OUT

#docker run -it --name scvelo -d -p 8883:8787 -v /mnt:/mnt -e PASSWORD=choupinou -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) #mparikhbroad/velocyto #genomicpariscentre/velocyto

#docker run -it --name scvelo -p 8883:8787 -v /mnt:/mnt -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -e PASSWORD=choupinou scvelo

docker run -itd --name scvelo -p 8883:8787 -v /mnt:/mnt -u $(id -u):$(id -g) -e PASSWORD=choupinou scvelo


gunzip $BARCODES.gz 

docker exec -i scvelo velocyto run \
	-b $BARCODES \
	-o $OUT \
    	$POSSORTED_GENOME \
	$ANNOTATIONS

#-m $MASK \
docker stop scvelo
docker rm scvelo

gzip $BARCODES



