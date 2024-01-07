# Dapars

echo Annotated_3UTR=/home/yangyb/yangyb/dapars_analysis_file/hg19_refseq_extracted_3UTR_2019.bed >> configure_file.txt

echo "Aligned_Wig_files=`ls *bedgraph|xargs -n1|tr "\n" ","|sed 's/,$//'`

Output_directory=Dapars_data/

Output_result_file=DaPars_data

Coverage_threshold=30

Num_Threads = 30

sequencing_depth_file=library.txt

" >> configure_file.txt

seq 1 22|xargs -n1 -I {} -P 22 python Dapars2_Multi_Sample.py configure_file.txt  chr{}

