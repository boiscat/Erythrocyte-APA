#mkdir pdui_res

head -n1 Dapars_data_chr1/DaPars_data_result_temp.chr1.txt >> pdui_res/pdui_result.txt

cat Dapars_data_chr*/tmp/Each_processor_3UTR_Result_1.txt >> pdui_res/pdui_result.txt

#done
