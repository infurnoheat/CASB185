#!/bin/bash
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate myconda

# define input arguments and variables
all_ids=$1
acc_id=$(head -n ${SGE_TASK_ID} ${all_ids} | tail -n 1)
av_len=$2
sd_len=$3

# run quant if not ran already
if ! [ -d "kallisto_data_2/${acc_id}" ]; then
	kallisto quant -i indeces/homo_sapiens/transcriptome.idx -t 8 -o kallisto_data_2/${acc_id} -b 100 --single -l $av_len -s $sd_len ${acc_id}/*gz
fi 
