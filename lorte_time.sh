#PBS -S /bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=28:cmbnode -l mem=60gb
#PBS -V
#PBS -N LoRTEv1.2_time
#PBS -o /lustre1/sh60271/LoRTE/0414_test
#PBS -e /lustre1/sh60271/LoRTE/0414_test
#PBS -M shhan@uga.edu
#PBS -m ae
#PBS -j oe

cd $PBS_O_WORKDIR

module load python/2.7.8
module load ncbiblast+/2.6.0

dir="/lustre1/sh60271/LoRTE/0414_test"
mkdir -p $dir
rm -rf $dir/*

output_dir="$dir/0413_27_cores"
mkdir -p $output_dir
rm -rf $output_dir/*

input_dir="/lustre1/sh60271/LoRTE/0411_test/Input_files"
python LoRTEv1_2_pa_time.py -r $input_dir/Genome_Dmel5.fa -L $input_dir/List_TE_Dmel5 -i $input_dir/Test_Pacbio_Dmel_1X_Coverage.fa -c $input_dir/Consensus_TE_Droso.fa -o $output_dir/output -n 27


# output_dir="$dir/0413_2_cores"
# mkdir -p $output_dir
# rm -rf $output_dir/*

# input_dir="/lustre1/sh60271/LoRTE/0411_test/Input_files"
# python LoRTEv1_2_pa_time.py -r $input_dir/Genome_Dmel5.fa -L $input_dir/List_TE_Dmel5 -i $input_dir/Test_Pacbio_Dmel_1X_Coverage.fa -c $input_dir/Consensus_TE_Droso.fa -o $output_dir/output -n 2