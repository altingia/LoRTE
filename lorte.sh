#PBS -S /bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=28:cmbnode -l mem=60gb
#PBS -V
#PBS -N LoRTEv1.2_sapelo
#PBS -o /lustre1/sh60271/LoRTE/0414_test2
#PBS -e /lustre1/sh60271/LoRTE/0414_test2
#PBS -M shhan@uga.edu
#PBS -m ae
#PBS -j oe

cd $PBS_O_WORKDIR

output_dir="/lustre1/sh60271/LoRTE/0414_test2"
mkdir -p $output_dir
rm -rf $output_dir/*

module load python/2.7.8
module load ncbiblast+/2.6.0

input_dir="/lustre1/sh60271/LoRTE/0411_test/Input_files"
python LoRTEv1_2_pa.py -r $input_dir/Genome_Dmel5.fa -L $input_dir/List_TE_Dmel5 -i $input_dir/Test_Pacbio_Dmel_1X_Coverage.fa -c $input_dir/Consensus_TE_Droso.fa -o $output_dir/output -n 27