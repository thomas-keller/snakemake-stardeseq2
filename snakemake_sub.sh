#! /bin/bash

#SBATCH --job-name=snakemake_gen
#SBATCH --time=120:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=8G
#SBATCH -o snakestar.e.%j.%N.txt
#SBATCH -e snakestar.o.%j.%N.txt

# this file is snakemake_sub.sh
#submit with
# sbatch --cpus-per-task=2 --mem=8g snakemake_sub.sh
#taken from biowulf example
#https://hpc.nih.gov/apps/snakemake.html

sbcmd="sbatch --cpus-per-task={threads}  --mem-per-cpu={cluster.mem-per-cpu-mb}"
sbcmd+=" --nodes={cluster.nodes} --ntasks={cluster.ntasks}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={cluster.output}"

#snakemake --unlock
#snakemake -pr --keep-going \
#             --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
#             --latency-wait 120 all1S

mkdir /work/t/tekeller/deseq2

cp ~/snakemake-stardeseq2/snakemake_sub.sh /work/t/tekeller/deseq2_ex
cp ~/snakemake-stardeseq2/cluster.json /work/t/tekeller/deseq2_ex
cp ~/snakemake-stardeseq2/config.yml /work/t/tekeller/deseq2_ex
cp ~/snakemake-stardeseq2/Snakefile /work/t/tekeller/deseq2_ex
cp ~/snakemake-stardeseq2/deseq_star.R /work/t/tekeller/deseq2_ex
cp ~/snakemake-stardeseq2/report/* /work/t/tekeller/deseq2_ex
cp ~/snakemake-stardeseq2/deseq.Rmd /work/t/tekeller/deseq2_ex
#cp ~/snakemake-stardeseq2/samples.json /work/t/tekeller/deseq2_ex

cd /work/t/tekeller/deseq2_ex

source activate snakemake
snakemake --unlock
snakemake -R deseq -j20 -pr --cluster "sbatch --time={cluster.time} --partition={cluster.partition} --mem={cluster.mem} --cpus-per-task={cluster.cpus} --out={cluster.output}" --cluster-config cluster.json --keep-going --latency-wait 900 --rerun-incomplete

#snakemake --report report.html
