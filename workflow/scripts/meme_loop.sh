ups_dir=$1
meme_mode=$2
motif_length=$3
tf_dir=$4


for fasta in `ls -d $ups_dir/*`; do 
  echo $fasta;
  cog=$(echo $fasta | rev | cut -d"/" -f1 | rev | cut -d'_' -f1,2);
  length=$(echo $fasta | rev | cut -d"/" -f1 | rev | cut -d'_' -f4 | cut -d'.' -f1);
  meme $fasta -dna -mod ${meme_mode} -nmotifs 4 -maxw ${motif_length} -oc ${tf_dir}/meme/meme_${cog}_${length}_${meme_mode}_${motif_length} -revcomp;
done