CELL=$1
NEGATIV=$2
#ORT="/home/cip/nf/snfrbert/test/"
ORT="/media/angua/Gwendolyn/Masterarbeit/"
EBWT=${ORT}"Bowtie/hg19.ebwt/hg19"
CHROMSIZE=${ORT}"Sequenzen/hg19.chrom.sizes"
GENOM=${ORT}"cisGenome/hg19"
INPUT=${ORT}"Sequenzen/"${CELL}
BOWTIEOUT=${ORT}"Bowtie/output/v0bestm1strata/"${CELL}
ALNOUT=${ORT}"cisGenome/aln/"${CELL}
BAROUT=${ORT}"cisGenome/bar/"${CELL}
SEQ2OUT=${ORT}"cisGenome/seqV2/"${CELL}
SEQ1OUT=${ORT}"cisGenome/seqV1/"${CELL}


#if [$1 != "HK2"]; then
#	echo "Genom angeben HK2 oder UMRC3"
#	exit 1;
#	if [$2 = ""]; then
#		echo "Negativ Kontrolle angeben: 19686 für HK2, 150 für UMRC3" 
#		exit 1;
#	fi
#fi


echo "Running Bowtie"
#fuer single strand
cd ${INPUT}
pwd
for sample in S_*; do
  name=${sample#*_}
  name=${name#*_}
  name=${name%_*}
  name=${name%_*}
#  bowtie -m1 --best --strata -v 0 -t ${EBWT} ${sample} ${BOWTIEOUT}/${name}."txt"
#  echo "bowtie -m1 --best --strata -v 0 -t "${EBWT}" "${sample}" "${BOWTIEOUT}/${name}."txt"
done

#fuer paired end
for sample in P_*_1_sequence.txt; do
  samp1=${sample%_*} 
  samp1=${samp1%_*}_"1_sequence.txt" 
  samp2=${sample%_*} 
  samp2=${samp2%_*}_"2_sequence.txt" 
  name=${sample%_*}
  name=${name%_*}
  name=${name#*_}
  name=${name#*_}
  name=${name#*_}
#  bowtie -m1 --best --strata -v 0 -t ${EBWT} -1 ${samp1} -2 ${samp2} ${BOWTIEOUT}/${name}."txt"
#  echo "bowtie -m1 --best --strata -v 0 -t "${EBWT}" -1 "${samp1}" -2 "${samp2}${BOWTIEOUT}/${name}."txt"
done

#nur single strand
for sample in *_1_sequence.txt; do
  name=${sample#*_}
  name=${name#*_}
  name=${name%_*}
  name=${name%_*}
#  bowtie -m1 --best --strata -v 0 -t ${EBWT} ${sample} ${BOWTIEOUT}/${name}."txt"
#  echo ${name}
#  echo "bowtie -m1 --best --strata -v 0 -t "${EBWT}" "${sample}" "${BOWTIEOUT}/${name}."txt"
done


echo " "
echo "Converting to aln"
cd ${BOWTIEOUT}
# convert to aln
pwd
for sample in *; do
  name=${sample%.*}
 # echo ${name}
 # ~/Documents/Masterarbeit/cisgenome/bin/file_bowtie2aln -i ${sample} -o ${ALNOUT}/${name}."aln"
#  echo "~/Documents/Masterarbeit/cisgenome/bin/file_bowtie2aln -i "${sample}" -o "${ALNOUT}/${name}."aln"
done
#mv ${NEGATIV}."aln" "control.aln"


#echo " "
#echo "Peak Detection new version"
#cd ${ALNOUT}
#for sample in [0-9]*.aln; do
#  name=${sample%.*}
  #echo `pwd`/${sample}" 1" > ${SEQ2OUT}"/samplelist"_${name}."txt"
  #echo `pwd`/"control.aln 0" >> ${SEQ2OUT}"/samplelist"_${name}."txt"
#done
#cd ${SEQ2OUT}
#for sample in samplelist*.txt; do
#  name=${sample%.*}
#  name=${name#*_}
 # ~/Documents/Masterarbeit/cisgenome/bin/seqpeak -i ${sample} -d . -o "output"_${name}
  #echo "~/Documents/Masterarbeit/cisgenome/bin/seqpeak -i "${sample}" -d . -o ""output"_${name}
#done

echo " " 
echo "Sorting aln files"
cd ${ALNOUT}
# sort aln file
pwd
for sample in *.aln; do
	echo ${sample}
 # ~/Documents/Masterarbeit/cisgenome/bin/tablesorter_str ${sample}
  # echo "~/Documents/Masterarbeit/cisgenome/bin/tablesorter_str "${sample}
done


echo " " 
echo "Converting to bar"
cd ${ALNOUT}
# convert to bar
pwd
for sample in *.sort; do
  name=${sample%.*}
  name=${name%.*}
  echo ${name}
  ~/Documents/Masterarbeit/cisgenome/bin/hts_aln2barv2 -i ${sample} -o ${BAROUT}/${name}."bar"
  #echo "~/Documents/Masterarbeit/cisgenome/bin/hts_aln2barv2 -i "${sample}" -o "${BAROUT}/${name}."bar"
done


echo " " 
echo "Computing FDR"
cd ${BAROUT}
pwd
for sample in *[0-9].bar; do
  name=${sample%.*}
	echo ${sample}
 ~/Documents/Masterarbeit/cisgenome/bin/hts_windowsummaryv2_2sample -i ${sample} -n control.bar -g ${GENOM}/chrlist.txt -l ${GENOM}/chrlen.txt -w 100 -o ${SEQ1OUT}/${name}_"sum"
#  echo  "~/Documents/Masterarbeit/cisgenome/bin/hts_windowsummaryv2_2sample -i "${sample}" -n control.bar -g "${GENOM}/"chrlist.txt -l "${GENOM}/"chrlen.txt -w 100 -o "${SEQ1OUT}/${name}_"sum" 
done   

echo " " 
echo "Peak Detection old version"
cd ${SEQ1OUT}
pwd
for sample in *sum; do
  name=${sample%_*}
  p0=`head -n 1 ${sample}`
  p0=${p0##*=}
 echo ${name}
  ~/Documents/Masterarbeit/cisgenome/bin/hts_peakdetectorv2_2sample -i ${BAROUT}/${name}.bar -n ${BAROUT}/control.bar -d . -o ${name} -f ${sample}.fdr -c 0.1 -m 20 -w 100 -s 25 -p ${p0} -br 1 -brl 30 -ssf 0
#  echo "~/Documents/Masterarbeit/cisgenome/bin/hts_peakdetectorv2_2sample -i "${BAROUT}/${name}."bar -n "${BAROUT}/"control.bar -d . -o "${name}" -f "${sample}."fdr -c 0.05 -m 40 -w 100 -s 25 -p "${p0}" -br 1 -brl 30 -ssf 1 -cf 20 -cr 20"
done

echo " "
echo "Compute L = half of the median length of the peaks and shift and recompute FDR"
# 35 nur test, halber median der peaklength aus dem .cod berechnen
for sample in *.cod; do
  name=${sample%.*}
  echo ${sample}
  number_samples=`wc -l ${sample} | cut -d " " -f 1` 
  mid=$(( (${number_samples} + (${number_samples}%2))/2 ))
  L=`cut -f 6 ${sample} | sed "1d" | sort -n | sed -n "${mid}p"`
  cd ${BAROUT}
  cp control.bar control.bar.bck
  ~/Documents/Masterarbeit/cisgenome/bin/hts_alnshift2bar -i "${name}".bar -s $L
  ~/Documents/Masterarbeit/cisgenome/bin/hts_alnshift2bar -i control.bar -s $L
  ~/Documents/Masterarbeit/cisgenome/bin/hts_windowsummaryv2_2sample -i ${name}.bar -n control.bar -g ${GENOM}/chrlist.txt -l ${GENOM}/chrlen.txt -w 100 -o ${SEQ1OUT}/${name}_"sum2" -z 1
  mv control.bar.bck control.bar
#  echo "~/Documents/Masterarbeit/cisgenome/bin/hts_alnshift2bar -i "${name}".bar -s" $L
#  echo "~/Documents/Masterarbeit/cisgenome/bin/hts_alnshift2bar -i control.bar -s" $L
#  echo  "~/Documents/Masterarbeit/cisgenome/bin/hts_windowsummaryv2_2sample -i "${name}".bar -n control.bar -g "${GENOM}/"chrlist.txt -l "${GENOM}/"chrlen.txt -w 100 -o "${SEQ1OUT}/${name}_"sum2" -z 1 
  cd ${SEQ1OUT}
done

echo " " 
echo "Recompute Peaks"
cd ${SEQ1OUT}
pwd
for sample in *sum2; do
  name=${sample%_*}
  echo ${name}
  p0=`head -n 1 ${sample}`
  p0=${p0##*=}
  ~/Documents/Masterarbeit/cisgenome/bin/hts_peakdetectorv2_2sample -i ${BAROUT}/${name}.bar -n ${BAROUT}/control.bar -d . -o "output"_${name} -f ${sample}.fdr -c 0.1 -m 20 -w 100 -s 25 -p ${p0} -br 1 -brl 30 -ssf 0 -z 1
 # echo "~/Documents/Masterarbeit/cisgenome/bin/hts_peakdetectorv2_2sample -i "${BAROUT}/${name}."bar -n "${BAROUT}/"control.bar -d . -o output"_${name}" -f "${sample}."fdr -c 0.05 -m 40 -w 100 -s 25 -p "${p0}" -br 1 -brl 30 -ssf 1 -cf 20 -cr 20 -z1"
done

echo " "
echo "Annotate peaklist with closest genes"
cd ${SEQ1OUT}
for sample in output*.cod; do
  name=${sample%.*}
  echo ${name}
#  echo "~/Documents/Masterarbeit/cisgenome/bin/refgene_getnearestgene -d "${GENOM}"/annotation/refFlat_sorted.txt -dt 1 -s human -i "${sample}" -o" ${name}"_annotated.cod -r 0 -up 8000 -down 2000"
#  ~/Documents/Masterarbeit/cisgenome/bin/refgene_getnearestgene -d ${GENOM}/annotation/refFlat_sorted.txt -dt 1 -s human -i ${sample} -o ${name}_annotated.cod -r 0 -up 100000 -down 100000
done

#cd ${SEQ2OUT}	
#for sample in *.cod; do
#  name=${sample%.*}
#  echo "~/Documents/Masterarbeit/cisgenome/bin/refgene_getnearestgene -d "${GENOM}"/annotation/refFlat_sorted.txt -dt 1 -s human -i "${sample}" -o" ${name}"_annotated.cod -r 0 -up 8000 -down 2000"
#  ~/Documents/Masterarbeit/cisgenome/bin/refgene_getnearestgene -d ${GENOM}/annotation/refFlat_sorted.txt -dt 1 -s human -i ${sample} -o ${name}_annotated.cod -r 0 -up 50000 -down 50000
#done

echo " "
echo "Convert to wig"
cd ${SEQ1OUT}
for sample in output_*[0-9].*.bar; do
	name=${sample%.*}
	echo ${name}
#	~/Documents/Masterarbeit/cisgenome/bin/affy_bar2wig2 -i ${sample} -o ${name}.wig
#	echo "~/Documents/Masterarbeit/cisgenome/bin/affy_bar2wig2 -i ${sample} -o ${name}.wig"
done

#cd ${SEQ2OUT}	
#for sample in *log2fc.bar; do
#	name=${sample%.*}
#	~/Documents/Masterarbeit/cisgenome/bin/affy_bar2wig2 -i ${sample} -o ${name}.wig
#	echo "~/Documents/Masterarbeit/cisgenome/bin/affy_bar2wig2 -i ${sample} -o ${name}.wig"
#done

#echo " "
#echo "Convert to BigWig"
#cd ${SEQ1OUT}
#for sample in *.wig; do
#	name=${sample%.*}
	#~/Documents/Masterarbeit/cisgenome/bin/wigToBigWig ${sample} ${CHROMSIZE} ${name}.BW 
#	echo "~/Documents/Masterarbeit/cisgenome/bin/wigToBigWig ${sample} ${CHROMSIZE} ${name}.BW"
#done

#cd ${SEQ2OUT}
#for sample in *.wig; do
#	name=${sample%.*}
	#~/Documents/Masterarbeit/cisgenome/bin/wigToBigWig ${sample} ${CHROMSIZE} ${name}.BW 
#	echo "~/Documents/Masterarbeit/cisgenome/bin/wigToBigWig ${sample} ${CHROMSIZE} ${name}.BW"
#done
