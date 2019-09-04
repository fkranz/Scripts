cd ~/Data/genomes/EPAS1/snps
python3 ~/Forschung/Scripts/edit-genome.py EPAS1-hg38.fa prime-wo-dups.csv
for i in rs*.fa; do bowtie-build EPAS1-hg38.fa,$i ebwt/${i%.fa}; done
for i in rs*.fa; do mkdir ~/Forschung/Med4/EPAS1-snps/${i%.fa}; done
cd ~/Forschung/Med4/rna-reads
for j in ~/Data/genomes/EPAS1/snps/rs*.fa; do
    temp=${j##*/}
    snp=${temp%*.fa}
    echo $snp
    for i in op120*; do
        bowtie -v 0 -m 1 -p8 ~/Data/genomes/EPAS1/snps/ebwt/$snp $i ~/Forschung/Med4/EPAS1-snps/$snp/${i%.*}.bow
    done
    cd ~/Forschung/Med4/EPAS1-snps/$snp/
    echo $snp >> ~/Forschung/Med4/EPAS1-snps/results
    for i in *; do
        echo $i >> ~/Forschung/Med4/EPAS1-snps/results
        numberEP=`grep EP $i | wc -l`
        numberRS=`grep rs $i | wc -l` 
        echo "SNP vs Normal:" $numberRS "/" $numberEP >> ~/Forschung/Med4/EPAS1-snps/results
        echo "" >> ~/Forschung/Med4/EPAS1-snps/results
        if [ $numberEP -ne 0 ]; then
            if [ $numberRS -gt 3 ]; then
                echo $snp >> ~/Forschung/Med4/EPAS1-snps/interesting-results
                echo $i >> ~/Forschung/Med4/EPAS1-snps/interesting-results
                echo "SNP vs Normal:" $numberRS "/" $numberEP >> ~/Forschung/Med4/EPAS1-snps/interesting-results
                echo "" >> ~/Forschung/Med4/EPAS1-snps/interesting-results
            fi
        fi
    done
    echo "" >> ~/Forschung/Med4/EPAS1-snps/results
    cd ~/Forschung/Med4/rna-reads
done