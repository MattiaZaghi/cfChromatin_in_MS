#for a given bam file test overlap with on- and off-target regions
#calculate enrichment for on-target reads

#note: currently works only for paired end reads - will need ot modify for single-end

bam=$1
targets=$2 #on.target.filt.bed
name=$3

printf 'name\tmark\ton_bp\toff_bp\ton_reads\toff_reads\tenrichment\n'
for on in ${targets}/*/on.target.filt.bed; do
off=`echo $on | sed 's/on.target/off.target/'`
mark=`echo $on | sed 's/\/on.target.filt.bed//' | sed 's/.*\///'`
	onbp=`awk 'BEGIN{sum=0}{sum = sum + $3-$2}END{print sum}' $on`
	offbp=`awk 'BEGIN{sum=0}{sum = sum + $3-$2}END{print sum}' $off` 
	onreads=`samtools view -c $bam -L $on`
	offreads=`samtools view -c $bam -L $off`                                              
	enr=`echo $onreads $onbp $offreads $offbp | awk '{print ($1/$2) / ($3/$4)}'`                
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $name $mark $onbp $offbp $onreads $offreads $enr
done
