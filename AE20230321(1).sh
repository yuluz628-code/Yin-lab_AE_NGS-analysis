export PATH=/data/junyiDuan/miniconda3/envs/AMP/bin:$PATH
cat reference*.txt | tr '-' '_'|tr '[a-z]' '[A-Z]' > REFERENCE_STANDARD_FORM.txt
fastp -w 6 -i *R1*.f*q.gz -I *R2*.f*q.gz -c -m --merged_out merg_rawR1.fq  --overlap_len_require 12
trimmomatic SE merg_rawR1.fq clean_merg_rawR1.fq \
	ILLUMINACLIP:/home/junyiDuan/miniconda3/envs/AMP/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true AVGQUAL:30
mkdir -p samples
for barcode in `tail -n+2 REFERENCE_STANDARD_FORM.txt | cut -d '	' -f 1`
do
sample_ID=`grep -w "$barcode" REFERENCE_STANDARD_FORM.txt |awk '{print $2}'`
BARCODE1=`echo ${barcode} |awk '{print $1}' | cut -d '+' -f 1`
revBARCODE2=`echo ${barcode}+ |awk '{print $1}' | cut -d '+' -f 2 |tr '[ATGC]' '[TACG]'|rev`
echo "$sample_ID barcode -> barcode1 ${BARCODE1} barcode2 $revBARCODE2"
grep  -A 2 -B 1 -P "^${BARCODE1}.+${revBARCODE2}$" clean_merg_rawR1.fq |grep -v -P '^--$' > samples/${barcode}.fq
mv samples/${barcode}.fq samples/${sample_ID}.fq
done

mkdir -p clean
for trim in ./samples/*fq
do
	no_path_sn=${trim##*/};sn=${no_path_sn%.*};
	temp_minlength=`grep -w "$sn" REFERENCE_STANDARD_FORM.txt |awk '{print $5}'`
	echo "sample $sn is being filtered, length <= $temp_minlength will be filtered"
	fastp -i  ${trim} -f 6 -t 5  -o ./clean/${trim##*/}  --length_required $temp_minlength
done

for sample_name in ./clean/*.fq
do
  no_path_sn=${sample_name##*/};sn=${no_path_sn%.*};notrim_sn=${sn#trim}

temp_reference_sequence=`grep -w "${notrim_sn}" REFERENCE_STANDARD_FORM.txt|awk '{print $3}'| sed 's/_//g'`
temp_reference_length=`echo "$temp_reference_sequence"|wc -c`
 echo  "$no_path_sn is being aligned, wt reference is $temp_reference_sequence"
 ((plot_table_length=${temp_reference_length}/2-3))
quantification_size=`echo "5"`
((a=$temp_reference_length/2-10))
((b=a+20))
mid_gRNA=`echo "$temp_reference_sequence"|cut -b ${a}-${b}`
echo "$no_path_sn is being aligned, wt reference is $temp_reference_sequence assumpted gRNA is$mid_gRNA "
CRISPResso -r1 ${sample_name} -an ${notrim_sn} -a  ${temp_reference_sequence}  -g ${mid_gRNA} -q 30  \
	--exclude_bp_from_left 0 --exclude_bp_from_right 0 --min_frequency_alleles_around_cut_to_plot 0.01 \
	-wc -10 -w $quantification_size --plot_window_size $plot_table_length  --min_identity_score 8 --ignore_substitutions

echo "******************************************************************************************************"
done

mkdir -p reports
mv CRI*/*Alleles_frequency_table_around_sgRNA*.txt ./
mv ./CRISPResso_on_*/*Alleles_frequency_table_around_sgRNA_*.pdf ./reports/
rm -r fastp* clean samples *.fq CR*

echo "***************************************************************************************************************"
echo "requantification start"
echo "***************************************************************************************************************"
for input_allele_frequency_table in ./*.Alleles_frequency_table_around_sgRNA*.txt
do
nopathsamplename=${input_allele_frequency_table##*/};samplename=${nopathsamplename%%.Alleles_frequency_table*}
MIN_FREQ=`grep -w "$samplename" REFERENCE_STANDARD_FORM.txt|awk '{print $6}'`
echo "min frq is $MIN_FREQ"
cat $input_allele_frequency_table |awk -vmin_freq="$MIN_FREQ" '$8>min_freq{print $0}'|awk '$NF="";{print $0}' > ${samplename}summary.txt
done
for summary in *summary.txt
do
        total_read_sum=`cat $summary |awk '{sum += $7};END {print sum}'`
echo "${summary}_report">> summary_of_all_sample.txt
echo "amplicon_sequence	reference_sequence	nill1	null2	null3	null4	reads_number	frequencey">> summary_of_all_sample.txt
cat $summary |while read line
        do
        per_line_sum=`echo $line |awk '{print $7}'`
        frequency=`echo "scale=4;${per_line_sum}/${total_read_sum}*100"| bc |awk '{printf "%.2f", $0}'`
        echo "$line $frequency" >> summary_of_all_sample.txt
        done
done
#mkdir -p reports
sed 's/ /\t/g' summary_of_all_sample.txt | sed 's/-//g'| cut -f 1,2,7,8 >summary_of_all_sample_tab.txt
rm *temp_allele_frequency_table.txt summary_of_all_sample.txt
rm *summary.txt standard_ref.txt
mv summary_of_all_sample_tab.txt ./reports
