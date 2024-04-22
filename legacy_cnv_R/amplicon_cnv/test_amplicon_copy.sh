# supply snappy primer counts file  list in script directory as only input - eg. to run ./test_amplicon_copy.sh list_of_primer_counts_files
## requires normal counts file in script directory. Example supplied. CAn be built from averaged normal (flat cnv) sample primer outputs. Needs proportion column ( individual  primer read proportion in respect to total primer reads for all primers) 
## output of primer normalization passed to Rscript plot_amplicon.r  here.
counts_file_list=$1

for line in $(cat $counts_file_list)
do
    echo $line
    base=$(basename  $line .counts.PRIMER)
    echo $base
    tail -n+8 $line > primerout1$base
    sum1=$(awk 'BEGIN{FS="\t"; OFS="\t"}; {sum += $3} END {print sum}' primerout1$base)
    echo $sum1
    sleep 4s
    awk -v var2="$sum1" 'BEGIN{FS="\t"; OFS="\t"} {print $0, $3=$3/var2, var2}' primerout1$base |  awk -v var="$base" 'BEGIN{FS="\t"; OFS="\t"} {print $0, var}' - | sort  -k1 - | tail -n +2 > outprimermetrics$base
    paste outprimermetrics$base normal_relative_copy.txt | awk 'BEGIN{FS="\t"; OFS="\t"} {print $0, $7=$7/$10}' > outprimermetrics${base}_normalized
    awk 'BEGIN{FS="\t"; OFS="\t"} {print $0, log($14)/log(2)}' outprimermetrics${base}_normalized > plotting_data_${base}
    Rscript plot_amplicon.r plotting_data_${base} ${base}.pdf
    rm outprimermetrics${base}*  primerout1$base
done
