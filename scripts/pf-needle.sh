#!/bin/bash
###You need SRR** ID list.
###You need the optional 'link_ena' of trouble samples, which is obtained from ENA archieve.
###Finding a bug that conda env activation will lead prefetch's version degrading, which I have encountered this before.
###Run in base env.
###Analysis dataset of mef2c with needle, which is dealed with prefetch and fasterq-dump in NGS env.
###Add "||" for the bug in salmon owning to the trimming in 20231126.

###If you run in a new env, alter the prefetch.2.11.0 to right name. 
###Note you need to set work directory, salmon index directory, gtf file path, prefix of library SRR IDs and prefix of species SRR IDs.  

###Just for 112 samples of wildtype  and tac mouse  in 20240903.
###Change codes for check 'quant.genes.sf' in 20240626.
### "if [[ $line != SRR587* ]]; then " not work for SRR1852* in 20240904.(#Wrong with '-r $WD/fastq/$line".fastq"' before 20240904!!!!!!!!!!!!!)

echo ==========================================$(date)==========================================
echo $PWD
echo ==========================================$(date)==========================================


##############################
###01 Preparation
##############################
#source activate && conda deactivate && conda activate NGS
source activate && conda activate NGS
prefetch.2.11.0 --version
fasterq-dump --version
ascp --version
fastqc --version
stringtie --version
hisat2 --version
samtools --version
rabbit_qc -v
salmon --version
needle --version


###############################################################################################
#You must set followings!
###############################################################################################
WD='/home/user_li/linqin_tmp/Analysis/RNA-seq/tac112/'
indexdir_hsa='/home/user_li/linqin_tmp/Genome/hsa/release107/index_human_r107_salmon_k31'
gtfpath_hsa='/home/user_li/linqin_tmp/Genome/hsa/release107/Homo_sapiens.GRCh38.107.gtf'

indexdir_mus='/home/user_li/linqin_tmp/Genome/mus/release107/index_mouse_r107_salmon_k31'
gtfpath_mus='/home/user_li/linqin_tmp/Genome/mus/release107/Mus_musculus.GRCm39.107.gtf'

indexdir_rat='/home/user_li/linqin_tmp/Genome/rat/release107/index_rat_r107_salmon_k31'
gtfpath_rat='/home/user_li/linqin_tmp/Genome/rat/release107/Rattus_norvegicus.mRatBN7.2.107.gtf'

###############################################################################################





cd $WD

if [[ -d $WD/fastq ]]; then
	echo fastq exit!
else
	mkdir fastq
fi

if [[ -d $WD/salmon_out ]]; then
	echo salmon_out exit!
else
	mkdir salmon_out
fi

if [[ -d $WD/salmon_out/quant_out ]]; then
	echo quant_out exit!
else
	mkdir salmon_out/quant_out
fi

if [[ -d $WD/log_rabbitQC ]]; then
	echo log_rabbitQC exit!
else
	mkdir log_rabbitQC
fi

if [[ -f $WD/md5_run ]]; then
	echo md5_run exit!
else
	echo $WD/md5_run NOT exit!
fi

#add 20240626
if [[ -f $WD/link_ena ]]; then
	echo link_ena exit!
else
	echo $WD/link_ena NOT exit! $$ touch $WD/link_ena
fi



##############################
###02 Download&Analysis
##############################

cd $WD

for line in $(cat $WD/runList)
do

	#if [[ -f $WD/log_rabbitQC/$line".html" && -f $WD/salmon_out/quant_out/$line/quant.genes.sf ]]; then
	if [[ -f $WD/salmon_out/quant_out/$line/quant.genes.sf ]]; then
		echo $line has DONE!
	else

		cd $WD/fastq
		echo =============================================
		echo $line
		echo =============================================
		# Ensure integraty of downloading sra files
		if [[ $( grep $line $WD/link_ena | wc -l ) -gt 0 ]]; then
			grep $line $WD/link_ena > $WD/link_tmp
			wget -nc -q --progress=dot:binary -P $WD/fastq/ -i $WD/link_tmp && echo wget $line OK! || echo wget $line failure!
			gunzip $WD/fastq/$line*".gz"  #wget download with .gz, gunzip will delete the original files
		else
			offset=1
			until [[ $offset == "" ]]
			do
				prefetch.2.11.0 --max-size 100G $line -O $WD/fastq/
				#Link is the information of $line accession
				link='https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc='$line
		
				#curl $link | grep "md5" | sed 's/.* "//g' | sed 's/".*//g' > $WD/md5_sra
				curl $link | grep "md5" | sed -n 's/.*"md5": "\([^"]*\)".*/\1/p' > $WD/md5_sra #change owning to the NCBI json format in 20240210
				md5sum $line".sra" | sed "s/ .*//g" > $WD/md5_tmp
				offset=$( diff $WD/md5_sra $WD/md5_tmp ) &&  echo diff OK!
	
			done
			cat $WD/md5_tmp >> $WD/md5_run
			#Add --threads 24 in 20231119
			fasterq-dump --threads 24 --split-files $WD/fastq/$line".sra" -O $WD/fastq && echo fasterq-dump $line OK! 
		fi

		#-i $WD/fastq/$line"_1.fastq.gz"
		## Trim reads by rabbitQC, and GSE135055 read length is about 125nt(not insert length); so I set --length_required 100'(default); (for salmon index k=31 when read length > 75!)
################
################
		###############################################################################################
		#You must set followings!
		###############################################################################################
		#Check the library, paired or single
		#$line != SRR587*
		if [[ $line == SRR1177* || $line == SRR1852* || $line == SRR5276* || $line == SRR6868* || $line == SRR9331* || $line == SRR9667* ]]; then
			rabbit_qc -c --detect_adapter_for_pe -w 24 -h $line".html" \
				-i $WD/fastq/$line"_1.fastq" \
				-I $WD/fastq/$line"_2.fastq" \
				-o $line".read1.fq" -O $line".read2.fq"
			date
			echo rabbitQC $line OK!

			###############################################################
			#cd $WD/needle_out/tmp/
			#needle minimiser $WD/fastq/$line*.fq -t 36 --paired

			#needle ibfmin $line*.minimiser -e 10 -t 36 -f 0.3 -c -o $line #Note two pair-end data need to set two "-e" 

			#needle estimate $WD/mef2c_hm.fasta -i $line -o $WD/needle_out/$line"_counts.txt" #must with ".fasta" suffix

			#echo ------------------------------------------$line needle OK!------------------------------------------
			#echo ------------------------------------------$line needle OK!------------------------------------------
			###############################################################

			cd $WD/salmon_out

			###############################################################################################
			#You must set followings!
			###############################################################################################
			#Check the species, just human or mouse here
			#str1="SRR831" #some data are from mouse
			#if [[ $line =~ $str1 ]]; then
			if [[ $line == kong1* ]]; then
				indexdir=$indexdir_hsa
				gtfpath=$gtfpath_hsa
			elif [[ $line == kong2* ]]; then
				indexdir=$indexdir_rat
				gtfpath=$gtfpath_rat	
			else
				indexdir=$indexdir_mus
				gtfpath=$gtfpath_mus
			fi
	
			salmon quant \
				-i $indexdir \
				-l A \
				-1 $WD/fastq/$line".read1.fq" \
				-2 $WD/fastq/$line".read2.fq" \
				-p 24 \
				-o $WD/salmon_out/quant_out/$line \
				-g $gtfpath \
				--numGibbsSamples 50 \
				--auxDir aux_info \
				--seqBias --gcBias -d --posBias --hardFilter --discardOrphansQuasi --writeUnmappedNames && echo salmon $line OK! 
			#changed in 20240626
			if [[ -f $WD/salmon_out/quant_out/$line/quant.genes.sf ]]; then
				echo $line quant with the fastq trimmed by rabbit_qc!
			else
				#Wrong with '-r $WD/fastq/$line".fastq"' before 20240904!!!!!!!!!!!!! 
				salmon quant \
					-i $indexdir \
					-l A \
					-1 $WD/fastq/$line"_1.fastq" \
					-2 $WD/fastq/$line"_2.fastq" \
					-p 24 \
					-o $WD/salmon_out/quant_out/$line \
					-g $gtfpath \
					--numGibbsSamples 50 \
					--auxDir aux_info \
					--seqBias --gcBias -d --posBias --hardFilter --discardOrphansQuasi --writeUnmappedNames && echo salmon using raw reads $line OK!
			fi

			# note: rabbitQC(ktrim) outputs with suffix, "read1.fa" and "read2.fq"; change "--dumpEq" to "-d";add "--hardFilter" and "--discardOrphansQuasi"; "-l ISR"  to "-l A"

################
################
		else
			rabbit_qc -c --detect_adapter_for_pe -w 24 -h $line".html" \
				-i $WD/fastq/$line".fastq" \
				-o $line".fq"
			date
			echo rabbitQC $line OK!

			cd $WD/salmon_out

			###############################################################################################
			#You must set followings!
			###############################################################################################
			#Check the species, just human or mouse here
			#str1="SRR831" #some data are from mouse
			#if [[ $line =~ $str1 ]]; then
			if [[ $line == kong1* ]]; then
				indexdir=$indexdir_hsa
				gtfpath=$gtfpath_hsa
			elif [[ $line == kong2* ]]; then
				indexdir=$indexdir_rat
				gtfpath=$gtfpath_rat	
			else
				indexdir=$indexdir_mus
				gtfpath=$gtfpath_mus
			fi
	
			salmon quant \
				-i $indexdir \
				-l A \
				-r $WD/fastq/$line".fq" \
				-p 24 \
				-o $WD/salmon_out/quant_out/$line \
				-g $gtfpath \
				--numGibbsSamples 50 \
				--auxDir aux_info \
				--seqBias --gcBias -d --posBias --hardFilter --discardOrphansQuasi --writeUnmappedNames && echo salmon $line OK!
			#changed in 20240626
			if [[ -f $WD/salmon_out/quant_out/$line/quant.genes.sf ]]; then
				echo $line quant with the fastq trimmed by rabbit_qc!
			else
				salmon quant \
					-i $indexdir \
					-l A \
					-r $WD/fastq/$line".fastq" \
					-p 24 \
					-o $WD/salmon_out/quant_out/$line \
					-g $gtfpath \
					--numGibbsSamples 50 \
					--auxDir aux_info \
					--seqBias --gcBias -d --posBias --hardFilter --discardOrphansQuasi --writeUnmappedNames && echo salmon using raw reads $line OK!
			fi

			# note: rabbitQC(ktrim) outputs with suffix, "read1.fa" and "read2.fq"; change "--dumpEq" to "-d";add "--hardFilter" and "--discardOrphansQuasi"; "-l ISR"  to "-l A"
		fi
################
################
		mv $WD/fastq/$line*"html" $WD/log_rabbitQC/
        	rm $WD/fastq/$line*

		echo mv-rm $line Done!
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		echo $line DONE!
		echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	fi
done


echo =================================================================================
echo $(date)
echo ALL DONE!
echo $(date)
echo =================================================================================
