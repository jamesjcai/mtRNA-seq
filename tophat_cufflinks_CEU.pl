#!/usr/bin/perl
@files=qw/ERR009096 ERR009097 ERR009099 ERR009153 ERR009102 ERR009103 ERR009104 ERR009105 ERR009106 ERR009107 ERR009108 ERR009109 ERR009110 ERR009111 ERR009113 ERR009114 ERR009115 ERR009117 ERR009118 ERR009119 ERR009121 ERR009122 ERR009123 ERR009124 ERR009126 ERR009127 ERR009129 ERR009130 ERR009131 ERR009132 ERR009133 ERR009134 ERR009135 ERR009136 ERR009137 ERR009139 ERR009140 ERR009141 ERR009142 ERR009143 ERR009144 ERR009145 ERR009146 ERR009147 ERR009149 ERR009150 ERR009151 ERR009152 ERR009154 ERR009155 ERR009156 ERR009157 ERR009159 ERR009160 ERR009163 ERR009164 ERR009165 ERR009166 ERR009167 ERR009168/;
@innerdist=qw/8 19 117 10 2 16 119 156 22 49 115 1 37 -5 120 3 23 -31 108 6 137 7 120 7 16 48 143 25 -20 29 10 152 9 -9 4 52 23 125 37 113 11 147 11 7 6 5 -16 48 16 133 52 -16 -16 1 -17 5 15 -7 25 15/;

for ($count = 0; $count < scalar(@files); $count++) {
	$fileid = $files[$count];
    $indist = $innerdist[$count];
    if (-e "transcriptome_data/known.gff"){
	    	system "tophat2 -p 8 --read-mismatches 3 -g 1 -r $indist --mate-std-dev 30 --segment-length 12 --transcriptome-index=transcriptome_data/known -o out/$fileid --no-novel-juncs idx/hg19ceu fq/$fileid.fastq";
	}else{
	       	system "tophat2 -p 8 --read-mismatches 3 -g 1 -r $indist --mate-std-dev 30 --segment-length 12 -G idx/gencode.v11.annotation_CEU_mt.gtf --transcriptome-index=transcriptome_data/known -o out/$fileid --no-novel-juncs idx/hg19ceu fq/$fileid.fastq";
	}
	system "cufflinks -p 8 -G idx/gencode.v11.annotation_CEU_mt.gtf -o cuf/$fileid out/$fileid/accepted_hits.bam";
}
