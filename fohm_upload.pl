#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename qw(basename dirname);
use POSIX qw( strftime);

#sftp -i .ssh/id_rsa_gensam se120@gensam-sftp.folkhalsomyndigheten.se

# ARGUMENTS: 1: pipeline output folder, 2: outputfolder for this script (files to be sent to fohm 3. metadatafile
my $in_dir = $ARGV[0];
##my $gisaid_log = $ARGV[1];
my $out_dir  = $ARGV[1];
my $metadata_csv = $ARGV[2]; 

my $SCRIPT_ROOT = dirname($0);

##my %config = read_config($SCRIPT_ROOT.'/fohm.config');
my %Locations = read_config($SCRIPT_ROOT.'/LocationsClean.config');
my %Adresses = read_config($SCRIPT_ROOT.'/AdressesClean.config');


##my %id_conversion_table = read_conversion_table($SCRIPT_ROOT.'/conversion.table');
##my %gisaid_ids = read_gisaid_ids($gisaid_log);
my %metadata = read_pat_metadata($metadata_csv);

##print $metadata{"midnight_barcode04"}{InternalLabId} . "\n";
##print $metadata{"midnight_barcode04"}{labcode} . "\n";
##print $metadata{"midnight_barcode04"}{LibID} . "\n";

##my $full_out_dir = $config{out_dir}.'/'.basename($in_dir);
##my $prefix = $full_out_dir.'/'.$config{region_code}."_".$config{lab_code};
##system("mkdir -p $full_out_dir");

##my @vcfs = glob "$in_dir/*.freebayes.vep.vcf";

# Get pangolin file
##my $pangolin_fn = "$in_dir/pangolin_all.csv";
my $pangolin_fn = "$in_dir/analysisReport.tsv";
die "No pangolin output!" unless -e $pangolin_fn;
my %Pangolin_data = read_PipelineResults($pangolin_fn);

##print $Pangolin_data{"midnight_barcode04"}{lineage} . "\n";

my $date = strftime '%Y-%m-%d', localtime;
   



opendir my $dir, $in_dir or die "Cannot open directory: $!"; 
my @files = readdir $dir;
closedir $dir;



## fix fohm



my $GlobRegion ="";
my $GlobLab = ""; 
foreach my $SampleID (keys % metadata){
    ##print $SampleID . ";" . $metadata{"midnight_barcode04"}{LibID} .  ";"  .  $metadata{$SampleID}{labcode}   .  "\n";
    my $fohm_prefix =  $out_dir . "/" .  $metadata{$SampleID}{LibID} .  "_"  . $metadata{$SampleID}{regioncode} . "_"  .  $metadata{$SampleID}{labcode}  ;
    my @filesR1 = glob( $in_dir . '/' . $SampleID . "*" . "R1"  ."*" . "fastq.gz" );
    if(0+@filesR1>1){die "multiple files match ID $SampleID"};
    if(0+@filesR1<1){die "no fastq-files assocated with $SampleID"};
    my $fq1 = $filesR1[0];
    my @filesR2 = glob( $in_dir . '/' . $SampleID . "*" . "R2"  ."*" . "fastq.gz" );
    if(0+@filesR2>1){die "multiple files match ID $SampleID"};
    if(0+@filesR2<1){die "no fastq-files assocated with $SampleID"};
    my $fq2 = $filesR2[0];
    print $fq1 . "\n";
    print $fq2 . "\n";
     system "ln -sf $fq1 ${fohm_prefix}_1.fastq.gz";
     system "ln -sf $fq2 ${fohm_prefix}_2.fastq.gz";    
    $GlobRegion = $metadata{$SampleID}{regioncode};
    $GlobLab = $metadata{$SampleID}{labcode};


} 

## gisaid - columns out of sync!!! FIX!!  DO NOT USE GISAID

open(GISAID, ">". $out_dir .'/'. $GlobRegion  ."_". $GlobLab  ."_".$date."_GISAIDsubmission.csv");
print GISAID "submitter,fn,covv_virus_name,covv_type,covv_passage,covv_collection_date,covv_location,covv_add_location,covv_host,covv_gender,covv_patient_age,covv_patient_status,covv_seq_technology,covv_orig_lab,covv_orig_lab_addr,covv_subm_lab,covv_subm_lab_addr,covv_subm_sample_id,covv_authors\n";

foreach my $SampleID (keys % metadata){
print GISAID join(",", ( $metadata{$SampleID}{Submitter},
                        ## basename($fasta_fn),
			 "fastaseq",
		         "",
                         $metadata{$SampleID}{Type},
			 
			 $metadata{$SampleID}{'Passage details/history'},##Passage details/history
			 $metadata{$SampleID}{'Collection date'},
			 
			 $Locations{ $metadata{$SampleID}{labcode} }, ## convert according to table  $Locations{}
			 $Adresses{ $metadata{$SampleID}{labcode} }, ## convert adress   $Adresses{}
			 $metadata{$SampleID}{Host},
			 $metadata{$SampleID}{Gender},
			 $metadata{$SampleID}{'Patient age'},
			 $metadata{$SampleID}{'Patient status'},
			 $metadata{$SampleID}{'Sequencing technology'},
			 $Locations{ $metadata{$SampleID}{'regioncode'} }, ## convert
			 $Adresses{ $metadata{$SampleID}{'regioncode'} }, ## convert
			 "Original",
		         "testauthours") , "\n"  ); # FIXME



}
close GISAID;
#                          date($collection_date),
#                          location($origin),
#                          "", # Additonal location (e.g. Cruise Ship, Convention, Live animal market)
# 		 $config{host},
#                          "", # Additional host information (e.g. Patient infected while traveling in...)
#                          gender($gender),
#                          age($age),
#                          patient_status(),
#                          specimen_source(),
#                          outbreak(),
#                          last_vaccinated(),
#                          treatment(),
# 		 $config{sequencing_technology}, # FIXME: Take as argument
# 		 $config{assembly_method},
#                          coverage($IN_DIR, $sample_id),
# 		 $config{originating_lab},
# 		 $config{originating_lab_address},
#                          "", # Sample ID given by the sample provider
# 		 $config{submitting_lab},
# 		 $config{submitting_lab_address},
#                          "", # Sample ID given by the submitting lab
#		 $config{authors})
##    ) . "\n" . );





### fohm komplettering NB : does not handle multiple centers in the same metadatafile - solve sync issue woth GISAID accessions

open(CSV, ">". $out_dir .'/'. $GlobRegion  ."_". $GlobLab  ."_".$date."_komplettering.csv");
print CSV "provnummer,urvalskriterium,pangolin,GISAID_accession\n";
foreach my $SampleID (keys % metadata){
print CSV $metadata{$SampleID}{LibID} . "," . $metadata{$SampleID}{SelectionCriteria} . "," . "" . "\n"; 

}
close CSV;

##

open(CSV, ">". $out_dir .'/'. $GlobRegion  ."_". $GlobLab  ."_".$date."_pangolin_classification_format3.txt");
print CSV "taxon\tlineage\tconflict\tambiguity_score\tscorpio_call\tscorpio_support\tscorpio_conflict\tversion\tpangolin_version\tpangoLEARN_version\tpango_version\tstatus\tnote\n";
foreach my $SampleID (keys % metadata){
    print CSV $metadata{$SampleID}{LibID} . "\t" . $Pangolin_data{$SampleID}{lineage}  . "\t" . $Pangolin_data{$SampleID}{conflict} . "\t" . $Pangolin_data{$SampleID}{ambiguity_score} . "\t". $Pangolin_data{$SampleID}{scorpio_call} . "\t". $Pangolin_data{$SampleID}{scorpio_support} . "\t". $Pangolin_data{$SampleID}{scorpio_conflict} . "\t". $Pangolin_data{$SampleID}{version} . "\t". $Pangolin_data{$SampleID}{pangolin_version} . "\t". $Pangolin_data{$SampleID}{pangoLEARN_version} . "\t". $Pangolin_data{$SampleID}{pango_version} . "\t". $Pangolin_data{$SampleID}{status} . "\t" . $Pangolin_data{$SampleID}{note} . "\n";
    ##taxon lineage conflict ambiguity_score scorpio_call scorpio_support scorpio_conflict version pangolin_version pangoLEARN_version pango_ver pango_version status note
}
close CSV;




##foreach my $vcf_fn ( @vcfs ) {
#     my ($sample_id) = (split /\./, basename($vcf_fn))[0];
#     next if $sample_id eq "NTC" or $sample_id =~ "No_Sample" or $sample_id eq "No" or $sample_id =~ /NegativeControl/;
#     ##my $mlu_id = $metadata{$sample_id}->{SID};
#     my $mlu_id = $metadata{$sample_id}->{InternalLabId};
    
#     # Parse QC data
#     my $qc_data;
#     if( -e "$in_dir/$sample_id.qc.csv") {
# 	$qc_data = read_qc_data("$in_dir/$sample_id.qc.csv");
#     } else {
# 	die "QC data not found for $sample_id!\n";
#     }

#     # Check if QC passed
#     next if $qc_data->{pct_N_bases} > 5;

#     ##print CSV "$mlu_id,".($metadata{$sample_id}->{Urval} or "Information saknas").",".$gisaid_ids{$id_conversion_table{$sample_id}}."\n";
#     print CSV "$mlu_id,".($metadata{$sample_id}->{SelectionCriteria} or "Information saknas").",". "GISAID!!!" ."\n";
    
#     my $fa_fn = "$in_dir/$sample_id.consensus.fa";
#     die "No fasta for $sample_id" unless -e $fa_fn;

#     my $fq1 = "$in_dir/${sample_id}_subsample_R1_001.fastq.gz";
#     my $fq2 = "$in_dir/${sample_id}_subsample_R2_001.fastq.gz";
#     die "Fastq files missing for $sample_id!" if( ! -e "$in_dir/${sample_id}_R1_001.fastq.gz" or ! -e "$in_dir/${sample_id}_R2_001.fastq.gz" );

#     copy_and_fix_fasta($fa_fn, "${prefix}_${mlu_id}.consensus.fasta", "${prefix}_${mlu_id}");
#     system "ln -sf $vcf_fn ${prefix}_${mlu_id}.vcf";
#     system "ln -sf $fq1 ${prefix}_${mlu_id}_1.fastq.gz";
#     system "ln -sf $fq2 ${prefix}_${mlu_id}_2.fastq.gz";
    
# }

##reformat_pangolin($pangolin_fn, "${prefix}_${date}_pangolin_classification_".($config{pangolin_format} == 2 ? "format2" : "").".txt", \%metadata);


#############################################################################################
#############################################################################################
#############################################################################################

sub copy_and_fix_fasta {
    my ($orig_file, $new_file, $new_id) = @_;
    open(my $orig_fh, $orig_file) or die "cannot read: $orig_file\n";
    open(my $new_fh, '>'.$new_file) or die "cannot create file: $new_file\n";
    while(<$orig_fh>) {
	if( /^>/ ) {
	    print $new_fh ">$new_id\n";
	}
	else {
	    print $new_fh $_;
	}
    }
    close $new_fh;
    close $orig_fh;
}

sub reformat_pangolin {
    my ($in_fn, $out_fn, $metadata) = @_;
    open(my $in, $in_fn);
    open(my $out, ">".$out_fn);
    my $header = <$in>;
    print $out $header;
    while(my $line = <$in>) {
	my ($old_id) = ($line =~ /^.*?_(.*?)\./);
	my $new_id;
	if($metadata->{$old_id}->{SID}) {
	    my $new_id =  $metadata->{$old_id}->{SID};
	} else {
	    print STDERR "No SID (MLU ID) found for $old_id! Removing from pangolin file\n";
	    next;
	}
    
	$line =~ s/^.*?_(.*?)\..*?,/$new_id,/; # Remove the stuff surrounding the ID.
	print $out $line;
    }
    close $in;
    close $out;
}
    
sub read_qc_data{
    my $fn = shift;
    my @data = read_csv($fn, ',');
    return $data[0];

}

sub read_csv {
    my $fn = shift;
    my $sep = shift;
    open (my $fh, $fn);
    chomp(my $header = <$fh>);
    $header =~ s/\r//;
    my @header = split /$sep/, $header;
    
    my @data;
    while(<$fh>) {
	chomp;
	s/\r//;
	my @a = split /$sep/;
	my %entry;
	for my $i (0..$#header) {
	    $entry{$header[$i]} = $a[$i];
	}
	push @data, \%entry;
    }
    return @data;
    
}

sub read_config {
    my $fn = shift;
    open(my $fh, $fn) or die;
    my %config;
    while(<$fh>) {
	chomp;
	my ($key, $value) = split /=/;
	$config{$key} = $value;
    }
    close $fh;
    return %config;
}

sub read_conversion_table {
    my $fn = shift;
    open(my $fh, $fn);
    my %table;
    while(<$fh>) {
	chomp;
	my ($pseudo, $real) = split /\t/;
	$table{$real} = $pseudo;
    }
    return %table;
}

sub read_gisaid_ids {
    my $fn = shift;
    open(my $fh, $fn);
    my %table;
    while(<$fh>) {
	next if /^submissions /;
	chomp;
	my ($name, $gisaid) = split '; ';
	my @a = split /\//, $name;
	my $pseudo_id = $a[2];
	$pseudo_id =~ /(0)*(\d+)$/;
	my $pseudo_num = $2;
	$table{$pseudo_num} = $gisaid;
    }
    return %table;
}

## sed $'1s/\xef\xbb\xbf//' < FohmAndGisaidMinimumExamplelim3.csv > FohmAndGisaidMinimumExamplelim3NoBom.csv

sub read_pat_metadata {
    my $fn = shift;
    ##$fn =~ s/^\x{FEFF}//;
    $fn =~ s/\N{U+FEFF}//;
    my @csv = read_csv("iconv -f iso-8859-1 -t UTF-8 '$fn'|", ";");
    my %csv;
    foreach my $entry (@csv) {
	##$csv{$entry->{Labbnummer}} = $entry;
	##print $entry->{InternalLabId} . "\n";
##	print "$_\n" for keys $entry;
	$csv{$entry->{InternalLabId}} = $entry;
    }
    return %csv;
}

sub read_PipelineResults {
    my $fn = shift;
    ##$fn =~ s/^\x{FEFF}//;
    $fn =~ s/\N{U+FEFF}//;
    my @tsv = read_csv("iconv -f iso-8859-1 -t UTF-8 '$fn'|", "\t");
    my %PipeRes;
    foreach my $entry (@tsv) {
	##$csv{$entry->{Labbnummer}} = $entry;
        ##print $entry->{taxon} . "\n";
	my $FixedID = $entry->{taxon};
	$FixedID =~ s/\/.*//;
	##print $FixedID  . "\n";
##      print "$_\n" for keys $entry;
        $PipeRes{$FixedID} = $entry;
    }
    return %PipeRes;
}
