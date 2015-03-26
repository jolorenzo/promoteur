#!/usr/bin/perl

=pod

=head1 NAME

[gff3fna2gff32bed2fasta.pl - short description]

=head1 SYNOPSIS

    qsub -q bioinfo.q -b yes -V -N format1 perl gff3fna2gff32bed2fasta.pl <gff3_input_file> <genome_sequence_file> <fna_output_file> -bed=all|gene|mRNA|polypeptide|... -type=gene|mRNA|polypeptide|... -begin=0|-x|+x -end=0|-x|+x

=head1 REQUIRES

[Perl5.004, POSIX, Win32] #+++

=head1 DESCRIPTION

#--- What this script can do, how to use it, ...
[Script description]. #+++

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Getopt::Long;
use Pod::Usage;
#--- when using try {...} catch XXX::YYY with {...} otherwise {...};
#+++ use Error qw(:try);
#--- when working with files
use Fatal qw/:void open close/;

use Bio::Tools::GFF;
use Bio::FeatureIO::bed;
use List::MoreUtils qw(uniq);
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::Fasta;
use JSON;
use Statistics::Descriptive;
use Bio::DB::SeqFeature;
use Env qw(HOME);
use Data::Dumper;

# Script global constants
##########################

=pod

=head1 CONSTANTS

B<$CONSTANT_NAME>: ([constant nature]) #+++

[constant description and use]. #+++

B<CONSTANT_NAME2>: (sub returning a [constant return nature]) #+++

[constant description and use]. #+++

#--- Example:
#--- B<PI>: (sub returning a real)
#---
#--- used by trigonometrical functions;
#--- can also be used to define an angle value.
#---
#--- ...
#---
#--- sub PI() { return (4 * atan2(1, 1)); }
#---

=cut

#+++ my [$CONSTANT_NAME] = ["value"];
#+++ sub [CONSTANT_NAME2]() { return ["value2"]; }


# Script global variables
##########################

=pod

=head1 VARIABLES

B<variable_name>: ([variable nature]) #+++

[variable description, use and default value]. #+++

#--- Example:
#--- B<$output_method>: (integer)
#---
#--- used to store current output method;
#--- 0=text on screen (default), 1=graphic on screen, 2=text in a file, 3=graphic to the printer.
#---
#---     ...
#---
#--- my $output_method = 0;
#---

=cut

my $gff3_file=$ARGV[0];  #gff3 input file
my $genome=$ARGV[1]; #genome sequence directory
my $fna_output_directory= $1 if $ARGV[2] =~ /(.*\/)[^\/]*\.?.*/i;
my $fna_output_file = $ARGV[2];
my $bed_output_directory= $1 if $ARGV[3] =~ /(.*\/)[^\/]*\.?.*/i;
my $bed_output_file = $ARGV[3];

my $title;
my $title2;
my $title3;

open my $FH_IN, "<", $gff3_file or die "Could not open file '$gff3_file' $!";
while (my $ligne = <$FH_IN>) {     # lit le fichier en entrée ligne par ligne
     if (($ligne !~ /^##/i) && ($ligne =~ /^[A-Z]{5}/i)){     
         my @mots = split('\t', $ligne);
         if (($mots[0] =~ /^([A-Z]{5})/i) && ($mots[1] !~ /manual_curation/i)){
        	$title = "$1-"."$mots[1]-sequence_feature" ;
         	$title2 = "$1-"."$mots[1]";
        	$title3 = "$1";
         }
     }     
}
close $FH_IN;

# Script global functions
##########################

=pod

=head1 FUNCTIONS

=head2 [subName] #+++

B<Description>: [function description]. #+++

B<ArgsCount>: [count of arguments] #+++

=over 4

=item variable_name: ([variable nature]) ([requirements]) #+++ see below

#--- requirement can be:
#--- (R)=required,
#--- (O)=optional
#--- (U)=optional and must be undef if omitted
[variable description]. #+++

=item variable_name2: ([variable nature]) ([requirements]) #+++

[variable description]. #+++

=back

B<Return>: ([return type]) #+++

[return description]. #+++

B<Exception>:

=over 4

[=item * exception_type:

description, case when it occurs.] #+++ see below

=back

#--- Example:
#---
#---=over 4
#---
#---=item * Range error:
#---
#---thrown when "argument x" is outside the supported range, for instance "-1".
#---
#---=back

B<Example>:

    [example code] #+++

=cut

#+++sub [subName] #+++
#+++{
#+++    my ([...]) = @_; #+++ add missing arguments
#--- if needed, parameters check:
#+++    # parameters check
#+++    if (0 != @_)
#+++    {
#+++        confess "usage: subName();"; #+++
#+++    }
#+++}
my ($type, $begin, $end, $bed);
sub sortGff3{

	my $gff = new Bio::Tools::GFF(
	-file => $gff3_file,
	-gff_version => 3
	);
	
	my $cpt = 0;
	my (%feature, $id, $mrna_id, %gene_refseq, %mrna, %polypeptide, %other, %other_plus, %cds, %transposable_gene, $transp_id);
	
	while(my $feature = $gff->next_feature) {
		if ($feature->primary_tag() eq "gene") {
			$cpt++;
			#last if $cpt ==20;
			($id) = $feature->get_tag_values("ID");	
			$feature{$id} = $feature;			
			push @{$gene_refseq{$feature->seq_id}} , $feature;
		}
		elsif ($feature->primary_tag() eq "pseudogene"){	
			push @{$gene_refseq{$feature->seq_id}} , $feature;		
		}
		elsif ($feature->primary_tag() eq "transposable_element"){	
			push @{$gene_refseq{$feature->seq_id}} , $feature;		
		}
		elsif ($feature->primary_tag eq "transposable_element_gene") {		 
			push @{$gene_refseq{$feature->seq_id}} , $feature;			
		}
		elsif ($feature->primary_tag eq "pseudogenic_transcript") {
			($mrna_id) = $feature->get_tag_values("ID");
			($id)  = $feature->get_tag_values("Parent"); 
			$mrna{$id}{$mrna_id} = $feature;
		}
		elsif ($feature->primary_tag eq "mRNA") {
			($mrna_id) = $feature->get_tag_values("ID");
			($id)  = $feature->get_tag_values("Parent"); 
			$mrna{$id}{$mrna_id} = $feature;
		}
		elsif ($feature->primary_tag =~ /polypeptide|protein/ ) { 
			($mrna_id) = $feature->get_tag_values("Derives_from"); 
			push @{$polypeptide{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "exon") { 	 
			($mrna_id) = $feature->get_tag_values("Parent");  
			push @{$other{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "pseudogenic_exon") {
			($mrna_id) = $feature->get_tag_values("Parent");  
			push @{$other{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "CDS") { 	 
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other{$mrna_id}}, $feature;
			push @{$cds{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag =~/.+_prime_UTR/i) { 	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other{$mrna_id}}, $feature;
		}
		elsif (($feature->primary_tag =~/.+RNA/i) && ($feature->primary_tag !~/^mRNA/i)) { 	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other{$mrna_id}}, $feature;
		}		
		elsif ($feature->primary_tag eq "transposon_fragment") { 	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other{$mrna_id}}, $feature;
		}		
	}
	$gff->close;
	
	my $in = new Bio::SeqIO(
		-file => $genome,
		-format => "fasta"
	);
	
	my %seq; 
	while(my $seqobj = $in->next_seq){
		$seq{$seqobj->display_id} = $seqobj;
	}
	$in->close;
	
	my $file_cds = $fna_output_directory.$title2."-CDS-genfam.fna";
	my $file_prot = $fna_output_directory.$title2."-polypeptide-genfam.faa";
	my $file_exon = $fna_output_directory.$title2."-exon-genfam.fna";
	my $seq_cds;
	my $seq_pep;
	my $seq_exon;
	$seq_cds = new Bio::SeqIO(
		-file => ">$file_cds",
		-format => "fasta"
	);
	$seq_pep = new Bio::SeqIO(
		-file => ">$file_prot",
		-format => "fasta"
	);
	$seq_exon = new Bio::SeqIO(
		-file => ">$file_exon",
		-format => "fasta"
	);
	
	my $outfile = $fna_output_directory.$title."-locus_tag-genfam.gff3";
	print "Print data to $outfile...\n";
	
	my $out = new Bio::Tools::GFF(
		-file => ">$outfile",
		-gff_version => 3
	);
	
	my %h_id = ( "ARATH" => "At", "BRADI" => "Bd", "GLYMA" => "Gm", "GOSRA"   => "Gr", "LOTJA"   => "Lj", "MEDTR"   => "Mt", "MUSAC"   => "Ma",
	"ORYSI"   => "Os", "ORYSJ"   => "Os", "POPTR"   => "Pt", "SOLLC"   => "Sl", "SORBI"   => "Sb", "THECC"   => "Tc", "VITVI"   => "Vv", 
	"MAIZE"   => "Zm", "MALDO"   => "Md", "MANES"   => "Me", "RICCO"   => "Rc", "SETIT"   => "Si", "SOLTU"   => "St");	

	my %locus;
	my %function;
	my $source_tag = "manual_curation";
	my $is_complete = 0; 
	my $is_incomplete=0;
	my $codon_start=0;
	my $codon_stop=0;
	
	foreach my $seq_id (sort {$a cmp $b} keys%gene_refseq){
		my $count = 0;
		my $seqobj;
		
		if (defined $seq{$seq_id}){
			$seqobj = $seq{$seq_id};
		}
		else{
			$seqobj = 0;
		}
		foreach my $gene (sort {$a->start <=> $b->start} @{$gene_refseq{$seq_id}}){
			if ($gene->primary_tag() eq "gene"){
				my ($name) = $gene->get_tag_values("ID");	
				$count++;		
				my $gene_id = sprintf( "_g%05d", $count * 10 );
				my $first_poly_id = sprintf( "_p%05d", $count * 10 );
				my $first_mrna_id = sprintf( "_t%05d", $count * 10 );
				my $start_gene = $gene->start;
				my $end_gene = $gene->end;
				my ($function) = $function{$name};
	
				my ($L5, $chr);
				if ($gene->seq_id =~ /^([A-Z]{5})(.+)/i) {
					$L5 = $1;
					$chr = $2;
				} elsif($gene->seq_id =~ /^([A-Z]{5})_scaffold(.+)/i){
					$L5 = $1;
					$chr = $2;
				} elsif($gene->seq_id =~ /^([A-Z]{5})_contig(.+)/i){
					$L5 = $1;
					$chr = $2;
				}
	
				my $L2 = $h_id{$L5};
	
	
				my ($Gene_id) =  $gene->get_tag_values("ID");
				my $locus_tag = $L2.$chr.$gene_id."_$L5";

				$locus{$Gene_id} = {   					
					locus_tag => $locus_tag,
					count => $count,
				};

				$gene->add_tag_value("locus_tag", $L2.$chr.$gene_id."_$L5");
	
				$out->write_feature($gene);
		
				my $count_mrna = 0;
				if (exists $mrna{$name}){
					foreach my $mrna (keys%{$mrna{$name}}){
						$count_mrna++;
						my $mrna_id = $first_mrna_id ."." . $count_mrna;
						my $poly_id = $first_poly_id ."." . $count_mrna;
						my $feat_mrna = $mrna{$name}{$mrna};
						my ($old_mrna_id) =  $feat_mrna->get_tag_values("ID");					
						my $start_poly = 100000000000000000;
						my $end_poly   = -1; 
						my $start_poly_exon = 100000000000000000;
						my $end_poly_exon   = -1; 
						my @cds;
						my @exon;
						my $protein_id;
						foreach my $other (@{$other{$mrna}}){		
							if ($other->primary_tag() eq "CDS") {
								($protein_id) = $other->get_tag_values("protein_id") if $other->has_tag("protein_id");
								$start_poly = $other->start if $other->start < $start_poly;
								$end_poly   = $other->end if $other->end > $end_poly;
								push @cds,$other;
							}			
							if ($other->primary_tag() eq "exon") { 
								push @exon,$other;
								$start_poly_exon = $other->start if $other->start < $start_poly_exon;
								$end_poly_exon   = $other->end if $other->end > $end_poly_exon;
							}	
						}
						my $strand = $feat_mrna->strand;
	
						if ($feat_mrna->seq_id =~ /^([A-Z]{5})(.+)/i){
							$L5 = $1;
							$chr = $2;
						}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_scaffold(.+)/i){
							$L5 = $1;
							$chr = $1;
						}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_contig(.+)/i){
							$L5 = $1;
							$chr = $2;
						}
			
						$feat_mrna->add_tag_value("locus_tag", $L2.$chr.$mrna_id."_$L5");
						$out->write_feature($feat_mrna);
				
						unless (@cds ) {
							@cds = @exon;
							$start_poly = $start_poly_exon;
							$end_poly   = $end_poly_exon;
						}
						if (exists $polypeptide{$mrna}){
							foreach my $poly (@{$polypeptide{$mrna}}){
	
								if ($poly->seq_id =~ /^([A-Z]{5})(.+)/i){
									$L5 = $1;
									$chr = $1;
								}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(.+)/i){
									$L5 = $1;
									$chr = $1;
								}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(.+)/i){
									$L5 = $1;
									$chr = $2;
								}
								$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
								$out->write_feature($poly);
							}
						}
						else {
							my $poly = new Bio::SeqFeature::Generic(
								-seq_id 	=> $seq_id,
								-source_tag => $source_tag,
								-primary_tag => 'polypeptide',
								-start       => $start_poly,
								-end         => $end_poly,
								-strand      => $strand,
								-tag 		 => {
									ID	=> $feat_mrna->get_tag_values("ID")."-protein",
									Name	=> $feat_mrna->get_tag_values("ID"),
									Derives_from => $feat_mrna->get_tag_values("ID")
								}
							); 
							if ($feat_mrna->has_tag("Dbxref")) {
								my @dbxref = $feat_mrna->get_tag_values("Dbxref");
								foreach my $dbxref (@dbxref) {
									$poly->add_tag_value("Dbxref",$dbxref);
								}
							}
							if ($poly->seq_id =~ /^([A-Z]{5})(.+)/i){
								$L5 = $1;
								$chr = $1;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(.+)/i){
								$L5 = $1;
								$chr = $1;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(.+)/i){
								$L5 = $1;
								$chr = $2;
							}
							$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
						
							$out->write_feature($poly);
						}
			
						foreach my $other (@{$other{$old_mrna_id}}){
							$out->write_feature($other);						
						}
						
						if ($seqobj != 0){
							if ($cds[0]->strand == 1) {
							@cds = sort{$a->start <=> $b->start} @cds;
							}
							else {
								@cds = sort{$b->start <=> $a->start} @cds;
							}
							my $cds;
							foreach my $feature (@cds) {	
								my $seqtrunc_cds = $seqobj->trunc($feature->start,$feature->end);	
								if ($feature->strand == -1) {
									$cds .= $seqtrunc_cds->revcom()->seq;
								}
								else {	
									$cds .= $seqtrunc_cds->seq;
								}
							}
							my $seqobj_cds = Bio::PrimarySeq->new (
								-display_id  => $old_mrna_id."_$title3",
								-seq         => $cds ,
								-desc		=> $function
							);  
					
							if ($exon[0]->strand == 1) {
								@exon = sort{$a->start <=> $b->start} @exon;
							}
							else {
								@exon = sort{$b->start <=> $a->start} @exon;
							}
							my $exon;
							foreach my $feature (@exon) {	
								my $seqtrunc_exon = $seqobj->trunc($feature->start,$feature->end);	
								if ($feature->strand == -1) {
									$exon .= $seqtrunc_exon->revcom()->seq;
								}
								else {	
									$exon .= $seqtrunc_exon->seq;
								}
							}
							my $seqobj_exon = Bio::PrimarySeq->new (
								-display_id  => $old_mrna_id."_$title3",
								-seq         => $exon ,
								-desc		=> $function
							);   
							if ($cds =~ /^ATG.*/ &&  $cds =~ /.*(TAG|TAA|TGA)$/){
								$is_complete++;
							}
							else {
								if ($cds =~ /^ATG.*/){
									$codon_stop++;
								}
								elsif ($cds =~ /.*(TAG|TAA|TGA)$/){
									$codon_start++;
								}
								else {
									$is_incomplete++;
									print $name ,"\t", $old_mrna_id,"\n";
								}
							}
							my $seqobj_pep = $seqobj_cds->translate();
							$seqobj_pep->display_id($old_mrna_id."_$title3");
							$seq_cds->write_seq($seqobj_cds);
							$seq_pep->write_seq($seqobj_pep);
							$seq_exon->write_seq($seqobj_exon);
						}
					}
				}
				elsif(exists $other{$name}){
					foreach my $other (@{$other{$name}}){
						$out->write_feature($other);						
					}
				}
				else {
					$count_mrna++;
					my $strand = $gene->strand;
					my $mrna = new Bio::SeqFeature::Generic(
						-seq_id 	=> $seq_id,
						-source_tag => $source_tag,
						-primary_tag => 'mRNA',
						-start       => $start_gene,
						-end         => $end_gene,
						-strand      => $strand,
						-tag 		 => {
							ID	    => "$Gene_id\.$count_mrna",
							Name	=> "$Gene_id\.$count_mrna",
							Parent  => $Gene_id
						}
					);							
			
					my $mrna_id = $first_mrna_id ."." . $count_mrna;
					my $poly_id = $first_poly_id ."." . $count_mrna;
					my $feat_mrna = $mrna;
					my ($old_mrna_id) =  $feat_mrna->get_tag_values("Parent");
					my $start_poly = 100000000000000000;
					my $end_poly   = -1; 
					my $start_poly_exon = 100000000000000000;
					my $end_poly_exon   = -1; 
					my @cds;
					my @exon;
					my $protein_id;
					foreach my $other (@{$other{$old_mrna_id}}){		
						if ($other->primary_tag() eq "CDS") {
							($protein_id) = $other->get_tag_values("protein_id") if $other->has_tag("protein_id");
							$start_poly = $other->start if $other->start < $start_poly;
							$end_poly   = $other->end if $other->end > $end_poly;
							push @cds,$other;
						}			
						if ($other->primary_tag() eq "exon") { 
							push @exon,$other;
							$start_poly_exon = $other->start if $other->start < $start_poly_exon;
							$end_poly_exon   = $other->end if $other->end > $end_poly_exon;
						}	
					}
		
					if ($feat_mrna->seq_id =~ /^([A-Z]{5})(.+)/i){
						$L5 = $1;
						$chr = $2;
					}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_scaffold(.+)/i){
						$L5 = $1;
						$chr = $1;
					}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_contig(.+)/i){
						$L5 = $1;
						$chr = $2;
					}

					$feat_mrna->add_tag_value("locus_tag", $L2.$chr.$mrna_id."_$L5");
					$out->write_feature($feat_mrna);
					unless (@cds ) {
						@cds = @exon;
						$start_poly = $start_poly_exon;
						$end_poly   = $end_poly_exon;
					}
					if (exists $polypeptide{$old_mrna_id}){
						foreach my $poly (@{$polypeptide{$old_mrna_id}}){
											
							if ($poly->seq_id =~ /^([A-Z]{5})(.+)/i){
								$L5 = $1;
								$chr = $1;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(.+)/i){
								$L5 = $1;
								$chr = $1;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(.+)/i){
								$L5 = $1;
								$chr = $2;
							}
							$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
							$out->write_feature($poly);
						}
					}
					else {
						my $poly = new Bio::SeqFeature::Generic(
							-seq_id 	=> $seq_id,
							-source_tag => $source_tag,
							-primary_tag => 'polypeptide',
							-start       => $start_poly,
							-end         => $end_poly,
							-strand      => $strand,
							-tag 		 => {
								ID	=> "$Gene_id\.$count_mrna"."-protein",
								Name	=> "$Gene_id\.$count_mrna",
								Derives_from => "$Gene_id\.$count_mrna",
							}
						); 
						if ($feat_mrna->has_tag("Dbxref")) {
							my @dbxref = $feat_mrna->get_tag_values("Dbxref");
							foreach my $dbxref (@dbxref) {
								$poly->add_tag_value("Dbxref",$dbxref);
							}
						}
						if ($poly->seq_id =~ /^([A-Z]{5})(.+)/i){
							$L5 = $1;
							$chr = $1;
						}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(.+)/i){
							$L5 = $1;
							$chr = $1;
						}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(.+)/i){
							$L5 = $1;
							$chr = $2;
						}
						$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
				
						$out->write_feature($poly);
					}
					foreach my $other (@{$other{$old_mrna_id}}){
						$out->write_feature($other);						
					}
					
					if ($seqobj != 0){
						if ($cds[0]->strand == 1) {
							@cds = sort{$a->start <=> $b->start} @cds;
						}
						else {
							@cds = sort{$b->start <=> $a->start} @cds;
						}
						my $cds;
						foreach my $feature (@cds) {	
							my $seqtrunc_cds = $seqobj->trunc($feature->start,$feature->end);	
							if ($feature->strand == -1) {
								$cds .= $seqtrunc_cds->revcom()->seq;
							}
							else {	
								$cds .= $seqtrunc_cds->seq;
							}
						}
						my $seqobj_cds = Bio::PrimarySeq->new (
							-display_id  => $old_mrna_id."_$title3",
							-seq         => $cds ,
							-desc		=> $function
						);  
				
						if ($exon[0]->strand == 1) {
							@exon = sort{$a->start <=> $b->start} @exon;
						}
						else {
							@exon = sort{$b->start <=> $a->start} @exon;
						}
						my $exon;
						foreach my $feature (@exon) {	
							my $seqtrunc_exon = $seqobj->trunc($feature->start,$feature->end);	
							if ($feature->strand == -1) {
								$exon .= $seqtrunc_exon->revcom()->seq;
							}
							else {	
								$exon .= $seqtrunc_exon->seq;
							}
						}
						my $seqobj_exon = Bio::PrimarySeq->new (
							-display_id  => $old_mrna_id."_$title3",
							-seq         => $exon ,
							-desc		=> $function
						);   
						if ($cds =~ /^ATG.*/ &&  $cds =~ /.*(TAG|TAA|TGA)$/){
							$is_complete++;
						}
						else {
							if ($cds =~ /^ATG.*/){
								$codon_stop++;
							}
							elsif ($cds =~ /.*(TAG|TAA|TGA)$/){
								$codon_start++;
							}
							else {
								$is_incomplete++;
								print $name ,"\t", $old_mrna_id,"\n";
							}
						}
						my $seqobj_pep = $seqobj_cds->translate();
						$seqobj_pep->display_id($old_mrna_id."_$title3");
						$seq_cds->write_seq($seqobj_cds);
						$seq_pep->write_seq($seqobj_pep);
						$seq_exon->write_seq($seqobj_exon);	
					}		
				}							
			}
			elsif  ($gene->primary_tag() eq "pseudogene"){
				my ($name) = $gene->get_tag_values("ID");
				$out->write_feature($gene);
				foreach my $mrna (keys%{$mrna{$name}}){
					my $feat_mrna = $mrna{$name}{$mrna};
					my ($old_mrna_id) =  $feat_mrna->get_tag_values("ID");	
					$out->write_feature($feat_mrna);						
					foreach my $other (@{$other{$old_mrna_id}}){
						$out->write_feature($other);						
					}				
				}
			}
			elsif ($gene->primary_tag() eq "transposable_element_gene"){
				my ($name) = $gene->get_tag_values("ID");
				$out->write_feature($gene);
				foreach my $mrna (keys%{$mrna{$name}}){
					my $feat_mrna = $mrna{$name}{$mrna};
					my ($old_mrna_id) =  $feat_mrna->get_tag_values("ID");	
					$out->write_feature($feat_mrna);						
					foreach my $other (@{$other{$old_mrna_id}}){
						$out->write_feature($other);						
					}				
				}
			}
			elsif (($gene->primary_tag() eq "transposable_element")){
				my ($name) = $gene->get_tag_values("ID");
				$out->write_feature($gene);						
				foreach my $other (@{$other{$name}}){
					$out->write_feature($other);						
				}				
			}
		}			
	}
	my $json = encode_json \%locus;
	my $file = $fna_output_directory.$title2."-locus_tag.json";
	unless(open FILE, '>'.$file) {
		die "Unable to create $file";
	}
	print FILE $json;
	close FILE;
	$out->close();
}

sub gff3tobed{

	my $new_file = $fna_output_directory.$title."-locus_tag-genfam.gff3";
	my $gffio = Bio::Tools::GFF -> new(-file =>$new_file , -gff_version => 3);

	my $strand=0;
	my $strand1='-';
	
	##### Liste les fichiers ayant l'extension .bed 
	my @list = glob("*.bed"); 
  
	### recupere le nombre de fichier 
	my $numberoffile = scalar(@list);
	
	for ( my $v = 0; $v < $numberoffile; $v++ ) {
		unlink $list[$v]; #supprime les fichiers bed
	}
	
	my $feat_start;
	my $feat_end;
	while(my $feature = $gffio->next_feature()) {
		if(($feature->primary_tag =~/gene/i) && ($bed =~/gene|all/i)){
			
			my $file = $fna_output_directory.$title."-gene-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
	
		if(($feature->primary_tag =~/mRNA/i) && ($bed =~/mRNA|all/i)){

			my $file = $fna_output_directory.$title."-mRNA-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
	
		if(($feature->primary_tag =~/five_prime_UTR/i) && ($bed =~/five_prime_UTR|all/i)){

			my $file = $fna_output_directory.$title."-five_prime_UTR-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
	
		if(($feature->primary_tag =~/three_prime_UTR/i) && ($bed =~/three_prime_UTR|all/i)){

			my $file = $fna_output_directory.$title."-three_prime_UTR-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
		
		if(($feature->primary_tag =~/CDS/i) && ($bed =~/CDS|all/i)){

			my $file = $fna_output_directory.$title."-CDS-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
	
		if(($feature->primary_tag =~/exon/i) && ($bed =~/exon|all/i)){

			my $file = $fna_output_directory.$title."-exon-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
	
		if(($feature->primary_tag =~/intron/i) && ($bed =~/intron|all/i)){

			my $file = $fna_output_directory.$title."-intron-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
	
		if(($feature->primary_tag =~/polypeptide|protein/i) && ($bed =~/polypeptide|protein|all/i)){

			my $file = $fna_output_directory.$title."-polypeptide-genfam.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
			close FILE;
		}
	}
	$gffio->close();
  
}

sub gff3tobed_before{

	my $new_file = $gff3_file;
	my $strand=0;
	my $strand1='-';
	
	my $gff = new Bio::Tools::GFF(
	-file => $new_file,
	-gff_version => 3
	);
	
	my $cpt = 0;
	my $cpt1 = 0;
	my $cpt2 = 0;
	my $cpt3 = 0;
	my $cpt4 = 0;
	my $cpt5 = 0;
	
	my (%gene, %mRNA, %polypeptide, %CDS, %exon, %five_prime_UTR, $mrna_id, $id);
	
	while(my $feature = $gff->next_feature) {
		if (($type =~ /gene/i) && ($feature->primary_tag() eq "gene")) {
			$cpt++;			
			($id) = $feature->get_tag_values("ID");
			$gene{$feature->seq_id}{$cpt} = $feature;			
		}
		elsif (($type =~ /mRNA/i) && ($feature->primary_tag eq "mRNA")) {
			$cpt1++;
			($mrna_id) = $feature->get_tag_values("ID");
			($id)  = $feature->get_tag_values("Parent");
			$mRNA{$id}{$mrna_id} = $feature;
			#$mRNA{$cpt1} = $feature;
		}
		elsif (($type =~ /polypeptide|protein/i) && ($feature->primary_tag =~ /polypeptide|protein/i)) { 
			$cpt2++;			 
			($mrna_id) = $feature->get_tag_values("Derives_from"); 
			push @{$polypeptide{$mrna_id}}, $feature;
			#$polypeptide{$cpt2} = $feature;
		}
		elsif (($type =~ /exon/i) && ($feature->primary_tag eq "exon")) { 	 
			$cpt3++;
			($mrna_id) = $feature->get_tag_values("Parent");  
			push @{$exon{$mrna_id}}, $feature;
			#$exon{$cpt3} = $feature;
		}
		elsif (($type =~ /CDS/i) && ($feature->primary_tag eq "CDS")) {
			$cpt4++;
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$CDS{$mrna_id}}, $feature;
			#$CDS{$cpt4} = $feature;
		}
		elsif (($type =~ /five_prime_UTR/i) && ($feature->primary_tag =~/five_prime_UTR/i)) { 	
			$cpt5++;
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$five_prime_UTR{$mrna_id}}, $feature;
			#$five_prime_UTR{$cpt5} = $feature;
		}
	}
	$gff->close;

	my $outfile = $bed_output_file;
	print "Print data to $outfile...\n";
	
	open(FILE, ">$outfile") or die "Could not open file '$outfile' $!";
	
	if ($type =~ /gene/i){ #à modifier
		foreach my $seq_id (sort {$a cmp $b} keys%gene){
			foreach my $gene (sort {$a cmp $b} keys%{$gene{$seq_id}}){			
				my ($feature) = $gene{$seq_id}{$gene};
				my $score = ".";					
				#Cas du premier gène du chromosome		
				if((not exists $gene{$seq_id}{$gene-1}) && ($feature->start+$begin > 0)){
					#print "1";
					
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}			
					print FILE join("\t",$feature->seq_id,$feature->start-1+$begin,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas du premier gène du chromosome
				elsif((not exists $gene{$seq_id}{$gene-1}) && ($feature->start+$begin <= 0)){
					#print "2";
				
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}
					print FILE join("\t",$feature->seq_id,$feature->start-$feature->start,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas du gène dans un autre gène (intron)
				elsif(($gene{$seq_id}{$gene-1}->start >= $feature->start+$begin) && ($gene{$seq_id}{$gene-1}->end >= $feature->start)){	#Exon à utiliser
					#print "3";
					
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}					
					print FILE join("\t",$feature->seq_id,$gene{$seq_id}{$gene-1}->start,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas du gène dans un autre gène (intron)
				elsif(($gene{$seq_id}{$gene-1}->start < $feature->start+$begin) && ($gene{$seq_id}{$gene-1}->end >= $feature->start)){ #Exon à utiliser
					#print "4";
					
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->start-1+$begin,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes identique en start et placé en début de chromosome
				elsif ((not exists $gene{$seq_id}{$gene-2}) && ($gene{$seq_id}{$gene-1}->start == $feature->start) && ($feature->start+$begin > 0)){
					#print "5";
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->start-1+$begin,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes identique en start et placé en début de chromosome
				elsif ((not exists $gene{$seq_id}{$gene-2}) && ($gene{$seq_id}{$gene-1}->start == $feature->start) && ($feature->start+$begin <= 0)){
					#print "6";
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->start-$feature->start,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes identique en start et région non codante entre le gene -2 et notre gène
				elsif ((exists $gene{$seq_id}{$gene-2}) && ($gene{$seq_id}{$gene-1}->start == $feature->start) && ($gene{$seq_id}{$gene-2}->end < $feature->start+$begin)){
					print "7";
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->start+$begin,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes identique en start et région non codante entre le gene -2 et notre gène
				elsif ((exists $gene{$seq_id}{$gene-2}) && ($gene{$seq_id}{$gene-1}->start == $feature->start) && ($gene{$seq_id}{$gene-2}->end >= $feature->start+$begin) ){
					#print "8";
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$gene{$seq_id}{$gene-2}->end,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes identique en start dans un autre gène, prendre l'exon
				elsif ((exists $gene{$seq_id}{$gene-2}) && ($gene{$seq_id}{$gene-1}->start == $feature->start) && ($gene{$seq_id}{$gene-2}->start < $feature->start+$begin) && ($gene{$seq_id}{$gene-2}->end > $feature->start)){
					#print "9";
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->start+$begin,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes identique en start dans un autre gène, prendre l'exon
				elsif ((exists $gene{$seq_id}{$gene-2}) && ($gene{$seq_id}{$gene-1}->start == $feature->start) && ($gene{$seq_id}{$gene-2}->end >= $feature->start+$begin) && ($gene{$seq_id}{$gene-2}->end > $feature->start)){
					#print "a";
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$gene{$seq_id}{$gene-2}->start,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes qui se suivent
				elsif(($gene{$seq_id}{$gene-1}->end >= $feature->start+$begin) && ($gene{$seq_id}{$gene-1}->end < $feature->start)){			
					#print "b";
					
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$gene{$seq_id}{$gene-1}->end,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}
				#Cas de 2 gènes qui se suivent				
				else{
					#print "c";
					
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->start-1+$begin,$feature->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
				}				
			}
		}
	}
	close FILE;
	
	system("bedtools getfasta -s -fi $genome -bed $outfile -fo $fna_output_file -name");
}

sub gff3tobed_after{

	my $new_file = $fna_output_directory.$title."-locus_tag-genfam.gff3";
	my $strand=0;
	my $strand1='-';
	
	##### Liste les fichiers ayant l'extension .bed 
	my @list = glob($fna_output_directory.$title."-".$type."_after-genfam.bed"); 
  
	### recupere le nombre de fichier 
	my $numberoffile = scalar(@list);
	
	for ( my $v = 0; $v < $numberoffile; $v++ ) {
		unlink $list[$v]; #supprime les fichiers
	}
	
	my $in = new Bio::SeqIO(
		-file => $genome,
		-format => "fasta"
	);
	
	my (%length, $seq_obj);
	while ($seq_obj = $in->next_seq){  
  		 $length{$seq_obj->primary_id} = $seq_obj->length();
	}
	$in->close;
	
	my $gff = new Bio::Tools::GFF(
	-file => $new_file,
	-gff_version => 3
	);
	
	my $cpt = 0;
	my $cpt1 = 0;
	my $cpt2 = 0;
	my $cpt3 = 0;
	my $cpt4 = 0;
	my $cpt5 = 0;
	
	my (%gene, %mRNA, %polypeptide, %CDS, %exon, %five_prime_UTR);
	
	while(my $feature = $gff->next_feature) {
		if (($type =~ /gene/i) && ($feature->primary_tag() eq "gene")) {
			$cpt++;
			$gene{$cpt} = $feature;			
		}
		elsif (($type =~ /mRNA/i) && ($feature->primary_tag eq "mRNA")) {
			$cpt1++;
			$mRNA{$cpt1} = $feature;
		}
		elsif (($type =~ /polypeptide/i) && ($feature->primary_tag eq "polypeptide")) { 
			$cpt2++; 
			$polypeptide{$cpt2} = $feature;
		}
		elsif (($type =~ /exon/i) && ($feature->primary_tag eq "exon")) { 	 
			$cpt3++;
			$exon{$cpt3} = $feature;
		}
		elsif (($type =~ /CDS/i) && ($feature->primary_tag eq "CDS")) {
			$cpt4++;
			$CDS{$cpt4} = $feature;
		}
		elsif (($type =~ /five_prime_UTR/i) && ($feature->primary_tag =~/five_prime_UTR/i)) { 	
			$cpt5++;
			$five_prime_UTR{$cpt5} = $feature;
		}
	}
	$gff->close;
	
	#my $var = join "\n", keys %gene;
	#print $var;
	
	my $gff1 = new Bio::Tools::GFF(
	-file => $new_file,
	-gff_version => 3
	);
	
	my $j = 0;	
	
	while(my $feature = $gff1->next_feature()) {		
		if ($type =~ /gene/i){
			if($feature->primary_tag =~/$type/i){
			$j++;
			print $feature->end;
				if((not exists $gene{$j+1}) && ($feature->end+$end <= $length{$feature->seq_id})){
					print "1";
					my $file = $fna_output_directory.$title."-".$type."_after-genfam.bed";
					unless(open FILE, '>>'.$file) {
						die "Unable to create $file";
					}
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}
					print FILE join("\t",$feature->seq_id,$feature->end+$begin,$feature->end+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
					close FILE;
				}
				elsif((not exists $gene{$j+1}) && ($feature->end+$end > $length{$feature->seq_id})){
					print "2";
					my $file = $fna_output_directory.$title."-".$type."_after-genfam.bed";
					unless(open FILE, '>>'.$file) {
						die "Unable to create $file";
					}
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}
					print FILE join("\t",$feature->seq_id,$feature->end+$begin,$length{$feature->seq_id},$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
					close FILE;
				}
				elsif(($gene{$j-1}->end > $feature->end+$begin) && ($gene{$j-1}->end > $feature->end)){			
					print "3";
					my $file = $fna_output_directory.$title."-".$type."_after-genfam.bed";
					unless(open FILE, '>>'.$file) {
						die "Unable to create $file";
					}
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->end+$begin,$feature->end+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
					close FILE;
				}
				elsif(($gene{$j+1}->start <= $feature->end+$end) && ($gene{$j-1}->end <= $feature->end)){			
					print "4";
					my $file = $fna_output_directory.$title."-".$type."_after-genfam.bed";
					unless(open FILE, '>>'.$file) {
						die "Unable to create $file";
					}
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->end+$begin,$gene{$j+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
					close FILE;
				}
				elsif(($gene{$j+1}->start > $feature->end+$end) && ($gene{$j-1}->end <= $feature->start)){			
					print "5";
					my $file = $fna_output_directory.$title."-".$type."_after-genfam.bed";
					unless(open FILE, '>>'.$file) {
						die "Unable to create $file";
					}
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->end+$begin,$feature->end+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
					close FILE;
				}
				elsif(($gene{$j+1}->start > $feature->end) && ($feature->end+$end >= $gene{$j+1}->start)){			
					print "6";
					my $file = $fna_output_directory.$title."-".$type."_after-genfam.bed";
					unless(open FILE, '>>'.$file) {
						die "Unable to create $file";
					}
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->end+$begin,$gene{$j+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
					close FILE;
				}
				else{
					print "7";
					my $file = $fna_output_directory.$title."-".$type."_after-genfam.bed";
					unless(open FILE, '>>'.$file) {
						die "Unable to create $file";
					}
					#change -1 to - and 1 to +
					if($feature->strand=~/-/){
						$strand1='-';
					}else{
						$strand1='+';
					}

					#gff2bed:
					print FILE join("\t",$feature->seq_id,$feature->end+$begin,$feature->end+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,".",$strand1)."\n";
					close FILE;
				}				
			}
		}
	}
	$gff1->close();
	
	system("bedtools getfasta -s -fi $genome -bed ".$fna_output_directory.$title."-".$type."_before-genfam.bed -fo $fna_output_file -name");
}

sub bedtools {

	my @tab_type = ("gene", "mRNA", "intron", "five_prime_UTR", "three_prime_UTR");
	
	#bedtools
	#usage system: system("bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>");
	if ($type =~/all/i){
		foreach my $VAR (@tab_type){
			my $fich = $fna_output_directory.$title."-".$VAR."-genfam.bed";
			if (-e $fich){
				system("bedtools getfasta -s -fi $genome -bed ".$fich." -fo ".$fna_output_directory.$title2."-".$VAR."-genfam.fna -name");
			}
		}			
	}
	elsif ($type =~/gene|polypeptide|mRNA|exon|intron|five_prime_UTR|three_prime_UTR|CDS/i){
		system("bedtools getfasta -s -fi $genome -bed ".$fna_output_directory.$title."-".$type."-genfam.bed -fo $fna_output_file -name");
	}
	
}

sub gff3_stat {
	my $file = $fna_output_directory.$title."-locus_tag-genfam.gff3";
	
	my $gff1 = new Bio::Tools::GFF(
		-file => $file,
		-gff_version => 3
	);
	
	my $gff = new Bio::Tools::GFF(
		-file => $file,
		-gff_version => 3
	);
	my (%stats, %mrna, %cds, %three, %five, %exon, %chr_length, $seq_obj);
	
	while(my $feature = $gff1->next_feature ){
		$chr_length{$feature->seq_id} = $feature->end;
	}
	$gff1->close; 
	
	my @ftypes = qw(gene transposable_element_gene mrna te cds exon five_prime_utr three_prime_utr intron intergenic genic);
	for my $t ( @ftypes ) {
		$stats{$t} = Statistics::Descriptive::Full->new;
	}
	my $id;
	my $mrna_id;
	my %gene;
	my %keep;	
	while(my $feature = $gff->next_feature) { 
		if ($feature->primary_tag() eq "mRNA") {
			($id) = $feature->get_tag_values("Parent");
			($mrna_id) = $feature->get_tag_values("ID");
			if (defined $gene{$id}){
				next;
			}
			else {
				$gene{$id} = 1;
				$keep{$mrna_id}  = 1;
				push @{$mrna{$feature->seq_id}{gene}},$feature;
			}
		}
		if ($feature->primary_tag() eq "transposable_element_gene") {	
			($id) = $feature->get_tag_values("ID");		
			if (defined $gene{$id}){
				next;
			}
			else {
				$gene{$id} = 1;
				push @{$mrna{$feature->seq_id}{transposable_element_gene}},$feature;
			}
		}	
		if ($feature->primary_tag() eq "CDS") {
			my ($parent) = $feature->get_tag_values("Parent");
			if (  $keep{$parent}){
				push @{$cds{$parent}},$feature;
			}
		}	
		if ($feature->primary_tag() eq "three_prime_UTR") {
			my ($parent) = $feature->get_tag_values("Parent");
			if (  $keep{$parent}){
				push @{$three{$parent}},$feature;
			}
		}	
		if ($feature->primary_tag() eq "five_prime_UTR") {
			my ($parent) = $feature->get_tag_values("Parent");
			if (  $keep{$parent}){
				push @{$five{$parent}},$feature;
			}
		}	
		if ($feature->primary_tag() eq "exon") {
			my ($parent) = $feature->get_tag_values("Parent");
			if (  $keep{$parent}){
				push @{$exon{$parent}},$feature;
			}
		}	
	}
	$gff->close;
	open(OUT,">".$fna_output_directory."intergenic.txt");
	my $genome = 0;
	my $mrna_sum = 0;
	my $te_sum  = 0;

	print join("\t","Chr","Length (bp)","Num Gene","Num TE"),"\n";
	foreach my $seq_id (sort {$a cmp $b} keys %chr_length){
		my $length 	= $chr_length{$seq_id};
		my @mrna 	= sort {$a->start <=> $b->start} @{$mrna{$seq_id}{gene}};
		
		my @mrna_te;
		if (defined $mrna{$seq_id}{transposable_element_gene}){
			@mrna_te = sort {$a->start <=> $b->start} @{$mrna{$seq_id}{transposable_element_gene}};
		}
		$genome	  += $length;
		$mrna_sum += scalar(@mrna);
		$te_sum   += scalar(@mrna_te);
		my @gene;
		push @gene , @mrna,@mrna_te; 
		print join("\t",$seq_id,$chr_length{$seq_id},scalar(@mrna),scalar(@mrna_te),scalar(@gene)),"\n";
		@gene = sort {$a->start <=> $b->start} @gene;
		my $total = 0; 
		my $lastfeature;
		my $cpt = 0;
		my $genic;
		my $intergenic;
		foreach my $feature (@gene) {
			$cpt++;
			my ($mrna_id) = $feature->get_tag_values("Name");
			if (defined $cds{$mrna_id}){
				foreach my $feature_cds (sort{$a->start <=> $b->start} @{$cds{$mrna_id}}) {
					my $cds_length = abs($feature_cds->start - $feature_cds->end);
					$stats{'cds'}->add_data($cds_length); 
					$total++;
				}
			}
			if (defined $five{$mrna_id}){
				foreach my $feature_five (sort{$a->start <=> $b->start} @{$five{$mrna_id}}) {
					my $length = abs($feature_five->start - $feature_five->end);
					$stats{'five_prime_utr'}->add_data($length); 
					$total++;
				}
			}
			if (defined $three{$mrna_id}){
				foreach my $feature_three (sort{$a->start <=> $b->start} @{$three{$mrna_id}}) {
					my $length = abs($feature_three->start - $feature_three->end);
					$stats{'three_prime_utr'}->add_data($length); 
					$total++;
				}
			}
			my $lastexon;
			
			if (defined $exon{$mrna_id}){
				for my $exon ( sort { $a->start  <=>   $b->start  }  @{$exon{$mrna_id}} ) {
					my $exonlen = abs($exon->start - $exon->end);
					$stats{'exon'}->add_data($exonlen);
					if( $lastexon ) {
						my $intronlen = abs($exon->start - ($lastexon->end)); 
						$stats{'intron'}->add_data($intronlen);
					}
					$lastexon = $exon;
				}
			}
			if ($lastfeature) {
				$intergenic = $feature->start - $lastfeature->end;
				if ($intergenic > 0) {
					$stats{'intergenic'}->add_data($intergenic);  
					print OUT join("\t",$cpt,"Inter",$seq_id,$lastfeature->end,$feature->start,$intergenic),"\n";  
				} 
				$genic = abs($feature->start - $feature->end);
				if ($feature->primary_tag() eq "mRNA") {
					$stats{'gene'}->add_data($genic); 
				}
				else {
					$stats{'te'}->add_data($genic); 
				}
				$stats{'genic'}->add_data($genic);  
				print OUT join("\t",$cpt,$feature->primary_tag,$seq_id,$feature->start,$feature->end,$genic),"\n"; 
				$lastfeature = $feature;
			}
			else {
				$intergenic = $feature->start;
				$stats{'intergenic'}->add_data($intergenic);   
				$genic = abs($feature->start - $feature->end);
				print OUT join("\t",$cpt,$feature->primary_tag,$seq_id,$feature->start,$feature->end,$genic),"\n";
				$stats{'genic'}->add_data($genic);  
				if ($feature->primary_tag() eq "mRNA") {
					$stats{'gene'}->add_data($genic); 
				}
				else {
					$stats{'te'}->add_data($genic); 
				}
				$lastfeature = $feature;
			}
		}
		$intergenic = abs($length - $gene[-1]->end ) ;
		$stats{'intergenic'}->add_data($intergenic);  
		print OUT join("\t",$cpt,"Intergenic",$seq_id,$gene[-1]->end +1,$length, $mrna_id),"\n\n";
 
	}
	print join("\t","All",$genome,$mrna_sum,$te_sum),"\n";
	for my $t ( qw (gene te intron exon cds five_prime_utr three_prime_utr intergenic genic ) ) {
		my $percent = 100 * $stats{$t}->sum / $genome;
		print join("\t",$t, $stats{$t}->count,$stats{$t}->sum,   $stats{$t}->mean,$percent),"\n";
   
	}
}

# Script options
#################

=pod

=head1 OPTIONS

#--- describes parameters given to the script
#+++ command line syntax
#--- requirement of the option and its parameter can be:
#--- required: name or nature inside <>
#--- optional: name or nature inside []
#--- alternative between 2 elements: elements separated by a |

=head2 Parameters

=over 4

=item B<[option_name]> ([option nature]): #+++

[option description]. #+++
Default: [option default value if one] #+++

=back
#--- Example:
#---
#--- Template.pl [-help | -man]
#---
#--- Template.pl [-debug [debug_level]] [-size <width> [height]]
#---
#--- =over 4
#---
#--- =item B<-help>:
#---
#--- Prints a brief help message and exits.
#---
#--- =item B<-man>:
#---
#--- Prints the manual page and exits.
#---
#--- =item B<-debug> (integer):
#---
#--- Executes the script in debug mode. If an integer value is specified, it will
#--- be the debug level. If "-debug" option was used without specifying a debug
#--- level, level 1 is assumed.
#--- Default: 0 (not in debug mode).
#---
#---=item B<-size> (positive_real) (positive_real):
#---
#--- Set the dimensions of the object that will be drawn. The first value is
#--- the width; the height is the second value if specified, otherwise it will
#--- assume height and width are equal to the first value.
#--- Default: width and height are set to 1.
#---
#---=back

=cut

# CODE START
#############

# options processing
my ($man, $help, $debug);

# parse options and print usage if there is a syntax error.
GetOptions("help|?"   => \$help,
           "man"      => \$man,
           "debug"    => \$debug, # a flag
           "bed=s"      => \$bed,
           "type|t=s"   => \$type,
           "begin=i" => \$begin,
           "end=i" => \$end)
#           "length=i" => \$length, # numeric
#           "file=s"   => \$data) # a string
    or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}


print "Looking for files in $gff3_file\n";
my @files = $gff3_file;
print "Found " . @files . " files to process...\n";

#sortGff3();

if (defined($bed)){
	sortGff3();
	gff3_stat();
	gff3tobed();
}

if (defined($type)){
	if ((defined($begin)) && ($begin == 0) && ($end == 0)){
		bedtools();
	}
	elsif ((defined($begin)) && ($begin =~ /^\-[1-9]+/) && ($end =~ /^\-[1-9]+|^0/)){
		gff3tobed_before();
	}
	elsif ((defined($begin)) && ($begin =~ /^[1-9]+|^0/) && ($end =~ /^[1-9]+/)){
		gff3tobed_after();
	}
	else{
		bedtools();
	}
}

# CODE END
###########


=pod

=head1 DIAGNOSTICS

#--- Give and explain here every error message the the script may generate.
#--- Use: (W)=warning (optional), (D)=deprecation (optional),
#---      (S)=severe warning (mandatory), (F)=fatal error (trappable),
#---      (X)=very fatal error (non-trappable).
#--- example:
#---
#--- =over 4
#---
#--- =item *
#---
#--- Negative radius;
#--- (F) Can not draw a circle with a negative radius.
#---
#---=back #+++

=head1 AUTHORS

Jonathan LORENZO (CIRAD), jonathan.lorenzo@cirad.fr

[author1_name (company), email]#+++

[author2_name (company), email]#+++

=head1 VERSION

Version [version.subversion.build] #+++

Date [DD/MM/YYY] #+++

=head1 SEE ALSO

#--- other documentation or objects related to this package
[perl(1), L<Geometry::Square>] #+++

=cut
