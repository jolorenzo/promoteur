#!/usr/bin/perl

=pod

=head1 NAME

[gff3fna2gff32bed2fasta.pl - short description]

=head1 SYNOPSIS

    qsub -q bioinfo.q -b yes -V -N format1 perl gff3fna2gff32bed2fasta.pl gff3_input_file genome_sequence_directory -bed -locus -type=gene|mRNA|polypeptide|... -begin=0|-x|+x -end=0|-x|+x

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
use List::MoreUtils qw(uniq);
use Bio::SeqIO;
use Bio::DB::Fasta;

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

my $title = $1 if $gff3_file =~ /[\\\/]?(.+)\.gff3/i;
my $title2 = $1 if $gff3_file =~ /[\\\/]?([A-Z]{5}.+)-sequence_feature-genfam\.gff3/i;


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
my ($type, $begin, $end);
sub sortGff3{

	my $gff = new Bio::Tools::GFF(
	-file => $gff3_file,
	-gff_version => 3
	);
	
	my $cpt = 0;
	my ($id, %feature, $mrna_id, %gene_refseq, %mrna, %polypeptide, %other, %cds, %seq, $seq_gene, $seq_cds, $product, $seq_pep);
	
	while(my $feature = $gff->next_feature) {
		if ($feature->primary_tag() eq "gene") {
			($id) = $feature->get_tag_values("ID");	
			$feature{$id} = $feature;
			$cpt++;
			push @{$gene_refseq{$feature->seq_id}} , $feature;
		}
		elsif ($feature->primary_tag eq "mRNA") {
			($mrna_id) = $feature->get_tag_values("ID"); 
			$mrna{$id}{$mrna_id} = $feature;
			#push @{$mrna_refseq{$feature->seq_id}} , $feature;
		}
		elsif ($feature->primary_tag eq "polypeptide") { 
			push @{$polypeptide{$mrna_id}}, $feature;
			my ($id) = $feature->get_tag_values("Derives_from"); 
		}
		elsif ($feature->primary_tag !~/polypeptide | gene | mRNA/i) { 	
			push @{$other{$mrna_id}}, $feature;
		}
	}
	$gff->close;
	
	my $outfile = "$title"."-locus_tag".".gff3";
	print "Print data to $outfile...\n";
	
	my $out = new Bio::Tools::GFF(
	-file => ">$outfile",
	-gff_version => 3
	);
	
	my %h_id = ( "ARATH" => "At", "BRADI" => "Bd", "GLYMA" => "Gm", "GOSRA"   => "Gr", "LOTJA"   => "Lj", "MEDTR"   => "Mt", "MUSAC"   => "Ma",
	"ORYSI"   => "Os", "ORYSJ"   => "Os", "POPTR"   => "Pt", "SOLLC"   => "Sl", "SORBI"   => "Sb", "THECC"   => "Tc", "VITVI"   => "Vv", 
	"MAIZE"   => "Zm", "MALDO"   => "Md", "MANES"   => "Me", "RICCO"   => "Rc", "SETIT"   => "Si", "SOLTU"   => "St");
	
	foreach my $seq_id (sort {$a cmp $b} keys%gene_refseq){
		#my $seqobj = $seq{$seq_id}; 
		my $count = 0;
		foreach my $gene (sort {$a->{start} <=> $b->{start}} @{$gene_refseq{$seq_id}}) {
			my ($name) = $gene->get_tag_values("Name");			 	
			#my $status = $gene->{status};
			#if (exists $feature{$name} ){
				#my $gene = $feature{$name}; 	
				$count++;
				my $gene_id = sprintf( "_g%05d", $count * 10 );
				my $first_poly_id = sprintf( "_p%05d", $count * 10 );
				my $first_mrna_id = sprintf( "_t%05d", $count * 10 );
				#my $seqobj_gene = $seqobj->trunc($gene->start,$gene->end);
				#$seqobj_gene->display_id($gene_id);
				#$seq_gene->write_seq($seqobj_gene); 
				#$gene->remove_tag("ID");
				#$gene->remove_tag("Name");
				#$gene->add_tag_value("ID",$gene_id);
				#$gene->add_tag_value("Name",$gene_id);
				#$gene->source_tag($source_tag);
				#$gene->add_tag_value("old_locus_tag",$name);
				
				my ($L5, $chr);
				if ($gene->seq_id =~ /^([A-Z]{5})\d+/i) {
					$L5 = $1;
				} elsif($gene->seq_id =~ /^([A-Z]{5})_scaffold\d+/i){
					$L5 = $1;
				}
				
				if ($gene->seq_id =~ /^[A-Z]{5}(\d+)/i){
					$chr = $1;
				}elsif($gene->seq_id =~ /^[A-Z]{5}_scaffold(\d+)/i){
					$chr = $1;
				}
				
				my $L2 = $h_id{$L5};
				
				$gene->add_tag_value("locus_tag", "$L2"."$chr"."$gene_id"."_$L5");
				$out->write_feature($gene);
				
				my $count_mrna = 0;
				foreach my $mrna (keys %{$mrna{$name}}){
					$count_mrna++;
					my $mrna_id = $first_mrna_id ."." . $count_mrna;
					my $poly_id = $first_poly_id ."." . $count_mrna;
					my $feat_mrna = $mrna{$name}{$mrna};
					my ($old_mrna_id) =  $feat_mrna->get_tag_values("ID");
				
					#$feat_mrna->remove_tag("ID");
					#$feat_mrna->source_tag($source_tag);
					#$feat_mrna->remove_tag("Name");
				
					#$feat_mrna->add_tag_value("old_locus_tag",$mrna);
					
					if ($feat_mrna->seq_id =~ /^([A-Z]{5})\d+/i) {
						$L5 = $1;
					} elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_scaffold\d+/i){
						$L5 = $1;
					}
				
					if ($feat_mrna->seq_id =~ /^[A-Z]{5}(\d+)/i){
						$chr = $1;
					}elsif($feat_mrna->seq_id =~ /^[A-Z]{5}_scaffold(\d+)/i){
						$chr = $1;
					}
					$feat_mrna->add_tag_value("locus_tag", "$L2"."$chr"."$mrna_id"."_$L5");
					#my $strand = $feat_mrna->strand;
					#$feat_mrna->remove_tag("Parent");
					#$feat_mrna->remove_tag("Note") if $feat_mrna->has_tag("Note");
					#$feat_mrna->add_tag_value("ID",$mrna_id);
					#$feat_mrna->add_tag_value("Name",$mrna_id);
					#$feat_mrna->add_tag_value("Note",$status);
					#$feat_mrna->add_tag_value("Parent",$gene_id);
					$out->write_feature($feat_mrna);
					
					if (exists $polypeptide{$old_mrna_id}){
						foreach my $poly (@{$polypeptide{$old_mrna_id}}){
							#$poly->remove_tag("ID");
							#$poly->source_tag($source_tag);
							#$poly->remove_tag("Name");
							#$poly->remove_tag("Derives_from");
							#$poly->add_tag_value("ID",$poly_id);
							#$poly->add_tag_value("Name",$poly_id);
							#$poly->add_tag_value("Derives_from",$mrna_id);
							#my ($note) = $poly->has_tag("note") ? $poly->get_tag_values("note") : join(" ~",$gene_id,"Hypothetical protein","unknown_gene","missing_functional_completeness");
							#my @data = (split(/~ /,$note));
							#my $first = shift @data;
							#unshift @data , $gene_id ;
							#my $new_note = join("~ ",@data); 
							#$poly->remove_tag("note") if $poly->has_tag("note");
							#$poly->add_tag_value("note",$new_note);
							
							if ($poly->seq_id =~ /^([A-Z]{5})\d+/i) {
								$L5 = $1;
							} elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold\d+/i){
								$L5 = $1;
							}
				
							if ($poly->seq_id =~ /^[A-Z]{5}(\d+)/i){
								$chr = $1;
							}elsif($poly->seq_id =~ /^[A-Z]{5}_scaffold(\d+)/i){
								$chr = $1;
							}
							$poly->add_tag_value("locus_tag", "$L2"."$chr"."$poly_id"."_$L5");
							$out->write_feature($poly);
						}
					}
					#else {	
					#	my $start_poly = 100000000000000000;
					#	my $end_poly   = -1; 
					#	my @cds;
					#	foreach my $other (@{$other{$old_mrna_id}}){		
					#		if ($other->primary_tag() eq "CDS") {
					#			$start_poly = $other->start if $other->start < $start_poly;
					#			$end_poly   = $other->end if $other->end > $end_poly;
					#			push @cds,$other;
					#		}	
					#	}	
					#	unless (@cds){
					#		#print $name ,"\t", $old_mrna_id,"\n";
					#	}
					#	else {
					#		if ($cds[0]->strand == 1) {
					#			@cds = sort{$a->start <=> $b->start} @cds;
					#		}
					#		else {
					#			@cds = sort{$b->start <=> $a->start} @cds;
					#		}
					#		my $cds;
					#		foreach my $feature (@cds) {	
					#			my $seqtrunc_cds = $seqobj->trunc($feature->start,$feature->end);	
					#			if ($feature->strand == -1) {
					#				$cds .= $seqtrunc_cds->revcom()->seq;
					#			}
					#			else {	
					#				$cds .= $seqtrunc_cds->seq;
					#			}
					#		}
					#		my $seqobj_cds = Bio::PrimarySeq->new (
					#			-display_id  => $mrna_id,
					#			-seq         => $cds ,
					#			-desc		=> $product
					#		);   
					#		if ($cds =~ /^ATG.*/ &&  $cds =~ /.*(TAG|TAA|TGA)$/){
					#			
					#			$is_complete++;
					#		}
					#		else {
					#			if ($cds =~ /^ATG.*/){
					#				$codon_stop++;
					#			}
					#			elsif ($cds =~ /.*(TAG|TAA|TGA)$/){
					#				$codon_start++;
					#			}
					#			else {
					#				$is_incomplete++;
					#				print $name ,"\t", $old_mrna_id,"\n";
					#			}
					#		}
					#		my $seqobj_pep = $seqobj_cds->translate();
					#		$seqobj_pep->display_id($poly_id);
					#		$seq_cds->write_seq($seqobj_cds);
					#		$seq_pep->write_seq($seqobj_pep);
					#		my ($product) = $feat_mrna->get_tag_values("Note");	
					#		my $note = join("~ ",$gene_id,$product,"unknown_gene","missing_functional_completeness");
					#		my $poly = new Bio::SeqFeature::Generic(
					#			-seq_id 	=> $seq_id,
					#			-source_tag => $source_tag,
					#			-primary_tag => 'polypeptide',
					#			-start       => $start_poly,
					#			-end         => $end_poly,
					#			-strand      => $strand,
					#			-tag 		 => {
					#				ID	=> $poly_id,
					#				Name	=> $poly_id,
					#				Derives_from => $mrna_id,
					#				inference => "{refseq}",
					#				owner  => "musa",
					#				alternative_splicing	=> "to_fill",
					#				annotator_comment =>"to_fill"	
					#			}
					#		); 
					#		if ($feat_mrna->has_tag("Dbxref")) {
					#			my @dbxref = $feat_mrna->get_tag_values("Dbxref");
					#			foreach my $dbxref (@dbxref) {
					#				$poly->add_tag_value("Dbxref",$dbxref);
					#			}
					#		}	
					#		$poly->add_tag_value("note",$note);
					#		$poly->add_tag_value("Ontology_term","PRODUCT:".$product);
					#		$poly->add_tag_value("Ontology_term" => "CC_functional_completeness:missing_functional_completeness");
					#		$poly->add_tag_value("Ontology_term" => "CC_evidence_code:ISS");
					#		$poly->add_tag_value("Ontology_term" => "CC_status:in_progress");
					#		$poly->add_tag_value("Ontology_term" => "CC_evidence:automatic");
					#		$poly->add_tag_value("Ontology_term" => "CC_EC_number:no_EC_number");
					#		$out->write_feature($poly);
					#	}
					#}	
					foreach my $other (@{$other{$old_mrna_id}}){		
						#$other->remove_tag("Parent");	
						#$other->source_tag($source_tag);
						#$other->add_tag_value("Parent",$mrna_id);
						$out->write_feature($other);	
					}	
				#}
			}					
		}
	}
	$out->close();
}

sub gff3tobed{

	my $new_file = "$title"."-locus_tag".".gff3";
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
	
	while(my $feature = $gffio->next_feature()) {
		if($feature->primary_tag =~/gene/i){
			
			my $file = "$title-gene.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			close FILE;
		}
	
		if($feature->primary_tag =~/mRNA/i){

			my $file = "$title-mRNA.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			close FILE;
		}
	
		if($feature->primary_tag =~/five_prime_UTR/i){

			my $file = "$title-five_prime_UTR.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start+1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{Parent}->[0],".",$strand1)."\n";
			close FILE;
		}
	
		if($feature->primary_tag =~/three_prime_UTR/i){

			my $file = "$title-three_prime_UTR.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start+1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{Parent}->[0],".",$strand1)."\n";
			close FILE;
		}
		
		if($feature->primary_tag =~/CDS/i){

			my $file = "$title-CDS.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start+1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{Parent}->[0],".",$strand1)."\n";
			close FILE;
		}
	
		if($feature->primary_tag =~/exon/i){

			my $file = "$title-exon.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start+1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{Parent}->[0],".",$strand1)."\n";
			close FILE;
		}
	
		if($feature->primary_tag =~/intron/i){

			my $file = "$title-intron.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start+1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{Parent}->[0],".",$strand1)."\n";
			close FILE;
		}
	
		if($feature->primary_tag =~/polypeptide/i){

			my $file = "$title-polypeptide.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start+1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1,$feature->end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			close FILE;
		}
	}
	$gffio->close();
  
}

sub gff3tobed_before{

	my $new_file = "$title"."-locus_tag".".gff3";
	my $gffio = Bio::Tools::GFF -> new(-file =>$new_file , -gff_version => 3);

	my $strand=0;
	my $strand1='-';
	
	##### Liste les fichiers ayant l'extension .bed 
	#my @list = glob("*.bed"); 
  
	### recupere le nombre de fichier 
	#my $numberoffile = scalar(@list);
	
	#for ( my $v = 0; $v < $numberoffile; $v++ ) {
	#	unlink $list[$v]; #supprime les fichiers
	#}
	
	while(my $feature = $gffio->next_feature()) {
		if($feature->primary_tag =~/$type/i){
			
			my $file = "$title-$type-before.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->start-$begin,$feature->start-$end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->start-1-$begin,$feature->start-1-$end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->start-1-$begin,$feature->start-1-$end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			close FILE;
		}
	}
	system("bedtools getfasta -s -fi $genome -bed $title-$type-before.bed -fo $title2-$type-before-genfam.fna");
}

sub gff3tobed_after{

	my $new_file = "$title"."-locus_tag".".gff3";
	my $gffio = Bio::Tools::GFF -> new(-file =>$new_file , -gff_version => 3);

	my $strand=0;
	my $strand1='-';
	
	##### Liste les fichiers ayant l'extension .bed 
	#my @list = glob("*.bed"); 
  
	### recupere le nombre de fichier 
	#my $numberoffile = scalar(@list);
	
	#for ( my $v = 0; $v < $numberoffile; $v++ ) {
	#	unlink $list[$v]; #supprime les fichiers
	#}
	
	while(my $feature = $gffio->next_feature()) {
		if($feature->primary_tag =~/$type/i){
			
			my $file = "$title-$type-after.bed";
			unless(open FILE, '>>'.$file) {
				die "Unable to create $file";
			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
				#print FILE join("\t",$feature->seq_id,$feature->end+1+$begin,$feature->end+1+$end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}else{
				$strand1='+';
				#print FILE join("\t",$feature->seq_id,$feature->end+$begin,$feature->end+$end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			}

			#gff2bed:
			print FILE join("\t",$feature->seq_id,$feature->end+$begin,$feature->end+$end,$feature->{_gsf_tag_hash}->{ID}->[0],".",$strand1)."\n";
			close FILE;
		}
	}
	system("bedtools getfasta -s -fi $genome -bed $title-$type-after.bed -fo $title2-$type-after-genfam.fna");
}

sub bedtools {

	#bedtools
	#usage system: system("bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>");
	system("bedtools getfasta -s -fi $genome -bed $title-$type.bed -fo $title2-$type-genfam.fna");

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
my ($man, $help, $debug, $bed);

# parse options and print usage if there is a syntax error.
GetOptions("help|?"   => \$help,
           "man"      => \$man,
           "debug"    => \$debug, # a flag
           "bed"      => \$bed,
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


if (defined($bed)){
	sortGff3();
	gff3tobed();
}

if (defined($type)){
	if (($begin == 0) && ($end == 0)){
		bedtools();
	}
	elsif (($begin =~ /\-[1-9]+/) && ($end =~ /\-[1-9]+/)){
		gff3tobed_before();
	}
	elsif (($begin =~ /\+[1-9]+/) && ($end =~ /\+[1-9]+/)){
		gff3tobed_after();
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
