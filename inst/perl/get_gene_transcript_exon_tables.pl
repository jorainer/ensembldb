#!/usr/bin/perl
#####################################
## version 0.3.4: * Add columns gene_id_version and tx_id_version to the gene
##                  and transcript tables.
## version 0.3.3: * Write the species' scientific name to the Organism metadata
##                  field.
## version 0.3.1: * Add ens_counts.txt with the total counts of genes, tx, exons
##                  and proteins for validation that all entries were added to
##                  the database.
##                * Extract gene descriptions and tx support level.
## version 0.3.0: * Change database layout by adding a dedicated entrezgene
##                  table.
## version 0.2.4: * Extract taxonomy ID and add that to  metadata table.
## version 0.2.3: * Add additional columns to the uniprot table:
##                  o uniprot_db: the Uniprot database name.
##                  o uniprot_mapping_type: method by which the Uniprot ID was
##                    mapped to the Ensembl protein ID.
## version 0.2.2: * Transform gene coordinates always to toplevel instead of
##                  try-and-error transformation to chromosome.
## version 0.2.1: * Get protein IDs and (eventually) Uniprot IDs.
## version 0.0.2: * get also gene_seq_start, gene_seq_end, tx_seq_start and tx_seq_end from the database!
##                * did rename chrom_start to seq_start.

## uses environment variable ENS pointing to the
## ENSEMBL API on the computer
use lib $ENV{ENS} || $ENV{PERL5LIB};
use IO::File;
use Getopt::Std;
use strict;
use warnings;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::Registry;
## unification function for arrays
use List::MoreUtils qw/ uniq /;
my $script_version = "0.3.4";
my $min_tsl_version = 87;   ## The minimal required Ensembl version providing support for the tsl method.

## connecting to the ENSEMBL data base
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
my $user = "anonymous";
my $host = "ensembldb.ensembl.org";
my $port = 5306;
my $pass = "";
my $registry = 'Bio::EnsEMBL::Registry';
my $ensembl_version="none";
my $ensembl_database="core";
my $species = "human";
my $slice;
my $coord_system_version="unknown";
## get all gene ids defined in the database...
my @gene_ids = ();

my %option=();
getopts("e:hH:P:p:U:s:",\%option);
if($option{ h }){
  ## print help and exit.
  print("\nget_gene_transcript_exon_tables version ".$script_version.".\n");
  print("Retrieves gene/transcript/exon annotations from Ensembl and stores it as tabulator delimited text files.\n\n");
  print("usage: perl get_gene_transcript_exon_tables -e:hH:P:U:s:\n");
  print("-e (required): the Ensembl version (e.g. -e 75). The function will internally check if the submitted version matches the Ensembl API version and database queried.\n");
  print("-H (optional): the hostname of the Ensembl database; defaults to ensembldb.ensembl.org.\n");
  print("-h (optional): print this help message.\n");
  print("-p (optional): the port to access the Ensembl database.\n");
  print("-P (optional): the password to access the Ensembl database.\n");
  print("-U (optional): the username to access the Ensembl database.\n");
  print("-s (optional): the species; defaults to human.\n");
  print("\n\nThe script will generate the following tables:\n");
  print("- ens_gene.txt: contains all genes defined in Ensembl.\n");
  print("- ens_entrezgene.txt: contains mapping between ensembl gene_id and entrezgene ID.\n");
  print("- ens_transcript.txt: contains all transcripts of all genes.\n");
  print("- ens_exon.txt: contains all (unique) exons, along with their genomic alignment.\n");
  print("- ens_tx2exon.txt: relates transcript ids to exon ids (m:n), along with the index of the exon in the respective transcript (since the same exon can be part of different transcripts and have a different index in each transcript).\n");
  print("- ens_chromosome.txt: the information of all chromosomes (chromosome/sequence/contig names). \n");
  print("- ens_protein.txt: the mapping between (protein coding) transcripts and protein IDs including also the peptide sequence.\n");
  print("- ens_protein_domain.txt: provides for each protein all annotated protein domains along with their start and end coordinates on the protein sequence.");
  print("- ens_uniprot.txt: provides the mapping between Ensembl protein IDs and Uniprot IDs (if available). The mapping can be 1:n.");
  print("- ens_metadata.txt\n");
  exit 0;
}

if(defined($option{ s })){
	$species=$option{ s };
}
if(defined($option{ U })){
	$user=$option{ U };
}
if(defined($option{ H })){
	$host=$option{ H };
}
if(defined($option{ P })){
	$pass=$option{ P };
}
if(defined($option{ p })){
  $port=$option{ p };
}
if(defined($option{ e })){
	$ensembl_version=$option{ e };
}else{
	die("The ensembl version has to be specified with the -e parameter (e.g. -e 75)");
}

my $api_version="".software_version()."";
if($ensembl_version ne $api_version){
    die "The submitted Ensembl version (".$ensembl_version.") does not match the version of the Ensembl API (".$api_version."). Please configure the environment variable ENS to point to the correct API.";
}
my $ensembl_version_num = $ensembl_version + 0;

print "Connecting to ".$host." at port ".$port."\n";

# $registry->load_registry_from_db(-host => $host, -user => $user,
# 				 -pass => $pass, -port => $port,
# 				 -verbose => "1");
$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -port => $port);
my $gene_adaptor = $registry->get_adaptor($species, $ensembl_database, "gene");
my $slice_adaptor = $registry->get_adaptor($species, $ensembl_database, "slice");
my $meta_container = $registry->get_adaptor($species, $ensembl_database,
					    'MetaContainer' );
## determine the species:
my $species_id = $gene_adaptor->db->species_id;
my $species_ens = $gene_adaptor->db->species;
## Determine the taxonomy ID:
my $taxonomy_id = 0;
$taxonomy_id = $meta_container->get_taxonomy_id();
my $common_name = $meta_container->get_common_name();
my $scientific_name = $meta_container->get_scientific_name();

my $infostring = "# get_gene_transcript_exon_tables.pl version $script_version:\nRetrieve gene models for Ensembl version $ensembl_version, species $species from Ensembl database at host: $host\n";

print $infostring;

## preparing output files:
open(GENE , ">ens_gene.txt");
print GENE "gene_id\tgene_name\tgene_biotype\tgene_seq_start\tgene_seq_end\tseq_name\tseq_strand\tseq_coord_system\tdescription\tgene_id_version\n";

open(TRANSCRIPT , ">ens_tx.txt");
print TRANSCRIPT "tx_id\ttx_biotype\ttx_seq_start\ttx_seq_end\ttx_cds_seq_start\ttx_cds_seq_end\tgene_id\ttx_support_level\ttx_id_version\n";

open(EXON , ">ens_exon.txt");
print EXON "exon_id\texon_seq_start\texon_seq_end\n";

open(ENTREZGENE, ">ens_entrezgene.txt");
print ENTREZGENE "gene_id\tentrezid\n";
# open(G2T , ">ens_gene2transcript.txt");
# print G2T "g2t_gene_id\tg2t_tx_id\n";

open(T2E , ">ens_tx2exon.txt");
print T2E "tx_id\texon_id\texon_idx\n";

open(PROTEIN, ">ens_protein.txt");
## print PROTEIN "tx_id\tprotein_id\tuniprot_id\tprotein_sequence\n";
print PROTEIN "tx_id\tprotein_id\tprotein_sequence\n";

open(UNIPROT, ">ens_uniprot.txt");
print UNIPROT "protein_id\tuniprot_id\tuniprot_db\tuniprot_mapping_type\n";

open(PROTDOM, ">ens_protein_domain.txt");
print PROTDOM "protein_id\tprotein_domain_id\tprotein_domain_source\tinterpro_accession\tprot_dom_start\tprot_dom_end\n";

open(CHR , ">ens_chromosome.txt");
print CHR "seq_name\tseq_length\tis_circular\n";

open(COUNTS, ">ens_counts.txt");
print COUNTS "gene\ttx\texon\tprotein\n";

##OK now running the stuff:
print "Start fetching data\n";
my %done_chromosomes=();
my %done_exons=();  ## to keep track of which exons have already been saved.
my $counta = 0;
my $count_gene = 0;
my $count_tx = 0;
my $count_exon = 0;
my $count_protein = 0;
@gene_ids = @{$gene_adaptor->list_stable_ids()};
foreach my $gene_id (@gene_ids){
  $counta++;
  $count_gene++;
  if(($counta % 2000) == 0){
    print "processed $counta genes\n";
  }
  my $orig_gene;

  $orig_gene = $gene_adaptor->fetch_by_stable_id($gene_id);
  if(defined $orig_gene){
    my $do_transform=1;
    ## Instead of transforming to chromosome we transform to 'toplevel',
    ## for genes encoded on chromosome this should be the chromosome, for others
    ## the most "top" level sequence.
    ## my $gene  = $orig_gene->transform("chromosome");
    my $gene  = $orig_gene->transform("toplevel");
    if(!defined $gene){
      ## gene is not on known defined chromosomes!
      $gene = $orig_gene;
      $do_transform=0;
    }
    my $coord_system = $gene->coord_system_name;
    my $chrom = $gene->slice->seq_region_name;
    my $strand = $gene->strand;

    ## check if we did already fetch some info for that chromosome
    if(exists($done_chromosomes{ $chrom })){
      ## don't do anything...
    }else{
      $done_chromosomes{ $chrom } = "done";
      my $chr_slice = $gene->slice->seq_region_Slice();
      my $name = $chr_slice->seq_region_name;
      my $length = $chr_slice->length;
      my $is_circular = $chr_slice->is_circular;
      print CHR "$name\t$length\t$is_circular\n";
      my $tmp_version = $chr_slice->coord_system()->version();
      if (defined $tmp_version and length $tmp_version) {
	$coord_system_version = $tmp_version;
      }
      # my $chr_slice_again = $slice_adaptor->fetch_by_region('chromosome', $chrom);
      # if(defined($chr_slice_again)){
      # 	$coord_system_version = $chr_slice_again->coord_system()->version();
      # }
    }

    ## get information for the gene.
    my $gene_external_name= $gene->external_name;
    if(!defined($gene_external_name)){
      $gene_external_name="";
    }
    my $gene_id_version = $gene->stable_id_version;
    my $gene_biotype = $gene->biotype;
    my $gene_seq_start = $gene->start;
    my $gene_seq_end = $gene->end;
    my $description = $gene->description;
    if(!defined($description)){
      $description = "NULL";
    }
    ## get entrezgene ID, if any...
    my $all_entries = $gene->get_all_DBLinks("EntrezGene");
    foreach my $dbe (@{$all_entries}){
      print ENTREZGENE "$gene_id\t".$dbe->primary_id."\n";
    }
    print GENE "$gene_id\t$gene_external_name\t$gene_biotype\t$gene_seq_start\t$gene_seq_end\t$chrom\t$strand\t$coord_system\t$description\t$gene_id_version\n";

    ## process transcript(s)
    my @transcripts = @{ $gene->get_all_Transcripts };
    ## ok looping through the transcripts
    foreach my $transcript (@transcripts){
      $count_tx++;
      if($do_transform==1){
	## just to be shure that we have the transcript in chromosomal coordinations.
	## $transcript = $transcript->transform("chromosome");
	$transcript = $transcript->transform("toplevel");
      }
      ##my $tx_start = $transcript->start;
      ##my $tx_end = $transcript->end;

      ## caution!!! will get undef if transcript is non-coding!
      my $tx_cds_start = $transcript->coding_region_start;
      if(!defined($tx_cds_start)){
	$tx_cds_start = "NULL";
      }
      my $tx_cds_end = $transcript->coding_region_end;
      if(!defined($tx_cds_end)){
	$tx_cds_end = "NULL";
      }
      my $tx_id = $transcript->stable_id;
      my $tx_id_version = $transcript->stable_id_version;
      my $tx_biotype = $transcript->biotype;
      my $tx_seq_start = $transcript->start;
      my $tx_seq_end = $transcript->end;
      my $tx_tsl = "NULL";
      if ($ensembl_version_num >= $min_tsl_version) {
	$tx_tsl = $transcript->tsl;
	if (!defined($tx_tsl)) {
	  $tx_tsl = "NULL";
	}
      }
      my $tx_description = $transcript->description;
      # if (!defined($tx_description)) {
      # 	$tx_description = "NULL";
      # }
      ## write info.
      print TRANSCRIPT "$tx_id\t$tx_biotype\t$tx_seq_start\t$tx_seq_end\t$tx_cds_start\t$tx_cds_end\t$gene_id\t$tx_tsl\t$tx_id_version\n";
##      print G2T "$gene_id\t$tx_id\n";

      ## Process proteins/translations (if possible)
      my $transl = $transcript->translation();
      if (defined($transl)) {
	$count_protein++;
	my $transl_id = $transl->stable_id();
	my $prot_seq = $transl->seq();
	## Check if we could get UNIPROT ID(s):
	my @unip = @{ $transl->get_all_DBLinks('Uniprot%') };
	if (scalar(@unip) > 0) {
	  foreach my $uniprot (@unip) {
	    my $unip_id = $uniprot->display_id();
	    # my $unip_acc = $uniprot->primary_id();
	    my $dbn = $uniprot->dbname();
	    $dbn =~ s/Uniprot\///g;
	    my $maptype = $uniprot->info_type();
	    print UNIPROT "$transl_id\t$unip_id\t$dbn\t$maptype\n";
	    ## print PROTEIN "$tx_id\t$transl_id\t$unip_id\t$prot_seq\n";
	  }
	}
	print PROTEIN "$tx_id\t$transl_id\t$prot_seq\n";
	my $prot_doms = $transl->get_all_DomainFeatures;
	while ( my $prot_dom = shift @{$prot_doms}) {
	  my $logic_name = $prot_dom->analysis()->logic_name();
	  my $prot_dom_id = $prot_dom->display_id();
	  my $interpro_acc = $prot_dom->interpro_ac();
	  my $prot_start = $prot_dom->start();
	  my $prot_end = $prot_dom->end();
	  print PROTDOM "$transl_id\t$prot_dom_id\t$logic_name\t$interpro_acc\t$prot_start\t$prot_end\n";
	}
      }
      ## process exon(s)
      ##my @exons = @{ $transcript->get_all_Exons(-constitutive => 1) };
      my @exons = @{ $transcript->get_all_Exons() };  ## exons always returned 5' 3' of transcript!
      my $current_exon_idx = 1;
      foreach my $exon (@exons){
	if($do_transform==1){
	  ## $exon->transform("chromosome");
	  $exon->transform("toplevel");
	}
	my $exon_start = $exon->start;
	my $exon_end = $exon->end;
	my $exon_id = $exon->stable_id;

	## write info, but only if we didn't already saved this exon (exon can be
	## part of more than one transcript).
	if(exists($done_exons{ $exon_id })){
	  ## don't do anything.
	}else{
	  $done_exons{ $exon_id } = 1;
	  $count_exon++;
	  print EXON "$exon_id\t$exon_start\t$exon_end\n";
	}
	## saving the exon id to this file that provides the n:m mappint; also saving
	## the index of the exon in the present transcript to that.
	print T2E "$tx_id\t$exon_id\t$current_exon_idx\n";

	$current_exon_idx++;
      }
    }
  }
}

## want to save:
## data, ensembl host, species, ensembl version, genome build?
open(INFO , ">ens_metadata.txt");
print INFO "name\tvalue\n";
print INFO "Db type\tEnsDb\n";
print INFO "Type of Gene ID\tEnsembl Gene ID\n";
print INFO "Supporting package\tensembldb\n";
print INFO "Db created by\tensembldb package from Bioconductor\n";
print INFO "script_version\t$script_version\n";
print INFO "Creation time\t".localtime()."\n";
print INFO "ensembl_version\t$ensembl_version\n";
print INFO "ensembl_host\t$host\n";
print INFO "Organism\t$scientific_name\n";
print INFO "taxonomy_id\t$taxonomy_id\n";
print INFO "genome_build\t$coord_system_version\n";
print INFO "DBSCHEMAVERSION\t2.1\n";

print COUNTS "$count_gene\t$count_tx\t$count_exon\t$count_protein\n";

close(INFO);

close(GENE);
close(TRANSCRIPT);
close(EXON);
close(ENTREZGENE);
##close(G2T);
close(T2E);
close(CHR);
close(PROTEIN);
close(PROTDOM);
close(UNIPROT);
close(COUNTS);
