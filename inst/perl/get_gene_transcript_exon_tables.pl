#!/usr/bin/perl
#####################################
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
my $script_version = "0.1.3";

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
  print("- ens_transcript.txt: contains all transcripts of all genes.\n");
  print("- ens_exon.txt: contains all (unique) exons, along with their genomic alignment.\n");
  print("- ens_tx2exon.txt: relates transcript ids to exon ids (m:n), along with the index of the exon in the respective transcript (since the same exon can be part of different transcripts and have a different index in each transcript).\n");
  print("- ens_chromosome.txt: the information of all chromosomes (chromosome/sequence/contig names). \n");
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

print "Connecting to ".$host." at port ".$port."\n";

# $registry->load_registry_from_db(-host => $host, -user => $user,
# 				 -pass => $pass, -port => $port,
# 				 -verbose => "1");
$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -port => $port);
my $gene_adaptor = $registry->get_adaptor($species, $ensembl_database, "gene");
my $slice_adaptor = $registry->get_adaptor($species, $ensembl_database, "slice");

## determine the species:
my $species_id = $gene_adaptor->db->species_id;
my $species_ens = $gene_adaptor->db->species;

my $infostring = "# get_gene_transcript_exon_tables.pl version $script_version:\nRetrieve gene models for Ensembl version $ensembl_version, species $species from Ensembl database at host: $host\n";

print $infostring;

## preparing output files:
open(GENE , ">ens_gene.txt");
print GENE "gene_id\tgene_name\tentrezid\tgene_biotype\tgene_seq_start\tgene_seq_end\tseq_name\tseq_strand\tseq_coord_system\n";

open(TRANSCRIPT , ">ens_tx.txt");
print TRANSCRIPT "tx_id\ttx_biotype\ttx_seq_start\ttx_seq_end\ttx_cds_seq_start\ttx_cds_seq_end\tgene_id\n";

open(EXON , ">ens_exon.txt");
print EXON "exon_id\texon_seq_start\texon_seq_end\n";

# open(G2T , ">ens_gene2transcript.txt");
# print G2T "g2t_gene_id\tg2t_tx_id\n";

open(T2E , ">ens_tx2exon.txt");
print T2E "tx_id\texon_id\texon_idx\n";

open(CHR , ">ens_chromosome.txt");
print CHR "seq_name\tseq_length\tis_circular\n";

##OK now running the stuff:
print "Start fetching data\n";
my %done_chromosomes=();
my %done_exons=();  ## to keep track of which exons have already been saved.
my $counta = 0;
@gene_ids = @{$gene_adaptor->list_stable_ids()};
foreach my $gene_id (@gene_ids){
  $counta++;
  if(($counta % 2000) == 0){
    print "processed $counta genes\n";
  }
  my $orig_gene;

  $orig_gene = $gene_adaptor->fetch_by_stable_id($gene_id);
  if(defined $orig_gene){
    my $do_transform=1;
    my $gene  = $orig_gene->transform("chromosome");
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
      my $chr_slice_again = $slice_adaptor->fetch_by_region('chromosome', $chrom);
      if(defined($chr_slice_again)){
	$coord_system_version = $chr_slice_again->coord_system()->version();
      }
      # if(defined $chr_slice){
      # 	my $name = $chr_slice->seq_region_name;
      # 	my $length = $chr_slice->length;
      # 	my $is_circular = $chr_slice->is_circular;
      # 	$coord_system_version = $chr_slice->coord_system()->version();
      # 	print CHR "$name\t$length\t$is_circular\n";
      # }else{
      # 	my $length = $gene->slice->seq_region_length();
      # 	print CHR "$chrom\t0\t0\n";
      # }
    }

    ## get information for the gene.
    my $gene_external_name= $gene->external_name;
    if(!defined($gene_external_name)){
      $gene_external_name="";
    }
    my $gene_biotype = $gene->biotype;
    my $gene_seq_start = $gene->start;
    my $gene_seq_end = $gene->end;
    ## get entrezgene ID, if any...
    my $all_entries = $gene->get_all_DBLinks("EntrezGene");
    my %entrezgene_hash=();
    foreach my $dbe (@{$all_entries}){
      $entrezgene_hash{ $dbe->primary_id } = 1;
    }
    my $hash_size = keys %entrezgene_hash;
    my $entrezid = "";
    if($hash_size > 0){
      $entrezid = join(";", keys %entrezgene_hash);
    }
    print GENE "$gene_id\t$gene_external_name\t$entrezid\t$gene_biotype\t$gene_seq_start\t$gene_seq_end\t$chrom\t$strand\t$coord_system\n";

    ## process transcript(s)
    my @transcripts = @{ $gene->get_all_Transcripts };
    ## ok looping through the transcripts
    foreach my $transcript (@transcripts){
      if($do_transform==1){
	## just to be shure that we have the transcript in chromosomal coordinations.
	$transcript = $transcript->transform("chromosome");
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
      my $tx_biotype = $transcript->biotype;
      my $tx_seq_start = $transcript->start;
      my $tx_seq_end = $transcript->end;
      ## write info.
      print TRANSCRIPT "$tx_id\t$tx_biotype\t$tx_seq_start\t$tx_seq_end\t$tx_cds_start\t$tx_cds_end\t$gene_id\n";
##      print G2T "$gene_id\t$tx_id\n";

      ## process exon(s)
      ##my @exons = @{ $transcript->get_all_Exons(-constitutive => 1) };
      my @exons = @{ $transcript->get_all_Exons() };  ## exons always returned 5' 3' of transcript!
      my $current_exon_idx = 1;
      foreach my $exon (@exons){
	if($do_transform==1){
	  $exon->transform("chromosome");
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
print INFO "Organism\t$species_ens\n";
print INFO "genome_build\t$coord_system_version\n";
print INFO "DBSCHEMAVERSION\t1.0\n";

close(INFO);

close(GENE);
close(TRANSCRIPT);
close(EXON);
##close(G2T);
close(T2E);
close(CHR);


