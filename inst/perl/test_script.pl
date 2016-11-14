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
my $script_version = "0.2.2";

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

my $gene_id = "ENSG00000109906";

$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -port => $port);
my $gene_adaptor = $registry->get_adaptor($species, $ensembl_database, "gene");
my $slice_adaptor = $registry->get_adaptor($species, $ensembl_database, "slice");

my $current_gene = $gene_adaptor->fetch_by_stable_id($gene_id);
print "Current gene: ".$current_gene->display_id()."\n";
my @transcripts = @{ $current_gene->get_all_Transcripts };

foreach my $transcript (@transcripts){
  print "Current tx: ".$transcript->display_id()."\n";
  my $transl = $transcript->translation();
  if (defined($transl)) {
    my $transl_id = $transl->stable_id();
    print "Current translation ".$transl_id."\n";
    my $attr = $transl->get_all_Attributes();
    foreach my $a (@{$attr}) {
      print "\tName: ", $a->name(), "\n";
      print "\tCode: ", $a->code(), "\n";
      print "\tDesc: ", $a->description(), "\n";
      print "\tValu: ", $a->value(), "\n";
    }
    my @unip = @{ $transl->get_all_DBLinks('Uniprot%') };
    if (scalar(@unip) > 0) {
      foreach my $uniprot (@unip) {
	my $unip_id = $uniprot->display_id();
	##print UNIPROT "$transl_id\t$unip_id\n";
	## OK, add also
	## o uniprot_db: $uniprot->dbname();
	## o uniprot_info: $uniprot->info_text();
	my $dbn = $uniprot->dbname();
	$dbn =~ s/Uniprot\///g;
	my $descr = $uniprot->description();
	my $infot = $uniprot->info_text();
	my $stat = $uniprot->status();
	my $infotype = $uniprot->info_type();
	## mapping type.
	print "uniprot: ".$unip_id."\n";
	print " dbname: ".$dbn."\n";
	print " info_text ".$infot."\n";
	## Defines the method by which this ID was mapped (Uniprot ID was
	## matched to the Ensembl protein ID).
	print " info_type ".$infotype."\n";
      }
    }
  }
}
