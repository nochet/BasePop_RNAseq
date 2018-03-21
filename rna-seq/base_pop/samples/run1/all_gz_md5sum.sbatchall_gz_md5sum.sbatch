#!/bin/env perl6
#SBATCH --partition BioCompute
#SBATCH -J all_gz_md5sum 
#SBATCH -e all_gz_md5sum.e_%j
#SBATCH -o all_gz_md5sum.o_%j
#SBATCH --nodes 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 10G
#SBATCH --time=02:00:00
# Find all files in the current directory that en in '.gz'
my $max-cpu = 24;
my @gz-files = dir( test => / '.gz' $/ );    
my $cpu = min( $max-cpu, @gz-files.elems);
$cpu max= 1;
my @files = dir( test => / '.gz' $/ );

# For each file, create md5sum file if it doesn't already exist
for @files.race(batch=>1, degree=>$cpu) -> $file
{
    my $md5file = "$file.md5sum";

    if ! $md5file.IO.e
    {
        shell("md5sum $file > $md5file");
    }
}
