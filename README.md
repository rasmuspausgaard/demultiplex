# NGS demultiplex and preprocessing script, KG Vejle

## General info:
This scripts requires at least 3 options: Path to runfolder, path to a samplesheet, and at least --DNA or --RNA parameter.
There is no need to manually edit the samplesheet: The script will automatically separate DNA from RNA samples, and demultiplex them in parallel, if both --DNA and --RNA are set.

PLEASE NOTE: The script will automatically perform preprocesssing and alignment of DNA samples. Fastq and aligned CRAM files will automatically be transferred to the long term storage location. This means no Fastq or aligned CRAM files will be found where the script is executed - only in the long term storage (dataArchive) location.

The resulting CRAM files will be available from the data archive location from each server as read-only locations:
../dataArchive/tank_kga_external_archive/alignedData/{genomeversion}/novaRuns/runfolder

The resulting FastQ files will be available from the data archive location from each server as read-only locations:
../dataArchive/tank_kga_external_archive/fastq_storage/novaRuns/runfolder

Data can be stored locally in the folder where the script is started instead by using the --localStorage option

By default, hg38 (v3) is used for alignment!

## Usage:

Main options:
  --help                print this help message

  --runfolder           Path to runfolder (required)

  --samplesheet         Path to samplesheet (required)

  --localStorage        If set, all output data (fastq and CRAM) will be stored locally in output folder where the script is started.
                            Default: Not set. Data are by default stored at the archive (see description above)

  --DNA                 Demultiplex DNA samples

  --RNA                 Demultiplex RNA samples

  --useBasesMask        manually set "--use-bases-mask" parameter for demultiplexing.
                        Default for demultiplexing DNA-only runs (I1=8 bp and I2=17 bp):
                        DNA samples: Y*,I8nnnnnnnnn,I8,Y
                        Default for demultiplexing both DNA- and RNA-samples (I1=10 bp and I2=19 bp):
                        DNA samples:"Y*,I8nnnnnnnnnnn,I8nn,Y*"
                        RNA samples: "Y*,I10nnnnnnnnn,I10,Y*"

  --skipAlign           Do not perform preprocessing and alignment (i.e. only demultiplexing)
                            Default: Not set (i.e. run preprocessing and alignment)

  --alignRNA            Align RNA samples (STAR 2-pass). Alignment will be run twice using the reccomended parameters for e.g. Arriba or STAR-fusion.
                            Default: Do not align RNA samples.

  --server              Choose server to run script from (lnx01 or kga01)
                            Default: kga01

  --genome               hg19 or hg38
                            Default: hg38v3
