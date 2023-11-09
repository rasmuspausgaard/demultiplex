#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
user="$USER"
date=new Date().format( 'yyMMdd' )


// Preset (default) parameters:
params.rundir                   ="${launchDir.baseName}"            // get basename of dir where script is started
params.outdir                   ='Results'                          // Default output folder.
params.genome                   ="hg38"                             // Default assembly
params.server                   ='kga01'                            // Default server

// Unset parameters: 
params.help                     =false
params.runfolder                =null   // required
params.samplesheet              =null   // required
params.DNA                      =null   // required
params.RNA                      =null   
params.skipAlign                =null
params.localStorage             =null
params.hg38v1                   =null
params.hg38v2                   =null
params.useBasesMask             =null
params.alignRNA                 =null
params.gatk                     =null

// Variables:
runID="${date}.${user}.demultiV1"
runtype="demultiAndPreprocessV1"


switch (params.server) {
    case 'lnx01':
        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
    break;
    case 'kga01':
        modules_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
    break;
}




def helpMessage() {
    log.info"""
    General info:
    This scripts requires at least 3 options: Path to runfolder, path to a samplesheet, and at least --DNA or --RNA parameter.
    There is no need to manually edit the samplesheet: The script will automatically separate DNA from RNA samples, and demultiplex them in parallel, if both --DNA and --RNA are set.
    
    PLEASE NOTE: The script will automatically perform preprocesssing and alignment of DNA samples. Fastq and aligned CRAM files will automatically be transferred to the long term storage location. This means no Fastq or aligned CRAM files will be found where the script is executed - only in the long term storage (dataArchive) location.
    
    The resulting CRAM files will be available from the data archive location from each server as read-only locations:
    ../dataArchive/tank_kga_external_archive/alignedData/{genomeversion}/novaRuns/runfolder
    
    The resulting FastQ files will be available from the data archive location from each server as read-only locations:
    ../dataArchive/tank_kga_external_archive/fastq_storage/novaRuns/runfolder

    Data can be stored locally in the folder where the script is started instead, by using the --localStorage option

    By default, hg38 (v3) is used for alignment! 
    
    NOTE:
    For NovaSeq runs, the index read lengths may vary, depending on the type of samples (DNA and/or RNA) being sequenced. 
    For DNA-only runs, the index read lengths are normally: I1 = 17bp, and I2 = 8 bp.
    For mixed runs (DNA and RNA samples), the index read lengths are normally: I1= 19bp, and I2 = 10bp.


    Usage:

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

    """.stripIndent()
}
if (params.help) exit 0, helpMessage()

/* TEST: BELOW ONLY IN MODULES FILE:
if (params.localStorage) {
aln_output_dir="${params.outdir}/"
fastq_dir="${params.outdir}/"
}
if (!params.localStorage) {
aln_output_dir="${tank_storage}/alignedData/${params.genome}/novaRuns/"
fastq_dir="${tank_storage}/fastq_storage/novaRuns/"
}
*/

//libParamLastLine-1="/data/shared/programmer/configfiles/library_param_last_line.v2.1.txt"
//libParamLastLine-2="/data/shared/programmer/configfiles/library_param_last_line.v2.2.txt"


channel
    .fromPath(params.runfolder)
    .map { it -> tuple (it.simpleName, it)}
    .set { runfolder_ch } 

channel
    .fromPath(params.runfolder)
    .map { it -> it.simpleName}
    .set {runfolder_simplename } 

if (!params.useBasesMask && !params.RNA) {
  dnaMask="Y*,I8nnnnnnnnn,I8,Y*"
}
if (!params.useBasesMask && params.RNA) {
  dnaMask="Y*,I8nnnnnnnnnnn,I8nn,Y*"
  rnaMask="Y*,I10nnnnnnnnn,I10,Y*"
}

if (params.useBasesMask) {
  dnaMask=params.useBasesMask
  params.DNA=true
}

log.info """\
=======================================
KGA Vejle demultiplex and Preprocess v1
=======================================
"""

channel
  .fromPath("${params.runfolder}/RunInfo.xml")
  .set {xml_ch}


channel.fromPath(params.samplesheet)
    .map { tuple(it.baseName,it) }
    .set {original_samplesheet}



include {      
         // Preprocess tools:
         prepare_DNA_samplesheet;
         bcl2fastq_DNA;
         prepare_RNA_samplesheet; 
         bcl2fastq_RNA;
         fastq_to_ubam;
         markAdapters;
         align;
         markDup_cram;
         } from "./modules/demulti_modules.nf" 


//from "${modules_dir}/demulti_modules.nf" 



workflow DEMULTIPLEX {
    take:
    original_samplesheet
    xml_ch
    main:
    if (params.DNA) {
        prepare_DNA_samplesheet(original_samplesheet)
        bcl2fastq_DNA(runfolder_ch, prepare_DNA_samplesheet.out, xml_ch)
    }
    if (params.RNA) {
        prepare_RNA_samplesheet(original_samplesheet)
        bcl2fastq_RNA(runfolder_ch, prepare_RNA_samplesheet.out, xml_ch)
    }
    emit: 
    dna_fastq=bcl2fastq_DNA.out.dna_fastq
}

workflow PREPROCESS {

    take:
    std_fq_input_ch
    
    main:
    fastq_to_ubam(std_fq_input_ch)
    markAdapters(fastq_to_ubam.out[0])
    align(markAdapters.out)
    markDup_cram(align.out)
    emit:
    finalAln=markDup_cram.out
    finalAln.view()
}

workflow {

    DEMULTIPLEX(original_samplesheet, xml_ch)
    
    if (params.DNA && !params.skipAlign){

        DEMULTIPLEX.out.dna_fastq.flatten()
        .filter {it =~/_R1_/}
        .map { tuple(it.baseName.tokenize('-').get(0)+"_"+it.baseName.tokenize('-').get(1),it) }
        .set { sampleid_R1 }

        DEMULTIPLEX.out.dna_fastq.flatten()
        .filter {it =~/_R2_/}
        .map { tuple(it.baseName.tokenize('-').get(0)+"_"+it.baseName.tokenize('-').get(1),it) }
        .set { sampleid_R2 }

        sampleid_R1.join(sampleid_R2)
        .combine(runfolder_simplename)
        .set { read_pairs_ch }

        read_pairs_ch
        .filter {it =~/AV1/}
        .set { av1_channel }

        read_pairs_ch
        .filter {it =~/WG4/ ||it =~/WG3/ ||it =~/EV8/}
        .set { rest_channel }

        av1_channel.concat(rest_channel)
        .set { std_fq_input_ch }

    PREPROCESS(std_fq_input_ch)
    }
}





   













/*

if (params.RNA && params.alignRNA) {

    rna_fq_out.flatten()
    .filter {it =~/_R1_/}
    .map { tuple(it.baseName.tokenize('-').get(0)+"_"+it.baseName.tokenize('-').get(1),it) }
    .into { RNA_sampleid_R1_ch1;RNA_sampleid_R1_ch2;RNA_sampleid_R1_ch3}

    rna_fq_out2.flatten()
    .filter {it =~/_R2_/}
    .map { tuple(it.baseName.tokenize('-').get(0)+"_"+it.baseName.tokenize('-').get(1),it) }
    .into { RNA_sampleid_R2_ch1;RNA_sampleid_R2_ch2;RNA_sampleid_R2_ch3}


    RNA_sampleid_R1_ch1.join(RNA_sampleid_R2_ch1)
	.combine(runfolder_simplename_value_ch2)
    .into {RNA_read_pairs_ch1;RNA_read_pairs_ch2;RNA_read_pairs_ch3}

    process STAR_align {
        errorStrategy 'ignore'
        tag "$sampleID"
        cpus 40
        maxForks 3
        publishDir "${aln_output_dir}/${runfolder_basename}/RNA/${sampleID}/", mode: 'copy', pattern: '*{Chimeric,ba,ba,cra,SJ,out}*'
    

        input:
        tuple val(sampleID), path(r1),path(r2), val(runfolder_basename) from RNA_read_pairs_ch1
        
        output:
        tuple val(sampleID), path("${sampleID}.STAR.Aligned.sortedByCoord.*.bam"),path("${sampleID}.STAR.Aligned.sortedByCoord.*.bai") into (split_cigar_input_bam, star_bam_ch2,star_bam_ch3,star_bam_ch4, preseq_input_bam,rseqc_input_bam,qualimapRNA_input_bam,qualimapBAMQC_input_bam,dupradar_input_bam, rnaseqc_input_bam, featurecounts_input_bam,htseq_count_input_bam, trinityfusion_input_bam)
        tuple val(sampleID),path("${sampleID}.STAR.Aligned.toTranscriptome.*.bam") into rsem_input_bam
    
        tuple val(sampleID), path("*.Chimeric.out.junction") into (trinityfusion_input_junction,trinity_splicing_input_junction)


        tuple val(sampleID), path("*.STAR.SJ.out.tab") into star_sjtab_out
        tuple val(sampleID), path("${sampleID}.forArriba.*")

        path("*.STAR.Log.*")
        path("*.forArriba.Log.out") into (trinity_collect_ch1,trinity_collect_ch2)
        script:
        """
        STAR --runThreadN ${task.cpus} \
            --genomeDir ${index_star} \
            --sjdbGTFfile ${gencode_gtf} \
            --readFilesIn ${r1} ${r2}\
            --outFileNamePrefix ${sampleID}.STAR. \
            --twopassMode Basic \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif --outSAMunmapped Within \
            --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMattrRGline ID:${date}.${user} LB:library PL:illumina PU:illumina SM:${sampleID} \
            --chimSegmentMin 12 \
            --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 \
            --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
            --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
            --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
            --outReadsUnmapped None \
            --alignInsertionFlush Right \
            --alignSplicedMateMapLminOverLmate 0 \
            --alignSplicedMateMapLmin 30 \
            --quantMode TranscriptomeSAM GeneCounts \
            --readFilesCommand zcat
    
        STAR --runThreadN ${task.cpus} \
            --genomeDir ${index_star} \
            --sjdbGTFfile ${gencode_gtf} \
            --readFilesIn ${r1} ${r2} \
            --outFileNamePrefix ${sampleID}.forArriba. \
            --twopassMode Basic \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif --outSAMunmapped Within \
            --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMattrRGline ID:${date}.${user} LB:library PL:illumina PU:illumina SM:${sampleID} \
            --chimSegmentMin 12 \
            --chimOutType WithinBAM \
            --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 \
            --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
            --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
            --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
            --outReadsUnmapped None \
            --readFilesCommand zcat
        
        bamtools index -in ${sampleID}.forArriba.Aligned.sortedByCoord.out.bam
        bamtools index -in ${sampleID}.STAR.Aligned.sortedByCoord.out.bam
        """
    }


}






 process markDup_v2 {
        errorStrategy 'ignore'
        maxForks 6
        tag "$sampleID"
        publishDir "${aln_output_dir}/${runfolder_basename}/", mode: 'copy', pattern: "*.BWA.MD.*"
        
        input:
        tuple val(sampleID), path(bam), val(runfolder_basename) from md_input1
        
        output:
        tuple val(sampleID), val(runfolder_basename), path("${sampleID}.${params.genome}.${genome_version}.BWA.MD.bam"), path("${sampleID}.${params.genome}.${genome_version}.BWA.MD*bai") into (final_bam_out,samblaster2)
        
        script:
        """
        samtools view -h ${bam} \
        | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=/data/TMP/TMP.${user}/ -o ${sampleID}.${params.genome}.${genome_version}.BWA.MD.bam /dev/stdin
        sambamba index ${sampleID}.${params.genome}.${genome_version}.BWA.MD.bam
        """
    }
    */
    