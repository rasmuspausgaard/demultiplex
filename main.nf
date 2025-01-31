#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
user="$USER"
date=new Date().format( 'yyMMdd' )


// Preset (default) parameters:
params.rundir                   ="${launchDir.baseName}"            // get basename of dir where script is started
params.outdir                   ='Results'                          // Default output folder.
params.genome                   ="hg38"                             // Default assembly
params.server                   ='lnx01'                            // Default server

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
params.nomail             			=null
params.keepwork 	            	=null

// Variables:
runID="${date}.${user}.demultiV1"
runtype="demultiAndPreprocessV1"


def helpMessage() {
    log.info"""
    General info:
    This scripts requires at least 3 options: Path to runfolder, path to a samplesheet, and at least --DNA or --RNA parameter.
    There is no need to manually edit the samplesheet: The script will automatically separate DNA from RNA samples, and demultiplex them in parallel, if both --DNA and --RNA are set.
    
    PLEASE NOTE: The script will automatically perform preprocesssing and alignment of DNA samples. Fastq and aligned CRAM files will automatically be transferred to the long term storage location. This means no Fastq or aligned CRAM files will be found where the script is executed - only in the long term storage (dataArchive) location.
    
    The resulting CRAM files will be available from the data archive location from each server as read-only locations:
    ../dataArchive/lnx02/alignedData/{genomeversion}/novaRuns/runfolder
    
    The resulting FastQ files will be available from the data archive location from each server as read-only locations:
    ../dataArchive/lnx02/fastq_storage/novaRuns/runfolder

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

      --server              Choose server to run script from (lnx01 or lnx02)
                                Default: lnx01

      --genome               hg19 or hg38
                                Default: hg38v3

      --keepwork            keep the workfolder generated by the nextflow script.
                                Default: not set - removes the Work folder generated by nextflow

      --nomail              Does not send a mail-message when completing a script
                                Default: not set - sends mail message if the user is mmaj or raspau and only if the script has been running longer than 20 minutes.

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
        .filter {it =~/AV1/||it =~/MV1/}
        .set { av1_channel }

        read_pairs_ch
        .filter {it =~/WG4/ ||it =~/WG3/ ||it =~/EV8/||it =~/LIB/}
        .set { rest_channel }

        av1_channel.concat(rest_channel)
        .set { std_fq_input_ch }

    PREPROCESS(std_fq_input_ch)
    }
}

/*
workflow.onComplete {

    // Extract the first six digits from the samplesheet name
    def samplesheetName = new File(params.samplesheet).getName()
    def samplesheetDate = samplesheetName.find(/\d{6}/)

    // Only send email if --nomail is not specified, duration is longer than 20 minutes, the script executed succesfully, and if the user is mmaj or raspau.
    if (!params.nomail && workflow.duration > 1200000 && workflow.success) {
        if (System.getenv("USER") in ["raspau", "mmaj"]) {
            
            def workDirMessage = params.keepwork ? "WorkDir             : ${workflow.workDir}" : "WorkDir             : Deleted"
            
            def body = """\
            Pipeline execution summary
            ---------------------------
            Demultiplexing of sequencing run: ${samplesheetDate}
            Duration            : ${workflow.duration}
            Success             : ${workflow.success}
            ${workDirMessage}
            OutputDir           : ${params.outdir ?: 'Not specified'}
            Exit status         : ${workflow.exitStatus}
            """.stripIndent()


            // Send email using the built-in sendMail function
            sendMail(to: 'Mads.Jorgensen@rsyd.dk,Rasmus.Hojrup.Pausgaard@rsyd.dk', subject: 'Demultiplexing pipeline Update', body: body)
        }
    }

    // Handle deletion of WorkDir based on --keepwork parameter
    if (!params.keepwork && workflow.duration > 1200000 && workflow.success) {
        println("Deleting work directory: ${workflow.workDir}")
        "rm -rf ${workflow.workDir}".execute()
    }
}
*/
