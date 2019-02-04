workflow bravoDataPrep {
    Boolean Offline
    File InputVCF
    File SamplesFile
    Int BufferSize
    String Assembly

    call computeAlleleCountsAndHistograms {
        input: inputVCF = InputVCF,
            samplesFile = SamplesFile,
    }
    if (Offline) {
    call variantEffectPredictorOffline {
        input: inputVCF = computeAlleleCountsAndHistograms.out,
            assembly = Assembly,
            bufferSize = BufferSize,
    }
    }
    if (!Offline) {
        call variantEffectPredictor {
            input: inputVCF = computeAlleleCountsAndHistograms.out,
            assembly = Assembly,
            bufferSize = BufferSize
        }
    }
}

task computeAlleleCountsAndHistograms {

    File inputVCF
    File samplesFile

    command {
        ComputeAlleleCountsAndHistograms -i ${inputVCF} -s ${samplesFile} -o test.vcf.gz
    }
    output {
        File out = "test.vcf.gz"
    }
    runtime {
        docker: 'statgen/bravo-pipeline:latest'
    }

}

task  variantEffectPredictorOffline {

    File inputVCF
    String assembly
    Int bufferSize

    command {
        vep -i ${inputVCF} \
        --assembly ${assembly} \
        --cache \
        --offline \
        --vcf \
        --sift b \
        --polyphen b \
        --ccds \
        --uniprot \
        --hgvs \
        --symbol \
        --numbers \
        --domains \
        --regulatory \
        --canonical \
        --protein \
        --biotype \
        --af \
        --af_1kg \
        --pubmed \
        --shift_hgvs 0 \
        --allele_number \
        --format vcf \
        --force \
        --buffer_size ${bufferSize} \
        --compress_output gzip \
        --no_stats \
        -o vep-out.gz
    }
# Need to change the output from the File type to a string
    output {
        File out = "vep-out.gz"
    }
    runtime {
        docker: 'ensemblorg/ensembl-vep:release_95.1'
    }

}

task  variantEffectPredictor {

    File inputVCF
    String assembly
    Int bufferSize

    command {
        vep -i ${inputVCF} \
        --assembly ${assembly} \
        --vcf \
        --sift b \
        --database \
        --polyphen b \
        --ccds \
        --uniprot \
        --hgvs \
        --symbol \
        --numbers \
        --domains \
        --regulatory \
        --canonical \
        --protein \
        --biotype \
        --af \
        --shift_hgvs 0 \
        --allele_number \
        --format vcf \
        --force \
        --buffer_size ${bufferSize} \
        --compress_output gzip \
        --no_stats \
        -o vep-out.gz
    }

    output {
        File out = "vep-out.gz"
    }
    runtime {
        docker: 'ensemblorg/ensembl-vep:release_95.1'
    }

}

task annotateVCF {

    command {

    }
    output {

    }
    runtime {
        docker: 'statgen/bravo-pipeline:latest'
    }

}
