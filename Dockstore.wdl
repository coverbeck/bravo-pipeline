workflow bravoDataPrep {
    File InputVCF
    File SamplesFile

    call computeAlleleCountsAndHistograms {
        input: inputVCF = InputVCF,
            samplesFile = SamplesFile,
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

task  variantEffectPredictor {

    File inputVCF
    String lofOptions

    command {
        vep -i ${inputVCF} \
        --plugin LoF ${lofOptions} \
        --assembly GRCh38 \
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
        --buffer_size 100000 \
        --compress_output gzip \
        --no_stats -o ${inputVCF}
    }
    output {
        File vcfOutput = "${inputVCF}-vep.gz"
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
