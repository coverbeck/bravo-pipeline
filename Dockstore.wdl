workflow bravoDataPrep {
    File inputVCF
    File samplesFile
    Int bufferSize
    String assembly
    String lofOptions
    File refDir
    File refFasta
    File cadScores

    call computeAlleleCountsAndHistograms {
        input: inputVCF = inputVCF,
            samplesFile = samplesFile,
    }
    call variantEffectPredictor {
        input: inputVCF = computeAlleleCountsAndHistograms.out,
            assembly = assembly,
            lofOptions = lofOptions,
            bufferSize = bufferSize,
            refDir = refDir,
            refFasta = refFasta
    }
    call annotateVCF {
        input: inputVCF = variantEffectPredictor.out,
            cadScores = cadScores
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
        docker: "statgen/bravo-pipeline:latest"
    }

}

task  variantEffectPredictor {
    File inputVCF
    String assembly
    String lofOptions
    Int bufferSize
    File refDir
    File refFasta

    command {
        vep -i ${inputVCF} \
        --plugin LoF,${lofOptions} \
        --dir_cache ${refDir} \
        --fasta ${refFasta} \
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
    output {
        File out = "vep-out.gz"
    }
    runtime {
        docker: "ensemblorg/ensembl-vep:release_95.1"
        cpu: "1"
        bootDiskSizeGb: "150"
    }

}

task annotateVCF {
    File inputVCF
    File cadScores

    command {
        add_cadd_scores.py -i ${inputVCF} -c ${cadScores} -o cad-out.vcf.gz
    }
    output {
        File out = "cad-out.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "150"
    }
}
