workflow bravoDataPrep {
    # Prepare VCF #
    File inputVCF
    File samplesFile
    Int bufferSize
    String assembly
    String lofOptions
    File refDir
    File refFasta
    File cadScores
    File cadIndex

    # Prepare percentiles #
    Array[String] infoField
    Int threads
    Int numberPercentiles
    String description

    ###############
    # Prepare VCF #
    ###############

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
    call addCaddScores {
        input: inputVCF = variantEffectPredictor.out,
            cadScores = cadScores,
            cadIndex = cadIndex
    }

    #######################
    # Prepare percentiles #
    #######################

    scatter (field in infoField) {
        call computePercentiles {
            input: inputVCF = addCaddScores.out,
                infoField = field,
                threads = threads,
                numberPercentiles = numberPercentiles,
                description = description
        }
        call indexVCF {
            input: variantPercentileVCF = computePercentiles.outVariantPercentile
        }
    }
}

task computeAlleleCountsAndHistograms {
    File inputVCF
    File samplesFile

    command {
        ComputeAlleleCountsAndHistograms -i ${inputVCF} -s ${samplesFile} -o computeAlleleCtHst.vcf.gz
    }
    output {
        File out = "computeAlleleCtHst.vcf.gz"
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
        --compress_output bgzip \
        --no_stats \
        -o variantEP.vcf.gz
    }
    output {
        File out = "variantEP.vcf.gz"
    }
    runtime {
        docker: "ensemblorg/ensembl-vep:release_95.1"
        cpu: "1"
        bootDiskSizeGb: "150"
    }

}

task addCaddScores {
    File inputVCF
    File cadScores
    File cadIndex

    command {
        add_cadd_scores.py -i ${inputVCF} -c ${cadScores} -o annotated.vcf.gz
    }
    output {
        File out = "annotated.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "150"
    }
}

task computePercentiles {
    File inputVCF
    String infoField
    Int threads
    Int numberPercentiles
    String description

    command <<<
        ComputePercentiles -i ${inputVCF} \
        -m ${infoField} \
        -t ${threads} \
        -p ${numberPercentiles} \
        -d ${description} \
        -o ${infoField}

        tabix ${infoField}.variant_percentile.vcf.gz
    >>>
    output {
        File outAllPercentiles = "${infoField}.all_percentiles.json.gz"
        File outVariantPercentile = "${infoField}.variant_percentile.vcf.gz"
        File outIndexedVariantPercentile = "${infoField}.variant_percentile.vcf.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: threads
        bootDiskSizeGb: "150"
    }
}

task addPercentiles {
    File inputVCF
    File inputVCFIndex
    Array[String] variantPercentiles

    command {
        add_percentiles.py -i ${inputVCF} -p ${sep=' ' variantPercentiles} -o percentiles.vcf.gz
    }
    output {
        File out = "percentiles.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "150"
    }

}
