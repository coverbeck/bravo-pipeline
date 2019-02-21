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

    # Prepare Coverage #
    Array[File] inputCramFiles
    Int chromosome

    # Prepare CRAM #
    File sampleLocationFile
    File sampleIndexLocationFile

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
    }
    call addPercentiles {
        input: inputVCF = addCaddScores.out,
            inputVCFIndex = computePercentiles.outVariantPercentileIndex,
            variantPercentiles = computePercentiles.outVariantPercentile

    }

    ####################
    # Prepare Coverage #
    ####################

    scatter (file in inputCramFiles) {
        call extractDepth {
            input: inputCramFile = file,
            chromosome = chromosome,
            refFasta = refFasta
        }
    }
    call aggrBasePair {
        input: inputFiles = extractDepth.outDepth,
            chromosome = chromosome
    }

    ################
    # Prepare Cram #
    ################

    call extractId {
        input: inputVCF = inputVCF,
            samplesFile = samplesFile
    }
    call prepareSequences {
        input: inputVCF = inputVCF,
            sampleLocationFile = sampleLocationFile,
            sampleIndexLocationFile = sampleIndexLocationFile
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
        File outVariantPercentileIndex = "${infoField}.variant_percentile.vcf.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: threads
        bootDiskSizeGb: "150"
    }
}

task addPercentiles {
    File inputVCF
    Array[File] inputVCFIndex
    Array[File] variantPercentiles

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

task extractDepth {
    File inputCramFile
    Int chromosome
    File refFasta
    String sample = basename(inputCramFile, ".bam")

    command {
        samtools view -q 20 -F 0x0704 -uh ${inputCramFile} ${chromosome} | \
        samtools calmd -uAEr - ${refFasta} | \
        bam clipOverlap --in -.ubam --out -.ubam | \
        samtools mpileup -f ${refFasta} -Q 20 -t DP - | \
        cut -f1-4 | \
        bgzip > ${chromosome}.${sample}.depth.gz \
        && tabix ${chromosome}.${sample}.depth.gz
    }
    output {
        File outDepth = "${chromosome}.${sample}.depth.gz"
        File outIndex = "${chromosome}.${sample}.depth.gz.tbi"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}

task aggrBasePair {
    Array[File] inputFiles
    Int chromosome
    # Not splitting by BP for now
    Int startBP = 0
    Int endBP = 999999999

    command {
        create_coverage.py -i ${write_lines(inputFiles)} aggregate -c ${chromosome} -s ${startBP} -e ${endBP} | \
        bgzip -c > ${chromosome}.${startBP}.${endBP}.json.gz
    }
    output {
        File outAggrBasePair = "${chromosome}.${startBP}.${endBP}.json.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}

task extractId {
    File inputVCF
    File samplesFile
    String sample = basename(inputVCF, ".vcf.gz")

    command {
        RandomHetHom -k 5 -e 1985 -i ${inputVCF} -s ${samplesFile} -o ${sample}.vcf.gz
    }
    output {
        File out = "${sample}.vcf.gz"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}

task prepareSequences {
    File inputVCF
    File sampleLocationFile
    File sampleIndexLocationFile

    # Need to read index and bam files into WDL task scope 
    Array[File] sampleFiles = read_lines(sampleLocationFile)
    Array[File] sampleFilesIndex = read_lines(sampleIndexLocationFile)

    command {
        prepare_sequences.py cram -i ${inputVCF} -c ${sep=' ' sampleFiles} -w 100 -o combined.cram
    }
    output {
        File out = "combined.cram"
    }
    runtime {
        docker: "statgen/bravo-pipeline:latest"
        cpu: "1"
        bootDiskSizeGb: "50"
    }
}
