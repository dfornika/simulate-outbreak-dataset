process simulate_variants {

  tag { seq_id + ' / ' + prefix + ' / ' + num_snps }

  publishDir "${params.outdir}", pattern: "${prefix}.{fa,vcf,variants.tsv}", mode: 'copy'

  input:
  tuple path(seq), val(iteration)

  output:
  tuple path("${prefix}.fa"), val(iteration), emit: seq
  tuple path("${prefix}.vcf"), path("${prefix}.variants.tsv"), emit: variants

  script:
  def generator = { String alphabet, int n ->
    new Random().with {
      (1..n).collect { alphabet[ nextInt( alphabet.length() ) ] }.join()
    }
  }
  prefix = generator((('A'..'F')+('0'..'9')).join(), 4)
  num_snps = Math.abs(new Random().nextInt() % params.max_snps_per_iteration.toInteger() + params.min_snps_per_iteration.toInteger())
  seq_id = seq.baseName
  iteration += 1
  """
  simuG -refseq ${seq} -snp_count ${num_snps} -prefix ${prefix}
  mv ${prefix}.simseq.genome.fa ${prefix}.fa
  sed -i '/^>/!s/.\\{70\\}/&\\n/g' ${prefix}.fa
  sed -i 's/^>.*/>${prefix}/g' ${prefix}.fa
  mv ${prefix}.refseq2simseq.map.txt ${prefix}.variants.tsv
  mv ${prefix}.refseq2simseq.SNP.vcf ${prefix}.vcf
  """
}

process simulate_reads {

    tag { seq_id }

    publishDir "${params.outdir}", pattern: "${seq_id}_R{1,2}.fastq.gz", mode: 'copy'

    input:
    tuple path(seq), val(iteration)

    output:
    tuple path("${seq_id}_R1.fastq.gz"), path("${seq_id}_R2.fastq.gz")

    script:
    seq_id = seq.baseName
    mean_fragment_length = 300
    stdev_fragment_length = 50
    read_length = 250
    fold_coverage = 5
    """
    echo "${seq_id}" > seq_id.txt
    art_illumina \
      --paired \
      --in ${seq} \
      --fcov ${fold_coverage} \
      --len ${read_length} \
      --mflen ${mean_fragment_length} \
      --sdev ${stdev_fragment_length} \
      --out ${seq_id}_R
    mv ${seq_id}_R1.fq ${seq_id}_R1.fastq
    mv ${seq_id}_R2.fq ${seq_id}_R2.fastq
    gzip ${seq_id}_R1.fastq
    gzip ${seq_id}_R2.fastq
    """
}

process select_ancestor {

  executor 'local'

  input:
  tuple path(seq), val(iteration)

  output:
  tuple path("selected.fa"), val(iteration)

  script:
  iteration += 1
  """
  cp \$(ls -1 *.fa | sort -R | head -n 1) selected.fa
  """
}