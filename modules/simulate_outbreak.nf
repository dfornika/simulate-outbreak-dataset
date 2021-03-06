process simulate_variants {

  tag { prefix + ' / ' + num_snps }

  publishDir "${params.outdir}/${prefix}", pattern: "${prefix}*.{fa,vcf,tsv,csv}", mode: 'copy'

  input:
    path(seq)

  output:
    path("*.fa"), emit: seq
    path("*.fasta"), emit: selected
    tuple path("${prefix}.vcf"), path("${prefix}.variants.tsv"), emit: variants
    path("${prefix}_mutation_info.csv"), emit: mutation_info

  script:
  def generator = { String alphabet, int n ->
    new Random().with {
      (1..n).collect { alphabet[ nextInt( alphabet.length() ) ] }.join()
    }
  }
  prefix = generator((('A'..'F')+('0'..'9')).join(), 4)
  num_snps = Math.abs(new Random().nextInt() % params.max_snps_per_iteration.toInteger() + params.min_snps_per_iteration.toInteger())
  """
  cp \$(ls -1 *.fa | sort -R | head -n 1) selected.fasta
  echo 'parent_id,num_snps' > ${prefix}_mutation_info.csv
  echo \$(head -n 1 selected.fasta | tr -d '>' | cut -d ' ' -f 1),${num_snps} >> ${prefix}_mutation_info.csv
  simuG -refseq selected.fasta -snp_count ${num_snps} -prefix ${prefix}
  mv ${prefix}.simseq.genome.fa ${prefix}.fa
  sed -i '/^>/!s/.\\{70\\}/&\\n/g' ${prefix}.fa
  sed -i 's/^>.*/>${prefix}/g' ${prefix}.fa
  mv ${prefix}.refseq2simseq.map.txt ${prefix}.variants.tsv
  mv ${prefix}.refseq2simseq.SNP.vcf ${prefix}.vcf
  """
}

process simulate_reads {

    tag { seq_id }

    publishDir "${params.outdir}/${seq_id}", pattern: "${seq_id}_R{1,2}.fastq.gz", mode: 'copy'

    input:
    path(seq)

    output:
    tuple path("${seq_id}_R1.fastq.gz"), path("${seq_id}_R2.fastq.gz")

    script:
    seq_id = seq.toRealPath().toFile().text.split("\\n")[0].split(" ")[0].replace(">", "")
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
