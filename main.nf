#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.recursion=true

include { simulate_variants } from './modules/simulate_outbreak.nf'
include { simulate_reads } from './modules/simulate_outbreak.nf'


workflow select_and_mutate {
  take:
    ch_ref

  main:
    simulate_variants(ch_ref)
    simulate_reads(simulate_variants.out.selected)

  emit:
    simulate_variants.out.seq
}

workflow {
  inputs_ch = channel.fromList( [file(params.ref)] * params.iterations )

  main:
    select_and_mutate.scan(inputs_ch)
}
