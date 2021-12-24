#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.recursion=true

include { simulate_variants } from './modules/simulate_outbreak.nf'
include { simulate_reads } from './modules/simulate_outbreak.nf'
include { select_ancestor } from './modules/simulate_outbreak.nf'


workflow mutate_and_select {
  take:
    ch_ref

  main:
    simulate_variants(ch_ref)
    simulate_reads(simulate_variants.out.selected)

  emit:
    simulate_variants.out.seq
}

workflow {
  ref = file(params.ref)

  main:
    mutate_and_select.recurse(ref).times(params.iterations).view()
    
    
}
