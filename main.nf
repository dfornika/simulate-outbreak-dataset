#!/usr/bin/env nextflow

import nextflow.extension.CH
import static nextflow.extension.DataflowHelper.newOperator
import groovyx.gpars.dataflow.operator.ChainWithClosure
import groovyx.gpars.dataflow.operator.CopyChannelsClosure

nextflow.enable.dsl = 2

include { simulate_variants } from './modules/simulate_outbreak.nf'
include { simulate_variants as simulate_variants_from_ancestor } from './modules/simulate_outbreak.nf'
include { simulate_reads } from './modules/simulate_outbreak.nf'
include { select_ancestor } from './modules/simulate_outbreak.nf'

// https://github.com/nextflow-io/nextflow/issues/1766#issuecomment-831286103
def attach(source, target) {
    newOperator([source.createReadChannel()], [target],
                new ChainWithClosure(new CopyChannelsClosure()))
}

workflow {
  ch_feedback = CH.create()
  condition = { it[1] >= params.iterations }
  ch_ref = Channel.fromPath(params.ref).combine(Channel.of(0)).mix(ch_feedback.until(condition))

  main:
    simulate_variants(ch_ref)
    simulate_reads(simulate_variants.out.seq)
    select_ancestor(simulate_variants.out.seq.map{ it -> it[0] }.buffer(size: params.buffer_size).collect().combine(simulate_variants.out.seq.map{ it -> it[1]}).view())
    attach(simulate_variants.out.seq, ch_feedback)
}
