add_targets("chunk.vcf")

rule to_vcf:
    input: DT("split.flag")
    output: T("chunk.vcf")
    log: **LFS('to_vcf')
    resources: mem_mb=1024
    params:
        source_dir=PP("source_dir"),
        chrom=P("chrom"),
        chunk=P("chunk")

    shell: '''
        (
           chrdir=`dirname {input}`
           time \
               haplothread2VCF.py \
                      {params.source_dir}/metaData.txt \
                      $chrdir/{params.chunk}-pos.txt \
                      $chrdir/{params.chunk}-hpth.txt \
                      {params.chrom} {output} \
                      > {log.O} 2> {log.E} \
        ) 2> {log.T}
        '''