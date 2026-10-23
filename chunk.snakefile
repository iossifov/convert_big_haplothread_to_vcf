add_targets("chunk.vcf.gz")

rule gzip:
    input: T("chunk.vcf")
    output: T("chunk.vcf.gz")
    log: **LFS('gzip')
    shell: "(time gzip {input} > {log.O} 2> {log.E}) 2> {log.T}"

rule to_vcf:
    input: DT("split.flag")
    output: temp(T("chunk.vcf"))
    log: **LFS('to_vcf')
    resources: mem_mb=1024
    params:
        source_dir=PP("source_dir"),
        chrom=P("chrom"),
        chunk=P("chunk")

    shell: '''
        (
           ulimit -n 10240
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