add_targets("split.flag")

rule split:
    output: T("split.flag")
    log: **LFS('split')
    params:
        bin_size=PP("positions_per_bin"),
        source_dir=PP("source_dir"),
        collection=PP("collection"),
        chrom=P("chrom")
    resources: mem_mb=100*1024
    shell: '''
        (time split_haplo.py \
                   {params.source_dir} {params.collection} \
                   {params.chrom} {params.bin_size} \
                   `dirname {output}` > {log.O} \
               && touch {output} \
               2> {log.E} \
        ) 2> {log.T}
    '''