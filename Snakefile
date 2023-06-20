from os import path, getcwd

# Run the pipeline like so:
# snakemake -pr --cores 10 --keep-going --rerun-incomplete --restart-times 3
# snakemake -s Snakefile --rulegraph |dot -Tpng > rulegraph.png
configfile: "config/config.yaml"

onsuccess:
    print("Workflow completed without any error")

onerror:
    print("An error occurred")


DBNAME=config['parameters']['blast']['dbname']
DB="13.MakeblastDB/{DB}".format(DB=DBNAME)
ASSEMBLER=config['ASSEMBLER']

RULES=["Count_seqs_pre_trim", "QC_pre_trim", "SummarizeQC_pre_trim", "Trimmomatic_trim_adaptors",
       "Trim_primers", "Trim_galore_trim_adaptors", "Trim_galore_trim_adaptors2", 
       "SummarizeQC_post_adaptor_trim", "Count_seqs_post_trim", "Assemble_rRNA", "Detect_chimeras", 
       "Filter_seqs", "MakeblastDB", "BlastSeqs", "SummarizeBlastResult"]

rule all:
    input:
        "02.Count_seqs_pre_trim/seqs_stat.txt",
        "04.SummarizeQC_pre_trim/multiqc_report.html",
        "08.SummarizeQC_post_adaptor_trim/multiqc_report.html",
        "09.Count_seqs_post_trim/seqs_stat.txt",
        multiext(DB, ".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"),
        expand("14.BlastSeqs/{sample}/{sample}.tsv", sample=config['samples']),
        expand("15.SummarizeBlastResult/{sample}/{sample}.txt", sample=config['samples'])



rule Make_logs_directories:
    output:
        directory("logs/Trim_galore_trim_adaptors/"),
        directory("logs/QC_pre_trim/"),
        directory("logs/SummarizeQC_post_adaptor_trim/"),
        directory("logs/Trim_primers/"),
        directory("logs/SummarizeQC_pre_trim/"),
        directory("logs/Trimmomatic_trim_adaptors/")
    threads: 1
    shell:
        """
         [ -d logs/ ] || mkdir -p logs/
         cd logs/
         for RULE in {RULES}; do
          [ -d ${{RULE}}/ ] || mkdir -p ${{RULE}}/
         done
        """


rule Count_seqs_pre_trim:
    input: expand(["01.raw_data/{sample}/{sample}_R1.fastq.gz", "01.raw_data/{sample}/{sample}_R2.fastq.gz"], sample=config['samples'])
    output: "02.Count_seqs_pre_trim/seqs_stat.txt"
    params:
        program=config['programs_path']['seqkit']
    shell:
        """
        # Get the stats on the sequences using seqkit
        {params.program} stats {input} > temp.txt

         # Sort the sequence statistics
         (sed -n '1p' temp.txt; awk 'NR>1{{print}}' temp.txt | \
           sort -V -k1,1) > {output} \
           && rm temp.txt
        """

rule QC_pre_trim:
    input:
        rules.Make_logs_directories.output,
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz"
    output:
        forward_html="03.QC_pre_trim/{sample}/{sample}_R1_fastqc.html",
        reverse_html="03.QC_pre_trim/{sample}/{sample}_R2_fastqc.html"
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=1
    log: "logs/QC_pre_trim/{sample}/{sample}.log"
    threads: 1
    shell:
        "{params.program} --outdir {params.out_dir} "
        "--threads {params.threads} {input.forward} {input.rev} > {log} 2>&1"


rule SummarizeQC_pre_trim:
    input:
        forward_html=expand("03.QC_pre_trim/{sample}/{sample}_R1_fastqc.html",
                            sample=config['samples']),
        rev_html=expand("03.QC_pre_trim/{sample}/{sample}_R2_fastqc.html",
                            sample=config['samples'])
    output: "04.SummarizeQC_pre_trim/multiqc_report.html"
    log: "logs/SummarizeQC_pre_trim/multiqc.log"
    params:
        program=config['programs_path']['multiqc'],
        in_dir=lambda w, input: path.dirname(input[0]).split("/")[0],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

          {params.program} \
              --interactive \
              -f {params.in_dir} \
              -o {params.out_dir}  > {log} 2>&1
        """



adaptors=config['parameters']['trimmomatic']['adaptors']
min_length=config['parameters']['trimmomatic']['min_len']
rule Trimmomatic_trim_adaptors:
    input:
        rules.Make_logs_directories.output,
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz",
        log_dirs=rules.Make_logs_directories.output
    output:
        r1="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        r2="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R1.unpaired.fastq.gz",
        r2_unpaired="05.Trimmomatic_trim_adaptors/{sample}/{sample}_R2.unpaired.fastq.gz"
    log:
        "logs/Trimmomatic_trim_adaptors/{sample}/{sample}.log"
    params:
        program=config['programs_path']['trimmomatic'],
        trimmer="ILLUMINACLIP:{adaptors}:2:30:10"
                " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20"
                " MINLEN:{min_length}".format(adaptors=adaptors,
                                          min_length=min_length)
    threads: 1
    resources:
        mem_mb=1024
    shell:
        "{params.program} PE "
        "-threads {threads} "
        "{input.forward} {input.rev} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer} > {log} 2>&1 "




# Trim Primers using cutadapt
rule Trim_primers:
    input:
        rules.Make_logs_directories.output,
        forward_reads=rules.Trimmomatic_trim_adaptors.output.r1,
        rev_reads=rules.Trimmomatic_trim_adaptors.output.r2
    output:
        forward_reads="06.Trim_primers/{sample}/{sample}_R1.fastq.gz",
        rev_reads="06.Trim_primers/{sample}/{sample}_R2.fastq.gz"
    log: "logs/Trim_primers/{sample}/{sample}.log"
    params:
        program=config['programs_path']['cutadapt'],
        forward_primer=config['parameters']['cutadapt']['forward_primer'],
        rev_primer=config['parameters']['cutadapt']['reverse_primer'],
        minimum_length=config['parameters']['cutadapt']['minimum_length'],
        quality_cutoff=config['parameters']['cutadapt']['quality_cutoff'],
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        # -a TGGAATTCTCGGGTGCCAAGG sequence here is the small RNA 3' adaptor
        # https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/RNA/Small-RNA/TruSeqSmallRNA.htm
        # -A CTGTCTCTTATACAC sequece here is the nextera transposa sequence 
        {params.program} \
              -g '{params.forward_primer}' \
              -G '{params.rev_primer}' \
              -o {output.forward_reads} \
              -p {output.rev_reads} \
              --minimum-length  {params.minimum_length} \
              --quality-cutoff  {params.quality_cutoff} \
             {input.forward_reads} {input.rev_reads} > {log} 2>&1

        """

# Trim adaptors and quality check using cutadapt and factqc with Trim galore
rule Trim_galore_trim_adaptors:
    input: 
        rules.Make_logs_directories.output,
        forward=rules.Trim_primers.output.forward_reads,
        rev=rules.Trim_primers.output.rev_reads
    output:
        forward_reads="07.Trim_galore_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        rev_reads="07.Trim_galore_trim_adaptors/{sample}/{sample}_R2.fastq.gz",
        forward_html="07.Trim_galore_trim_adaptors/{sample}/{sample}_R1_fastqc.html",
        rev_html="07.Trim_galore_trim_adaptors/{sample}/{sample}_R2_fastqc.html"
    log: "logs/Trim_galore_trim_adaptors/{sample}/{sample}.log"
    threads: 1
    params:
        program=config['programs_path']['trim_galore'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    shell:
        """ 
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
           -o {params.out_dir} \
           --fastqc  \
           --paired {input.forward} {input.rev} > {log} 2>&1

         #Rename the files
         #Fastq files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1.fq.gz {params.out_dir}/{wildcards.sample}_R1.fastq.gz
        mv {params.out_dir}/{wildcards.sample}_R2_val_2.fq.gz {params.out_dir}/{wildcards.sample}_R2.fastq.gz
         #HTML files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.html {params.out_dir}/{wildcards.sample}_R1_fastqc.html
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.html {params.out_dir}/{wildcards.sample}_R2_fastqc.html
         #Zip files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.zip {params.out_dir}/{wildcards.sample}_R1_fastqc.zip
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.zip {params.out_dir}/{wildcards.sample}_R2_fastqc.zip
        """

# Trim adaptors and quality check using cutadapt and factqc with Trim galore a second time
rule Trim_galore_trim_adaptors2:
    input: 
        rules.Make_logs_directories.output,
        forward=rules.Trim_galore_trim_adaptors.output.forward_reads,
        rev=rules.Trim_galore_trim_adaptors.output.rev_reads
    output:
        forward_reads="07.Trim_galore_trim_adaptors2/{sample}/{sample}_R1.fastq.gz",
        rev_reads="07.Trim_galore_trim_adaptors2/{sample}/{sample}_R2.fastq.gz",
        forward_unpaired="07.Trim_galore_trim_adaptors2/{sample}/{sample}_R1_unpaired.fastq.gz",
        rev_unpaired="07.Trim_galore_trim_adaptors2/{sample}/{sample}_R2_unpaired.fastq.gz",
        forward_html="07.Trim_galore_trim_adaptors2/{sample}/{sample}_R1_fastqc.html",
        rev_html="07.Trim_galore_trim_adaptors2/{sample}/{sample}_R2_fastqc.html"
    log: "logs/Trim_galore_trim_adaptors2/{sample}/{sample}.log"
    threads: 1
    params:
        program=config['programs_path']['trim_galore'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    shell:
        """ 
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
           -o {params.out_dir} \
           --fastqc  \
           --retain_unpaired \
           --paired {input.forward} {input.rev} > {log} 2>&1

         #Rename the files
         # Fastq files
        # Paired
        mv {params.out_dir}/{wildcards.sample}_R1_val_1.fq.gz {params.out_dir}/{wildcards.sample}_R1.fastq.gz
        mv {params.out_dir}/{wildcards.sample}_R2_val_2.fq.gz {params.out_dir}/{wildcards.sample}_R2.fastq.gz

        # Unpaired
        mv {params.out_dir}/{wildcards.sample}_R1_unpaired_1.fq.gz {params.out_dir}/{wildcards.sample}_R1_unpaired.fastq.gz
        mv {params.out_dir}/{wildcards.sample}_R2_unpaired_2.fq.gz {params.out_dir}/{wildcards.sample}_R2_unpaired.fastq.gz
         #HTML files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.html {params.out_dir}/{wildcards.sample}_R1_fastqc.html
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.html {params.out_dir}/{wildcards.sample}_R2_fastqc.html
         #Zip files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.zip {params.out_dir}/{wildcards.sample}_R1_fastqc.zip
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.zip {params.out_dir}/{wildcards.sample}_R2_fastqc.zip
        """



# Aggregate and summarize the quality check using multiqc
rule SummarizeQC_post_adaptor_trim:
    input:
        forward_html=expand("07.Trim_galore_trim_adaptors2/{sample}/{sample}_R1_fastqc.html",
                            sample=config['samples']),
        rev_html=expand("07.Trim_galore_trim_adaptors2/{sample}/{sample}_R2_fastqc.html",
                            sample=config['samples'])
    output: "08.SummarizeQC_post_adaptor_trim/multiqc_report.html"
    log: "logs/SummarizeQC_post_adaptor_trim/multiqc.log"
    params:
        program=config['programs_path']['multiqc'],
        in_dir=lambda w, input: path.dirname(input[0]).split("/")[0],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

          {params.program} \
              --interactive \
              -f {params.in_dir} \
              -o {params.out_dir}  > {log} 2>&1
        """




rule Count_seqs_post_trim:
    input: expand(["07.Trim_galore_trim_adaptors2/{sample}/{sample}_R1.fastq.gz", "07.Trim_galore_trim_adaptors2/{sample}/{sample}_R2.fastq.gz"], sample=config['samples'])
    output: "09.Count_seqs_post_trim/seqs_stat.txt"
    params:
        program=config['programs_path']['seqkit']
    shell:
        """
        # Get the stats on the sequences using seqkit
        {params.program} stats {input} > temp.txt

         # Sort the sequence statistics
         (sed -n '1p' temp.txt; awk 'NR>1{{print}}' temp.txt | \
           sort -V -k1,1) > {output} \
           && rm temp.txt
        """

# Assemble the genome using spades
# --careful =  minimize the number of mismatches in the contigs
# -o = output folder
# -1 = forward reads
# -2 = rev reads
# --s1, --s2 = unpaired reads 1 and 2, respectively
# -k = kmer selection
rule Assemble_rRNA:
    input:
        forward=rules.Trim_galore_trim_adaptors2.output.forward_reads,
        rev=rules.Trim_galore_trim_adaptors2.output.rev_reads,
        forward_unpaired=rules.Trim_galore_trim_adaptors2.output.forward_unpaired,
        rev_unpaired=rules.Trim_galore_trim_adaptors2.output.rev_unpaired
    output:
        "10.Assemble_rRNA/{sample}/scaffolds.fasta" if ASSEMBLER == "spades" else "10.Assemble_rRNA/{sample}/final.contigs.fa"
    log: "logs/Assemble_rRNA/{sample}/{sample}_rRNA.log"
    threads: 5
    params:
        program=config['programs_path']['spades'],
        out_dir=lambda w, output: path.dirname(output[0]),
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib']
    shell:
        """
        # module load SPAdes/3.13.0
        set +u; {params.conda_activate}; set -u
        
        if [ {ASSEMBLER} == "spades" ]; then
            spades.py --careful -o {params.out_dir} \
                -1 {input.forward} -2 {input.rev} > {log} 2>&1 #\
                #--s1 {input.forward_unpaired} --s2 {input.rev_unpaired} \
                #    > {log} 2>&1
        else 
            if [ -d {params.out_dir} ]; then rm -rf {params.out_dir} ; fi 
            megahit --continue -o {params.out_dir} \
              -1 {input.forward} -2 {input.rev} > {log} 2>&1
        fi
        """

# Detect Chimeric asemblies using usearch's -uchime_ref algorithm
rule Detect_chimeras:
    input:
        database=config['FASTA'], # SILVA databases at 99% clustering
        query=rules.Assemble_rRNA.output
    output: 
        tab="11.Detect_chimera/{sample}/{sample}.chimeras.tsv",
        fasta="11.Detect_chimera/{sample}/{sample}.chimeras.fasta",
        txt= "11.Detect_chimera/{sample}/{sample}.chimeras.txt"
    log: "logs/Detect_chimera/{sample}.log"
    threads: 2
    params:
        program=config['programs_path']['usearch']
    shell:
        """ 
        {params.program} -uchime2_ref {input.query} \
                         -db {input.database} \
                         -uchimeout {output.tab} \
                         -chimeras {output.fasta} \
                         -strand plus \
                         -mode sensitive > {log} 2>&1

        # If no chimera was detected then create a dummy chimera file 
        [ -s {output.fasta} ] || echo -ne ">EMPTY\\nATTTATATATAT\\n" >  {output.fasta}

        # get the Chemic sequence name
        grep '^>' {output.fasta} | sed 's/>//g' > {output.txt}
        """
# Remove chimeras and retain only read longer Xbp e.g. 1200bp
rule Filter_seqs:
    input: 
        tab=rules.Detect_chimeras.output.txt,
        query=rules.Assemble_rRNA.output
    output: "12.Filter_seqs/{sample}/{sample}.fasta"
    log: "logs/Filter_seqs/{sample}/{sample}.log"
    threads: 1
    params:
        conda_activate=config['conda']['bioinfo']['env'],
        length=1000
    shell:
        """
        set +u; {params.conda_activate}; set -u
        seqkit grep -v -f {input.tab} {input.query} | \
         bioawk -c fastx 'length($seq) >= {params.length} {{ printf ">%s\\n%s\\n", $name,$seq }}' \
         >  {output} 2> {log}
        """

# Make blast based SILVA database
rule MakeblastDB:
    input:
        fasta=config['FASTA'],
        log_dir=rules.Make_logs_directories.output
    output: multiext(DB, ".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto")
    threads: 2
    log: "logs/MakeblastDB/MakeblastDB.log"
    params:
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib'],
        program=config['programs_path']['makeblastdb'],
        dbtype=config['parameters']['blast']['dbtype'],
        out_dir=lambda w, output: path.dirname(output[0])
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        makeblastdb \
            -in {input.fasta} \
            -out {params.out_dir}/{DBNAME} \
            -dbtype {params.dbtype} > {log} 2>&1
        """

rule BlastSeqs:
    input:
        DB=rules.MakeblastDB.output,
        query=rules.Filter_seqs.output
    output: "14.BlastSeqs/{sample}/{sample}.tsv"
    threads: 10
    log: "logs/BlastSeqs/{sample}/{sample}.log"
    params:
        conda_activate=config['conda']['bioinfo']['env'],
        PERL5LIB=config['conda']['bioinfo']['perl5lib'],
        program=config['programs_path']['blastn'],
        parallel=config['programs_path']['parallel'],
        #out_dir=lambda wildcards,input: path.dirname(input.query), 
        threads=2,
        jobs=4,
        use_parallel="False",
        blocksize="100k",
        outfmt=config['parameters']['blast']['outformat'],
        db_dir=lambda wildcards, input: path.dirname(input.DB[0])
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        cat {input.query} | \
            {params.parallel} \
              --jobs {params.jobs} \
              --block {params.blocksize} \
              --recstart '>' \
              --pipe {params.program} \
              -outfmt '{params.outfmt}' \
              -max_target_seqs 1 \
              -num_threads {params.threads} \
              -db {params.db_dir}/{DBNAME} \
              -query - > {output} 2> {log}
        """

rule SummarizeBlastResult:
    input: rules.BlastSeqs.output
    output: "15.SummarizeBlastResult/{sample}/{sample}.txt"
    threads: 1
    log: "logs/SummarizeBlastResult/{sample}/{sample}.log"
    shell:
       """
       # Percent Identity greter than or equal to 90
       # And a query coverage >= 80
       awk 'BEGIN{{FS="\t"}} $4>=98&&$15>=80{{print $3}}' {input} | \
    sort | \
    uniq -c | \
    sort -nr  > {output} 2> {log}
       """

