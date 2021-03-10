TMPDIR = "/localscratch/msauria"
SYSTEM="linux.x86_64"
MAXTHREADS=100
MAXMEM=1024
ANNOTATION_FILE="data/t2t_cenAnnotation.v2.021621.bed"
GENOMES = ["chm13v1", "GRCh38p13"]
GENOME_SIZE = {"chm13v1": "3.03e9", "GRCh38p13": "2.79e9"}
SAMPLEFILE="data/encode_samples.txt"
KMERS=["50", "75", "80", "85", "90", "95", "100"]
CHROMS=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
        'chr21', 'chr22', 'chrX']

CONCAT_DICT = {}
EXPERIMENT_DICT = {}
lib_count = {}
data = []
for line in open(SAMPLEFILE):
    if line.startswith('#'):
        continue
    line = line.rstrip().split(',')
    data.append(line)
    lib_count.setdefault(line[1], 0)
    lib_count[line[1]] += 1
    EXPERIMENT_DICT.setdefault(line[0], [])
    EXPERIMENT_DICT[line[0]].append(line[1])

for libID, count in lib_count.items():
    if count > 1:
        CONCAT_DICT[libID] = [[], []]

ARRAYS=set()
for line in open(ANNOTATION_FILE):
    a = line.split()[3].split('_')[0].upper()
    if a != "CENSAT":
        ARRAYS.add(a)
ARRAYS=list(ARRAYS)
ARRAYS.sort()

CONTROL_EXPERIMENT_DICT = {}
NONCONTROL_EXPERIMENT_DICT = {}
NONCONTROL_LIBRARIES = set()
CONTROL_LIBRARIES = set()
SLINK_DICT = {}
EXPERIMENT_CONTROL_DICT = {}
EXPERIMENT_CONTROL_LIB_DICT = {}
LABEL_DICT = {}
LABEL_CONTROL_DICT = {}
broad = set(['H3K9me3', 'H3K27me3', 'H3K36me3'])
BROAD_TRACKS = set()
NARROW_TRACKS = set()
READS = []
CTs = set()
MARK_DICT = {}
MARKS = set()

for line in data:
    expID, libID, read1, read2, control, target, ct = line[:7]
    MARKS.add(target)
    ct = ct.replace(' ', '_').replace('(', '').replace(')', '')
    CTs.add(ct)
    label = "{}_{}".format(ct, target)
    if libID in CONCAT_DICT:
        CONCAT_DICT[libID][0].append(read1)
        CONCAT_DICT[libID][1].append(read2)
    else:
        SLINK_DICT[libID] = [read1, read2]
    READS += [read1, read2]
    if control == 'NA':
        CONTROL_EXPERIMENT_DICT.setdefault(expID, []) 
        CONTROL_EXPERIMENT_DICT[expID].append(libID)
        CONTROL_LIBRARIES.add(libID)
    else:
        if target in broad:
            BROAD_TRACKS.add(label)
        else:
            NARROW_TRACKS.add(label)
        NONCONTROL_EXPERIMENT_DICT.setdefault(expID, [])
        NONCONTROL_EXPERIMENT_DICT[expID].append(libID)
        EXPERIMENT_CONTROL_LIB_DICT[expID] = EXPERIMENT_DICT[control]
        EXPERIMENT_CONTROL_DICT[expID] = control
        LABEL_CONTROL_DICT[label] = "{}_Control".format(ct)
        MARK_DICT.setdefault(ct, set())
        MARK_DICT[ct].add(target)
    LABEL_DICT[label] = expID

for key in MARK_DICT.keys():
    MARK_DICT[key] = list(MARK_DICT[key])

EXPERIMENTS = list(EXPERIMENT_DICT.keys())
LIBRARIES = list(CONCAT_DICT.keys()) + list(SLINK_DICT.keys())
CONTROL_LIBRARIES = list(CONTROL_LIBRARIES)
TRACKS=list(LABEL_DICT.keys())
TRACKS.sort()
NONCONTROL_TRACKS = []
CONTROL_TRACKS = []
for label in LABEL_DICT.keys():
    if not label.endswith("Control"):
        NONCONTROL_TRACKS.append(label)
    else:
        CONTROL_TRACKS.append(label)
NONCONTROL_TRACKS.sort()
CONTROL_TRACKS.sort()
MARKS = list(MARKS)
CTs = list(CTs)
CONTROL_EXPERIMENTS = list(CONTROL_EXPERIMENT_DICT.keys())
BROAD_TRACKS = list(BROAD_TRACKS)
NARROW_TRACKS = list(NARROW_TRACKS)

CT_MARK_SET = {}
for ct, marks in MARK_DICT.items():
    CT_MARK_SET[ct] = "|".join(marks)

CONTROL_EXPERIMENT_SET = "|".join(CONTROL_EXPERIMENTS)
NONCONTROL_EXPERIMENT_SET = "|".join(list(NONCONTROL_EXPERIMENT_DICT.keys()))
EXPERIMENT_SET = "|".join(EXPERIMENTS)
TRACK_SET = "|".join(TRACKS)
NONCONTROL_TRACK_SET = "|".join(NONCONTROL_TRACKS)
LIBRARY_SET = "|".join(LIBRARIES)
CONTROL_LIBRARY_SET = "|".join(CONTROL_LIBRARIES)
CONCAT_SET = "|".join(list(CONCAT_DICT.keys()))
SLINK_SET = "|".join(list(SLINK_DICT.keys()))
BROAD_TRACK_SET = "|".join(BROAD_TRACKS)
NARROW_TRACK_SET = "|".join(NARROW_TRACKS)
GENOME_SET = "|".join(GENOMES)
CHROM_SET="|".join(CHROMS)
KMER_SET="|".join(KMERS)
ARRAY_SET = "|".join(ARRAYS)
MARK_SET = "|".join(MARKS)
CT_SET = "|".join(CTs)
READ_SET = "|".join(READS)


rule all:
    input:
        expand("results/{file}",
               file=["chm13v1_array_enrichments.tsv",
                     "macs2_peak_overlap.txt",
                     "control_mapping_stats.txt",
                     "mapping_stats.txt"]),
        expand("Kmer_BigWigs/{genome}_{kmer}mer.bw",
               kmer=KMERS, genome=GENOMES)


######## Get software #########

rule download_software:
    output:
        "bin/{sw}"
    wildcard_constraints:
        sw="liftOver|bedToBigBed|wigToBigWig"
    params:
        sw="{sw}",
        system=SYSTEM
    log:
        "logs/Encode/download/{sw}.log"
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/admin/exe/{params.system}/{params.sw} -O {output}
        """

rule build_KMC:
    output:
        "bin/kmc",
        "bin/kmc_genome_counts"
    log:
        "logs/build_kmc.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        git clone https://github.com/msauria/KMC.git
        git checkout kmer_mapping
        cd KMC && make
        cp KMC/bin/kmc bin/
        cp KMC/bin/kmc_genome_counts bin/
        """


######## Get chm13 data #########

rule preprocess_annotations:
    output:
        "data/t2t_cenSatAnnotation.bed",
        "data/t2t_cenSatAnnotation2.bed"
    params:
        anno=ANNOTATION_FILE
    log:
        "logs/preprocess/cenSatAnnotation.log"
    shell:
        """
        grep -v "censat" {params.anno} > {output[1]}
        awk 'BEGIN{{ OFS="\t" }} {{ split($4,a,"_"); $4=toupper(a[1]); print $0 }}' {output[1]} > {output[0]}
        """

rule download_liftover:
    output:
        "data/hg38.t2t-chm13-v1.0.over.chain.gz"
    log:
        "logs/download/liftover.log"
    shell:
        """
        wget -O {output} http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/hg38Lastz/hg38.t2t-chm13-v1.0.over.chain.gz
        """

######## Get genomic data #########

rule get_GRCh38_fasta:
    output:
        fa="fasta/GRCh38p13.fasta",
        sizes="fasta/GRCh38p13.chrom.sizes"
    log:
        "logs/fasta/GRCh38p13_download.log"
    shell:
        """
        wget -O {output.fa}.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gzip -d -c {output.fa}.gz > {output.fa}
        wget -O {output.sizes} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
        """

rule get_chm13v1_fasta:
    output:
        fa="fasta/chm13v1.fasta",
        sizes="fasta/chm13v1.chrom.sizes"
    log:
        "logs/fasta/chm13v1_download.log"
    shell:
        """
        wget -O {output.fa}.gz http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/genome/t2t-chm13-v1.0.fa.gz
        gzip -d -c {output.fa}.gz > {output.fa}
        cat {output.fa} | awk 'BEGIN{{ OFS=""; ORS=""; tab="\t"; nl="\n"; tl=0 }}\
                               {{ if ( NR == 1 ){{ split($1,A,">"); print A[2],tab; }} \
                                  else {{ if ( $1 ~ /^>/ ){{ split($1,A,">"); print tl,nl,A[2],tab; tl = 0;}} \
                                          else {{ tl = tl + length; }} }} }}\
                               END{{ print tl,nl }}' | sort -r -k2,2n > {output.sizes}
        """


######## Get Encode data #########

rule download_data:
    output:
        "fastq/{readID}.fastq.gz"
    params:
        readID="{readID}"
    wildcard_constraints:
        readID=READ_SET
    log:
        "logs/download/{readID}.log"
    shell:
        """
        wget https://www.encodeproject.org/files/{params.readID}/@@download/{params.readID}.fastq.gz -O {output}
        """

######## Concatenate fastq files #########

rule concat:
    input:
        lambda wildcards: expand("fastq/{read}.fastq.gz",
                                 read=CONCAT_DICT[wildcards.libID][int(wildcards.side)-1])
    output:
        "concat/{libID}.{side}.fastq.gz"
    wildcard_constraints:
        libID=CONCAT_SET,
        genome=GENOME_SET,
        side="1|2"
    log:
        "logs/concat/{libID}.{side}.log"
    shell:
        """
        zcat {input} | gzip > {output}
        """

rule symlink:
    input:
        lambda wildcards: expand("fastq/{read}.fastq.gz",
                                 read=SLINK_DICT[wildcards.libID])
    output:
        "concat/{libID}.1.fastq.gz",
        "concat/{libID}.2.fastq.gz"
    wildcard_constraints:
        libID=SLINK_SET,
        genome=GENOME_SET
    log:
        "logs/symlink/{libID}.log"
    shell:
        """
        ln -s ${{PWD}}/{input[0]} {output[0]}
        ln -s ${{PWD}}/{input[1]} {output[1]}
        """

rule count_fastq:
    input:
        "concat/{libID}.1.fastq.gz"
    output:
        "concat/{libID}.fastq.count"
    wildcard_constraints:
        libID=LIBRARY_SET
    log:
        "logs/count/{libID}.log"
    shell:
        """
        echo "$(expr $(eval zcat {input} | wc -l) / 4)" > {output}
        """


######## Bowtie2 Alignment #########

rule bowtie2_build:
    input:
        "fasta/{genome}.fasta"
    output:
        multiext(
            "bt2_indices/{genome}",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
        )
    params:
        extra=""  # optional parameters
    threads: 8
    log:
        "logs/bowtie2_build.{genome}.log"
    wrapper:
        "v0.69.0/bio/bowtie2/build"

rule bowtie2:
    input:
        sample=["concat/{libID}.1.fastq.gz", "concat/{libID}.2.fastq.gz"],
        index="bt2_indices/{genome}.1.bt2"
    output:
        "mapped/{libID}.{genome}.bam"
    params:
        index="bt2_indices/{genome}",  # prefix of reference genome index (built with bowtie2-build)
        extra="--no-discordant --no-mixed --very-sensitive \
               --no-unal --omit-sec-seq --xeq --reorder"  # optional parameters
    threads: 100  # Use at least two threads
    wildcard_constraints:
        genome=GENOME_SET,
        libID=LIBRARY_SET
    log:
        "logs/bowtie2/{libID}.{genome}.log"
    wrapper:
        "v0.69.0/bio/bowtie2/align"


######## Create KMC DBs #########

rule create_kmc_db:
    input:
        fa="fasta/{genome}.fasta",
        kmc="bin/kmc"
    output:
        "KMC_db/{genome}_{kmer}.kmc_pre",
        "KMC_db/{genome}_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="KMC_db/{genome}_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=KMER_SET,
        genome=GENOME_SET
    log:
        "logs/kmc/{genome}_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci2 -t{threads} \
            {input.fa} {params.prefix} {params.tmpdir}
        """


######## Create Kmer BigWigs #########

rule create_kmer_wiggle:
    input:
        fa="fasta/{genome}.fasta",
        db="KMC_db/{genome}_{kmer}.kmc_pre",
        kmc="bin/kmc_genome_counts"
    output:
        temp("Kmer_BigWigs/{genome}_{kmer}mer.wig")
    params:
        prefix="KMC_db/{genome}_{kmer}",
        chroms=CHROMS,
        tmpdir=TMPDIR
    wildcard_constraints:
        genome=GENOME_SET,
        kmer=KMER_SET
    threads:
        len(CHROMS)
    log:
        "logs/wiggle/{genome}_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        tmpdir=$(mktemp -d -p {params.tmpdir});
        for C in {params.chroms}; do \
            {input.kmc} -ch$C {params.prefix} {input.fa} > $tmpdir/$C.wig & \
        done; wait
        wigs=""
        for C in {params.chroms}; do wigs="$wigs $tmpdir/$C.wig"; done;
        cat $wigs > {output}
        rm -r $tmpdir
        """

rule create_kmer_bigwig:
    input:
        wig="Kmer_BigWigs/{genome}_{kmer}mer.wig",
        wigtobigwig="bin/wigToBigWig",
        sizes="fasta/{genome}.chrom.sizes"
    output:
        "Kmer_BigWigs/{genome}_{kmer}mer.bw"
    wildcard_constraints:
        genome=GENOME_SET,
        kmer=KMER_SET
    log:
        "logs/bigwigs/kmer_{genome}_{kmer}.log"
    shell:
        """
        {input.wigtobigwig} -clip {input.wig} {input.sizes} {output}
        """


######## BAM filtering #########

rule samtools_flag_mapped:
    input:
        "mapped/{libID}.{genome}.bam"
    output:
        "mapped/{libID}.{genome}.stats"
    wildcard_constraints:
        genome=GENOME_SET,
        libID=LIBRARY_SET
    log:
        "logs/samtools/stats/mapped.{libID}.{genome}.log"
    wrapper:
        "0.68.0/bio/samtools/flagstat"

rule samtools_flag_exp:
    input:
        "{folder}/{expID}.{genome}.bam"
    output:
        "{folder}/{expID}.{genome}.stats"
    wildcard_constraints:
        genome=GENOME_SET,
        expID=EXPERIMENT_SET
    log:
        "logs/samtools/stats/{folder}.{expID}.{genome}.log"
    wrapper:
        "0.68.0/bio/samtools/flagstat"

rule samtools_index:
    input:
        "{folder}/{expID}.{genome}.bam"
    output:
        "{folder}/{expID}.{genome}.bam.bai"
    wildcard_constraints:
        genome=GENOME_SET,
        expID=EXPERIMENT_SET
    log:
        "logs/samtools/index/{folder}.{expID}.{genome}.log"
    wrapper:
        "0.68.0/bio/samtools/index"

rule samtools_filter:
    input:
        "mapped/{libID}.{genome}.bam"
    output:
        "filtered/{libID}.{genome}.bam",
        "filtered/{libID}.{genome}.stats"
    wildcard_constraints:
        genome=GENOME_SET,
        libID=LIBRARY_SET
    log:
        "logs/samtools/stats/mapped.{libID}.{genome}.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -b -F 1804 -f 2 -q 2 {input[0]} > {output[0]}
        samtools flagstats {output[0]} > {output[1]}
        """
# 1804 = 1024+512+256+8+4
# Remove flag bits 4, 8, 256, 512, and 1024 (unmapped segment, unmapped next segment, secondary, failed filter, duplicate)
# Require flag bit 2 (both sides properly aligned)

rule remove_duplicates:
    input:
        "filtered/{libID}.{genome}.bam"
    output:
        bam="dedup/{libID}.{genome}.bam",
        metrics="dedup/{libID}.{genome}.metrics.txt"
    params:
        "VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=queryname REMOVE_DUPLICATES=true"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    wildcard_constraints:
        libID=LIBRARY_SET,
        genome=GENOME_SET
    resources:
        mem_mb=MAXMEM
    log:
        "logs/picard/dedup/{libID}.{genome}.log"
    wrapper:
        "0.68.0/bio/picard/markduplicates"

rule concat_bam:
    input:
        lambda wildcards: expand("{{folder}}/{libID}.{{genome}}.bam",
                                 libID=EXPERIMENT_DICT[wildcards.expID])
    output:
        "{folder}_concat/{expID}.{genome}.bam"
    params:
        prefix="{folder}_concat/{expID}.{genome}"
    wildcard_constraints:
        expID=EXPERIMENT_SET,
        genome=GENOME_SET,
        folder="dedup|filtered|mapped"
    threads:
        8
    log:
        "logs/samtools/merge/{folder}.{expID}.{genome}.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        INPUTS=({input})
        if [ ${{#INPUTS[*]}} \> 1 ]; then \
            samtools merge -n --threads {threads} {params.prefix}.tmp {input}
        else \
            cp {input[0]} {params.prefix}.tmp
        fi;
        samtools sort -@ {threads} -o {output} {params.prefix}.tmp
        rm {params.prefix}.tmp
        """

rule filter_by_kmers_presort:
    input:
        "dedup_concat/{expID}.{genome}.bam"
    output:
       temp("dedup_kmer/{expID}.{genome}.namesorted.bam")
    wildcard_constraints:
        genome=GENOME_SET,
        expID=EXPERIMENT_SET
    params:
        tmp_dir=TMPDIR,
        extra="-n"
    threads: 8
    log:
        "logs/samtools/sort/dedup_concat.{expID}.{genome}.log"
    wrapper:
        "0.68.0/bio/samtools/sort"

rule filter_by_kmers:
    input:
        bam="dedup_kmer/{expID}.{genome}.namesorted.bam",
        bigwigs=expand("Kmer_BigWigs/{{genome}}_{kmer}mer.bw", kmer=KMERS)
    output:
        temp("dedup_kmer/{expID}.{genome}.filtered.bam")
    params:
        kmers=KMERS,
        genome="{genome}",
        prefix="{expID}.{genome}",
        kmer_prefix="Kmer_BigWigs/{genome}",
        tmpdir=TMPDIR
    wildcard_constraints:
        genome=GENOME_SET,
        expID=EXPERIMENT_SET
    threads: 5
    log:
        "logs/kmer_filter/{expID}.{genome}.log"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        tmpdir=$(mktemp -d -p {params.tmpdir})
        kmers=({params.kmers})
        SIZES=${{kmers[0]}}
        cp {params.kmer_prefix}_${{kmers[0]}}mer.bw $tmpdir
        for I in $(seq 1 $(expr ${{#kmers[*]}} - 1)); do
            SIZES="${{SIZES}},${{kmers[$I]}}"
            cp {params.kmer_prefix}_${{kmers[$I]}}mer.bw $tmpdir
        done
        mv {params.tmpdir}/{params.prefix}.namesorted.bam ${{tmpdir}}/
        cp {input} ${{tmpdir}}/{params.prefix}.namesorted.bam
        bin/filter_by_unique_kmers.py {input} "${{tmpdir}}/{params.genome}_*mer.bw" \
            $SIZES ${{tmpdir}}/{params.prefix}.tmp
        mv ${{tmpdir}}/{params.prefix}.tmp {output}
        rm -r $tmpdir
        """

rule filter_by_kmers_postsort:
    input:
        "dedup_kmer/{expID}.{genome}.filtered.bam"
    output:
       "dedup_kmer/{expID}.{genome}.bam"
    wildcard_constraints:
        genome=GENOME_SET,
        expID=EXPERIMENT_SET
    params:
        tmp_dir=TMPDIR,
        extra=""
    threads: 8
    log:
        "logs/samtools/sort/dedup_kmer.{expID}.{genome}.log"
    wrapper:
        "0.68.0/bio/samtools/sort"


######## Encode BigWigs #########

rule bam_coverage:
    input:
        bam=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam",
                                     expID=LABEL_DICT[wildcards.label]),
        index=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam.bai",
                                       expID=LABEL_DICT[wildcards.label])
    output:
        "Encode_BigWigs/{label}.{genome}.bw"
    wildcard_constraints:
        label=TRACK_SET,
        genome=GENOME_SET
    threads:
        16
    log:
        "logs/deeptools/bamCoverage/{label}.{genome}.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output} -p {threads} -bs 1
        """

rule bigwig_compare:
    input:
        treat="Encode_BigWigs/{label}.{genome}.bw",
        control=lambda wildcards: expand("bigwigs/{control}.{{genome}}.bw",
                                         control=LABEL_CONTROL_DICT[wildcards.label])
    output:
        "Encode_BigWigs/{label}.{genome}_normed.bw"
    wildcard_constraints:
        label=NONCONTROL_TRACK_SET,
        genome=GENOME_SET
    threads:
        16
    log:
        "logs/deeptools/bigwigCompare/{label}.{genome}.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        bigwigCompare -p {threads} -bs 50 -b1 {input.treat} -b2 {input.control} \
            -o {output} --skipZeroOverZero -of bigwig
        """


######## Macs2 #########

rule macs2_callpeak:
    input:
        treatment=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam",
                                           expID=LABEL_DICT[wildcards.label]),
        control=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam",
                                         expID=EXPERIMENT_CONTROL_DICT[LABEL_DICT[wildcards.label]]),
        treatment_index=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam.bai",
                                                 expID=LABEL_DICT[wildcards.label]),
        control_index=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam.bai",
                                               expID=EXPERIMENT_CONTROL_DICT[LABEL_DICT[wildcards.label]])
    output:
        multiext("macs2/{label}.{genome}",
                 "_peaks.xls",
                 "_peaks.narrowPeak")
    params:
        lambda wildcards: expand("-f BAMPE -g {size}", size=GENOME_SIZE[wildcards.genome])[0]
    wildcard_constraints:
        label=NARROW_TRACK_SET,
        genome=GENOME_SET
    log:
        "logs/macs2/{label}.{genome}.log"
    wrapper:
        "0.66.0/bio/macs2/callpeak"

rule macs2_callpeak_broad:
    input:
        treatment=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam",
                                           expID=LABEL_DICT[wildcards.label]),
        control=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam",
                                         expID=EXPERIMENT_CONTROL_DICT[LABEL_DICT[wildcards.label]]),
        treatment_index=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam.bai",
                                                 expID=LABEL_DICT[wildcards.label]),
        control_index=lambda wildcards: expand("dedup_kmer/{expID}.{{genome}}.bam.bai",
                                               expID=EXPERIMENT_CONTROL_DICT[LABEL_DICT[wildcards.label]])
    output:
        multiext("macs2/{label}.{genome}",
                 "_peaks.xls",
                 "_peaks.broadPeak")
    params:
        lambda wildcards: expand("-f BAMPE -g {size}", size=GENOME_SIZE[wildcards.genome])[0]
    wildcard_constraints:
        label=BROAD_TRACK_SET,
        genome=GENOME_SET
    log:
        "logs/macs2/{label}.{genome}.log"
    wrapper:
        "0.66.0/bio/macs2/callpeak"

rule bedToBigBed:
    input:
        peaks="macs2/{label}.{genome}_peaks.narrowPeak",
        sizes="fasta/{genome}.chrom.sizes",
        bedtobigbed="bin/bedToBigBed"
    output:
        "Encode_BigBeds/{label}.{genome}_peaks.bb"
    params:
        tmp="Encode_BigBeds/{label}.{genome}.tmp"
    wildcard_constraints:
        label=NARROW_TRACK_SET,
        genome=GENOME_SET
    log:
        "logs/bigbeds/{label}.{genome}.log"
    shell:
        """
        awk -v OFS="\t" '{{$5=$5>1000?1000:$5}} {{print}}' {input.peaks} > {params.tmp};
        {input.bedtobigbed} -type=bed6+4 {params.tmp} {input.sizes} {output};
        rm -f {params.tmp}
        """

rule bedToBigBed_broad:
    input:
        peaks="macs2/{label}.{genome}_peaks.broadPeak",
        sizes="fasta/{genome}.chrom.sizes",
        bedtobigbed="bin/bedToBigBed"
    output:
        "Encode_BigBeds/{label}.{genome}_peaks.bb"
    params:
        tmp="Encode_BigBeds/{label}.{genome}.tmp"
    wildcard_constraints:
        label=BROAD_TRACK_SET,
        genome=GENOME_SET
    log:
        "logs/bigbeds_{label}.{genome}.log"
    shell:
        """
        awk -v OFS="\t" '{{$5=$5>1000?1000:$5}} {{print}}' {input.peaks} > {params.tmp};
        {input.bedtobigbed} -type=bed6+3 {params.tmp} {input.sizes} {output};
        rm -f {params.tmp}
        """

rule liftOver_peaks:
    input:
        peaks="macs2/{track}.GRCh38p13_peaks.narrowPeak",
        liftover="bin/liftOver",
        bedtobigbed="bin/bedToBigBed",
        sizes="fasta/chm13v1.chrom.sizes",
        chain="data/hg38.t2t-chm13-v1.0.over.chain.gz"
    output:
        lifted="macs2/{track}.GRCh38p13-chm13v1_peaks.narrowPeak",
        unlifted="macs2/{track}.GRCh38p13-chm13v1_peaks.unlifted",
        bb="Encode_BigBeds/{track}_GRCh38p13LO_peaks.bb"
    params:
        track="{track}"
    wildcard_constraints:
        track=NARROW_TRACK_SET
    log:
        "logs/liftover/{track}.log"
    shell:
        """
        {input.liftover} -minMatch=0.2 -bedPlus=3 {input.peaks} {input.chain} \
            tmp_{params.track} {output.unlifted}
        sort -k1,1 -k2,2n tmp_{params.track} > tmp2_{params.track}
        awk 'BEGIN{{ OFS="\t" }} NR==FNR {{ sizes[$1]=$2; next }} \
             {{ if (sizes[$1] > $3) print $0; else if (sizes[$1] > $2){{ $3=sizes[$1]; print $0; }} }}' \
             {input.sizes} tmp2_{params.track} > {output.lifted}
        {input.bedtobigbed} -type=bed3+6 {output.lifted} {input.sizes} {output.bb}
        rm tmp_{params.track} tmp2_{params.track}
        """

rule liftOver_peaks_broad:
    input:
        peaks="macs2/{track}.GRCh38p13_peaks.broadPeak",
        liftover="bin/liftOver",
        bedtobigbed="bin/bedToBigBed",
        sizes="fasta/chm13v1.chrom.sizes",
        chain="data/hg38.t2t-chm13-v1.0.over.chain.gz"
    output:
        lifted="macs2/{track}.GRCh38p13-chm13v1_peaks.broadPeak",
        unlifted="macs2/{track}.GRCh38p13-chm13v1_peaks.unlifted",
        bb="Encode_BigBeds/{track}_GRCh38p13LO_peaks.bb",
    params:
        track="{track}"
    wildcard_constraints:
        track=BROAD_TRACK_SET
    log:
        "logs/liftover/{track}.log"
    shell:
        """
        {input.liftover} -minMatch=0.2 -bedPlus=3 {input.peaks} {input.chain} \
            tmp_{params.track} {output.unlifted}
        sort -k1,1 -k2,2n tmp_{params.track} > tmp2_{params.track}
        awk 'BEGIN{{ OFS="\t" }} NR==FNR {{ sizes[$1]=$2; next }} \
             {{ if (sizes[$1] > $3) print $0; else if (sizes[$1] > $2){{ $3=sizes[$1]; print $0; }} }}' \
             {input.sizes} tmp2_{params.track} > {output.lifted}
        {input.bedtobigbed} -type=bed3+6 {output.lifted} {input.sizes} {output.bb}
        rm tmp_{params.track} tmp2_{params.track}
        """


######## Analysis #########

rule get_control_mapping_stats:
    input:
        expand("{folder}/{libID}.{genome}.stats",
               folder=['mapped','filtered'],
               genome=GENOMES, libID=CONTROL_LIBRARIES),
        expand("dedup_concat/{expID}.{genome}.stats",
               genome=GENOMES, expID=CONTROL_EXPERIMENTS),
        expand("dedup_kmer/{expID}.{genome}.stats",
               genome=GENOMES, expID=CONTROL_EXPERIMENTS)
    output:
        "results/control_mapping_stats.txt"
    params:
        genomes=GENOMES,
        folders=['mapped','filtered','dedup', 'kmer'],
        expID=CONTROL_EXPERIMENTS,
        samples=SAMPLEFILE
    log:
        "logs/analysis/control_mapping_stats.log"
    shell:
        """
        genomes=({params.genomes})
        folders=({params.folders})
        expID=({params.expID})
        temp="expID,group,mark"
        for G in ${{genomes[*]}}; do
            temp="${{temp}},${{G}}"
        done
        echo -e "$temp" > {output}
        for E in ${{expID[*]}}; do
            total=$(eval grep ${{E}} {params.samples} | grep -v ",${{E}}," | \
                    awk 'BEGIN{{ FS=","; ORS=""; COUNT=0 }}{{ COUNT=COUNT+$9 }}END{{ print COUNT }}')
            echo "${{E}},total,${{total}},${{total}}" >> {output}
            LIBS=($(eval grep ${{E}} {params.samples} | grep -v ",${{E}}," | \
                    awk 'BEGIN{{ FS=","; OFS=""; ORS=" " }}{{ print $2 }}'))
            for F in mapped filtered; do
                temp="${{E}},${{F}}
                for G in ${{genomes[*]}}; do
                    C=0
                    for L in ${{LIBS[*]}}; do
                        let C=C+$(eval tail -n +5 ${{F}}/${{L}}.${{G}}.stats | \
                                  head -n 1 | awk 'BEGIN{{ ORS="" }}{{ print $1 }}')/2
                    done
                    temp="${{temp}},$C"
                done
                echo $temp >> output
            done
            temp="${{E}},${{F}}
            for G in ${{genomes[*]}}; do
                C=$(eval tail -n +5 dedup_concat/${{E}}.${{G}}.stats | \
                    head -n 1 | awk 'BEGIN{{ ORS="" }}{{ print $1 }}')
                let C=C/2
                temp="${{temp}},$C"
            done
            echo $temp >> {output}
            temp="${{E}},${{F}}
            for G in ${{genomes[*]}}; do
                C=$(eval tail -n +5 dedup_kmer/${{E}}.${{G}}.stats | \
                    head -n 1 | awk 'BEGIN{{ ORS="" }}{{ print $1 }}')
                let C=C/2
                temp="${{temp}},$C"
            done
            echo $temp >> {output}
        done
        """

rule get_mapping_stats:
    input:
        expand("{folder}/{libID}.{genome}.stats",
               folder=['mapped','filtered'],
               genome=GENOMES, libID=LIBRARIES),
        expand("dedup_concat/{expID}.{genome}.stats",
               genome=GENOMES, expID=EXPERIMENTS),
        expand("dedup_kmer/{expID}.{genome}.stats",
               genome=GENOMES, expID=EXPERIMENTS)
    output:
        "results/mapping_stats.txt"
    params:
        genomes=GENOMES,
        folders=['mapped','filtered','dedup', 'kmer'],
        expID=EXPERIMENTS,
        samples=SAMPLEFILE
    log:
        "logs/analysis/mapping_stats.log"
    shell:
        """
        genomes=({params.genomes})
        folders=({params.folders})
        expID=({params.expID})
        temp=",,,"
        for G in ${{genomes[*]}}; do
            for F in ${{folders[*]}}; do
                temp="${{temp}},${{G}}"
            done
        done
        echo -e "$temp" > {output}
        temp="Experiment,Cell type,Target,Total Reads"
        for G in ${{genomes[*]}}; do
            for F in ${{folders[*]}}; do
                temp="${{temp}},${{F}}"
            done
        done
        echo -e "$temp" >> {output}
        for E in ${{expID[*]}}; do
            temp=$(eval grep ${{E}} {params.samples} | grep -v ",${{E}}," | \
                   head -n 1 | awk 'BEGIN{{ FS=","; OFS=","; ORS=""; }}{{ print $1,$7,$6 }}')
            temp="${{temp}},$(eval grep ${{E}} {params.samples} | \
                              grep -v ",${{E}}," | \
                              awk 'BEGIN{{ FS=","; ORS=""; COUNT=0 }}{{ COUNT=COUNT+$9 }}END{{ print COUNT }}')"
            LIBS=($(eval grep ${{E}} {params.samples} | grep -v ",${{E}}," | \
                    awk 'BEGIN{{ FS=","; OFS=""; ORS=" " }}{{ print $2 }}'))
            for G in ${{genomes[*]}}; do
                for F in mapped filtered; do
                    C=0
                    for L in ${{LIBS[*]}}; do
                        let C=C+$(eval tail -n +5 ${{F}}/${{L}}.${{G}}.stats | \
                                  head -n 1 | awk 'BEGIN{{ ORS="" }}{{ print $1 }}')/2
                    done
                    temp="${{temp}},$C"
                done
                C=$(eval tail -n +5 dedup_concat/${{E}}.${{G}}.stats | \
                    head -n 1 | awk 'BEGIN{{ ORS="" }}{{ print $1 }}')
                let C=C/2
                temp="${{temp}},$C"
                C=$(eval tail -n +5 dedup_kmer/${{E}}.${{G}}.stats | \
                    head -n 1 | awk 'BEGIN{{ ORS="" }}{{ print $1 }}')
                let C=C/2
                temp="${{temp}},$C"
            done
            echo -e $temp >> {output}
        done
        """

rule count_reads:
    input:
        bam="dedup_kmer/{expID}.{genome}.bam",
        index="dedup_kmer/{expID}.{genome}.bam.bai",
        anno="data/t2t_cenSatAnnotation.bed"
    output:
        "region_scores/{expID}.{genome}_counts.npy"
    wildcard_constraints:
        expID=EXPERIMENT_SET,
        genome=GENOME_SET
    log:
        "logs/count_reads/{expID}.{genome}.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        bin/count_region_reads.py {input.bam} {input.anno} {output}
        """

rule score_read_counts:
    input:
        counts=expand("region_scores/{expID}.{{genome}}_counts.npy",
                      expID=EXPERIMENTS),
        stats=expand("dedup_kmer/{expID}.{{genome}}.stats",
                      expID=EXPERIMENTS),
        anno="data/t2t_cenSatAnnotation2.bed"
    output:
        "results/{genome}_array_enrichments.tsv"
    params:
        samples=SAMPLEFILE,
        datadir="region_scores"
    wildcard_constraints:
        genome=GENOME_SET
    log:
        "logs/analysis/{genome}_array_enrichments.log"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        bin/score_read_counts.py {params.samples} {input.anno} {params.datadir} {output}
        """

rule liftOver_intersection:
    input:
        c_narrow=expand("macs2/{label}.chm13v1_peaks.narrowPeak",
                        label=NARROW_TRACKS),
        c_broad=expand("macs2/{label}.chm13v1_peaks.broadPeak",
                       label=BROAD_TRACKS),
        G_narrow=expand("macs2/{label}.GRCh38p13_peaks.narrowPeak",
                        label=NARROW_TRACKS),
        G_broad=expand("macs2/{label}.GRCh38p13_peaks.broadPeak",
                       label=BROAD_TRACKS),
        GLO_narrow=expand("macs2/{label}.GRCh38p13-chm13v1_peaks.narrowPeak",
                          label=NARROW_TRACKS),
        GLO_broad=expand("macs2/{label}.GRCh38p13-chm13v1_peaks.broadPeak",
                         label=BROAD_TRACKS)
    output:
        "results/macs2_peak_overlap.txt"
    params:
        broad=BROAD_TRACKS,
        narrow=NARROW_TRACKS
    log:
        "logs/analysis/liftOver_intersection.log"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        BROAD=({params.broad})
        NARROW=({params.narrow})
        echo -e "track\tchm13v1\tGRCh38p13\tlifted\toverlap" > {output}
        for mark in ${{BROAD[*]}}; do
            G="macs2/${{mark}}.GRCh38p13_peaks.broadPeak"
            C="macs2/${{mark}}.chm13v1_peaks.broadPeak"
            L="macs2/${{mark}}.GRCh38p13-chm13v1_peaks.broadPeak"
            bedtools intersect -u -a ${{L}} -b ${{C}} > intersected
            echo -e "${{mark}}\t$(eval wc -l ${{C}} | \
                awk '{{ print $1 }}')\t$(eval wc -l ${{G}} | \
                awk '{{ print $1 }}')\t$(wc -l ${{L}} | \
                awk '{{ print $1 }}')\t$(eval wc -l intersected | \
                awk '{{ print $1 }}')" >> {output}
        done
        for mark in ${{NARROW[*]}}; do
            G="macs2/${{mark}}.GRCh38p13_peaks.narrowPeak"
            C="macs2/${{mark}}.chm13v1_peaks.narrowPeak"
            L="macs2/${{mark}}.GRCh38p13-chm13v1_peaks.narrowPeak"
            bedtools intersect -u -a ${{L}} -b ${{C}} > intersected
            echo -e "${{mark}}\t$(eval wc -l ${{C}} | \
            awk '{{ print $1 }}')\t$(eval wc -l ${{G}} | \
            awk '{{ print $1 }}')\t$(wc -l ${{L}} | \
            awk '{{ print $1 }}')\t$(eval wc -l intersected | \
            awk '{{ print $1 }}')" >> {output}
        done
        rm intersected
        """
