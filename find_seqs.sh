 mkdir Merge_reads Fastq2Fasta CountSeqs TabToFasta MakeblastDB BlastSeqs SummarizeBlastResult

find 15.SummarizeBlastResult/ -type f -size +0 -exec basename -s '.txt' {} \; | sort > turtle-seq-infile.txt

mkdir found/{sequences,blast_results,seqnames}

# Filter the blast files for samples for which turtle sequences were detected
for f in $(ls -1 01.raw_data/);do awk 'BEGIN{FS=OFS="\t"} $4>95&&$15>=80{print $1,$4,$15}' 14.BlastSeqs/${f}/${f}.tsv > found/blast_results/${f}.tsv; done && \
find found/blast_results/ -type f -size 0 -delete

cd found/
basename -s '.tsv' blast_results/* | sort > turtle-seq-infile.txt
# Write the sequence names to file for filtering
for f in $(cat turtle-seq-infile.txt);do awk 'BEGIN{FS=OFS="\t"} {print $1}' blast_results/${f}.tsv > seqnames/${f}.txt; done

# Retrieve the turtle sequences
cat turtle-seq-infile.txt | parallel --jobs 10 'seqkit grep -f seqnames/{}.txt /home/jeffbrady/Projects/Turtle/12.merge_reads/{}/{}.fasta.gz -o sequences/{}.fasta'
