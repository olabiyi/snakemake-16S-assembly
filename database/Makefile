
LSU_URL="https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz"
SSU_URL="https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
BIOAWK="/home/jeffbrady/anaconda3/envs/bioinfo/bin/bioawk"

.PHONY: clean
all: ref.fasta

LSU.fasta.gz:
	# Dowload the large subunit rRNA sequence from SILVA
	wget -O LSU.fasta.gz ${LSU_URL}

SSU.fasta.gz:
	# Dowload the small subunit rRNA sequence from SILVA
	wget -O SSU.fasta.gz ${SSU_URL}

temp.rna.fasta:
	# Combine the subunits
	zcat LSU.fasta.gz SSU.fasta.gz > temp.rna.fasta

ref.fasta: temp.rna.fasta
	# Convert RNA to DNA
	${BIOAWK} -c fastx '{print ">"$$name" "$$comment; gsub(/U/,"T",$$seq); print $$seq}' temp.rna.fasta > ref.fasta
	
clean:
	rm -rf LSU.fasta.gz SSU.fasta.gz temp.rna.fasta    
