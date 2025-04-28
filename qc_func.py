################################################################################
#	This file contains various helper functions & clases used in run_qc.py.
################################################################################
#FASTA deserialization
from Bio.SeqIO.FastaIO import SimpleFastaParser
#GC content calculation
from Bio.SeqUtils import gc_fraction
#GFF deserialization & processing
from pandas import DataFrame, read_csv
#Calling external processes: BLAST, BUSCO etc.
from subprocess import run, DEVNULL
#Dependency checks
from shutil import which
#Printing to stderr
from sys import stderr
#Proper management of paths
from pathlib import Path
#Stats
from statistics import mean, median
#Parsing of BUSCO's JSON output
import json
#Type hints
from typing import *
#Splice site distribution
from collections import Counter
#Efficient sorting
from operator import itemgetter


################################################################################
#	TYPES
################################################################################
#Simple type alias for a deserialized FASTA file
#See get_fasta() for the deserialization function
Fasta: TypeAlias = dict[str:str]

#Type alias for a deserialized GFF/GTF file
#See get_gxf() for deserialization function & gxf_*() for various operations
Gff: TypeAlias = DataFrame


################################################################################
#	MISC CONVENIENCE FUNCTIONS
################################################################################
#Terminate script with custom message
def error(msg: str):
	raise SystemExit(msg)

#Print message to stderr
#Same semantics as builtin print()
def eprint(*args, **kwargs):
	print(*args, file=stderr, **kwargs)

#Check if a given executable exists in PATH, and exit if not
def require_exec(exec: str):
	if not which(exec):
		error(f"Executable {exec} not found in PATH")

#Run a command - check the exit code, discard stderr and redirect stdout to file
#Same semantics as subprocess.run()
def exec(*args, **kwargs):
	run(*args, stderr=DEVNULL, check=True, shell=True, **kwargs)

#Returns the modification time of a file, or -Infinity if it does not exist
def get_mtime(path: Path) -> float:
	return path.stat().st_mtime if path.exists() else float("-inf")


################################################################################
#	DESERIALIZATION
################################################################################
#Deserialize FASTA file
def get_fasta(path: Path) -> Fasta:
	with open(path) as fd:
		return dict(SimpleFastaParser(fd))

#Deserialize GFF/GTF file
def get_gff(path: Path) -> Gff:
	return read_csv(path, sep='\t', header=None, comment='#')

#Deserialize a BUSCO JSON file, based on the path the results directory
#Expects only one JSON present in the directory
def get_busco(path: Path) -> dict:
	#Build path to the JSON
	json_path: Path = next((path/"busco").glob("short_summary.specific.*.json"))
	with open(json_path) as fd:
		return json.load(fd)

#Deserialize a BLAST results file in the tabular format
def get_blast(path: Path) -> DataFrame:
	return read_csv(path, sep='\t', header=None, comment='#')


################################################################################
#	EXTERNAL PROCESS CALLING
################################################################################
#Extracts all features of type `ft' from GFF file `gff',
#based on FASTA file `fasta', to FASTA file `outfile'
def run_xtractore(gff: Path, fasta: Path, ft: str, outfile: Path):
	exec(f"xtractore {gff} {fasta} -t {ft} -o {outfile}")

#Make BLAST database based on `fasta' and save it under `outfile'
def run_makeblastdb(fasta: Path, outfile: Path):
	exec(f"makeblastdb -dbtype nucl -in {fasta} -out {outfile}",
		stdout=DEVNULL)

#Run BUSCO with a given protein FASTA file as input, using a given lineage,
#saving results to a given directory
def run_busco(fasta: Path, lineage: str, outdir: Path):
	exec(f"busco -i {fasta} -l {lineage} -m prot --augustus --force -o busco" + 
		 f" --out_path {outdir.absolute()} --download_path {outdir.absolute()}",
		stdout=DEVNULL)

#Run BLAST with a preset output format, given paths to input, database & output
def run_blastn(fasta: Path, db: Path, outfile: Path):
	exec(f"blastn -query {fasta} -db {db} -out {outfile}" +
		  " -max_target_seqs 1 -evalue 1" +
		  " -outfmt '6 qseqid sseqid evalue length slen qlen pident'")

#Subset features of GFF file `in1' such that only ones which overlap with any
#features in `in2' are present, and write this file to `out'
def run_intersectBed(in1: Path, in2: Path, out: Path):
	exec(f"intersectBed -u -a {in1} -b {in2} > {out}")


################################################################################
#	GTF/GFF PROCESSING & STATISTICS
################################################################################
#Harmonize GFF & GTF formats by converting GTF-specific feature types to
#Sequence Ontology-complaint terms, and removing the attribute column
def gff_harmonize(gff: Gff):
	gff[2] = gff[2].replace("transcript", "mRNA")
	gff[2] = gff[2].replace("Selenocysteine", "selenocysteine")
	gff[8] = '.'

#Get number of features of given type present in a GTF/GFF file
def gff_ftcnt(gff: Gff, ft: str) -> int:
	return sum(gff[2] == ft)

#Get average length of features in a GFF type of a given type
def gff_ftlen(gff: Gff, ft: str) -> float:
	#Subset GFF to rows containing features of type `ft'
	subset: Gff = gff[gff[2] == ft]
	#Calculate length of each feature in `subset'
	lengths: DataFrame = subset.apply(lambda row: row[4] - row[3] + 1, axis=1)
	#Calculate mean
	return mean(lengths)

#Serialize a GFF file to a given path
def gff_save(gff: Gff, path: Path):
	gff.to_csv(path, sep='\t', header=False, index=False)


################################################################################
#	BLAST RESULTS PROCESSING & STATISTICS
################################################################################
#Get number of query sequences that have at least one hit in the results
def blast_qcnt(blast: DataFrame) -> int:
	return len(blast[0].unique())

#Get average E-value
def blast_eval(blast: DataFrame) -> float:
	return median(blast[2])

#Get average percent of identity
def blast_pident(blast: DataFrame) -> float:
	return mean(blast[6])

#Get average query coverage as percentage
def blast_qcov(blast: DataFrame) -> float:
	return mean(blast[3] / blast[5]) * 100

#Get average subject coverage as percentage
def blast_scov(blast: DataFrame) -> float:
	return mean(blast[3] / blast[4]) * 100


################################################################################
#	FASTA PROCESSING & STATISTICS
################################################################################
#Get average GC content across all sequences
def fasta_gccont(fasta: Fasta) -> float:
	return mean( gc_fraction(seq) for seq in fasta.values() ) * 100

#Get number of sequences from `fasta1' which are also present in `fasta2'
def fasta_matchcnt(fasta1: Fasta, fasta2: Fasta) -> int:
	fasta2_seqs: set[str] = set(fasta2.values())
	return sum( seq in fasta2_seqs for seq in fasta1.values() )

#Get the Jaccard index between the two sets of sequences
def fasta_jaccard(fasta1: Fasta, fasta2: Fasta) -> int:
	intersect: int = len( set(fasta1.values()) & set(fasta2.values()) )
	return intersect / ( len(fasta1) + len(fasta2) - intersect )

#Get distribution of splice sites across all sequences
#If `format' is True, returns the distribution formatted into a string; if
#False, returns a dictionary (keys are splice sites, values are their frequency)
def fasta_ssdist(fasta: Fasta, format: bool = True) -> dict[str:float]|str:
	#Count splice sites across sequences
	ss_cnts: Counter = Counter( seq[:2]+seq[-2:] for seq in fasta.values() )
	#Number of sequences
	total: int = len(fasta)
	
	#Convert absolute counts to relative frequencies
	ss_freqs: dict[str:float] = { 
			ss:(100*cnt/total) for ss, cnt in ss_cnts.items()
		}
	
	if not format:
		return freqs
	else:
		#Sort `ss_freqs' by descending frequency
		ss_freqs = sorted(ss_freqs.items(), key=itemgetter(1), reverse=True)
		#Format `ss_freqs' to a string
		return ';'.join( f"{ss}:{freq}" for ss, freq in ss_freqs )
