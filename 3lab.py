from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import AlignIO
import subprocess
from Bio.Align import AlignInfo
from Bio.Align import substitution_matrices
import io

Entrez.email = "neringamaja123@gmail.com"

def fetch_sequence(accession):
    with Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text") as handle:
        seq_record = SeqIO.read(handle, "fasta")
    return seq_record

def blast_sequence(sequence):
    result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence.seq)
    blast_records = NCBIXML.parse(result_handle)
    return blast_records

def fetch_and_align_sequences(hsa_sequence, blast_records, max_sequences, threshold):
    sequences = [hsa_sequence] 
    for blast_record in blast_records:
        for alignment in blast_record.alignments[:max_sequences]:
            for hsp in alignment.hsps:
                if hsp.expect < 0.01:
                    with Entrez.efetch(db="protein", id=alignment.accession, rettype="gb", retmode="text") as handle:
                        seq_record = SeqIO.read(handle, "genbank")
                        taxonomy = seq_record.annotations.get('taxonomy', [])
                        if 'Mammalia' in taxonomy:
                            sequences.append(seq_record)
    
    sequences = [record for record in sequences if len(record.seq) >= int(len(hsa_sequence.seq) * threshold)]
    SeqIO.write(sequences, "mammalian_sequences.fasta", "fasta")
    
    input_file = "C:\\Users\\planktonelis\\Desktop\\3 kursas\\Bioinformatika\\3lab\\Albumin_antibody_design\\mammalian_sequences.fasta"
    output_file = "C:\\Users\\planktonelis\\Desktop\\3 kursas\\Bioinformatika\\3lab\\Albumin_antibody_design\\aligned_sequences.fasta"

    executable = "C:\\Program Files\\mafft-win\\mafft.bat"

    mafft_cmd = [executable, "--auto", input_file]
    subprocess.run(mafft_cmd, stdout=open(output_file, "w"))

    alignment = AlignIO.read(output_file, "fasta")
    return alignment

def analyze_alignment(alignment, window_size):
    human_sequence = alignment[0].seq
    unique_fragment = None
    unique_fragment_score = float('inf')
    conserved_fragment = None
    conserved_fragment_score = 0  

    for i in range(len(human_sequence) - window_size + 1):
        window = human_sequence[i:i+window_size]
        window_score, match_count = calculate_window_score(window, alignment, i)

        if match_count == len(alignment) and window_score > conserved_fragment_score:
            conserved_fragment = window
            conserved_fragment_score = window_score

        if window_score < unique_fragment_score:
            unique_fragment = window
            unique_fragment_score = window_score

    return unique_fragment, conserved_fragment

def calculate_window_score(window, alignment, window_start):
    match_count = 0
    window_score = 0
    for record in alignment:
        sequence = record.seq
        exact_match = True
        for i, aa in enumerate(window):
            if aa != sequence[i + window_start] or aa == '-':
                exact_match = False
                break
        if exact_match:
            match_count += 1
            window_score += 1
    return window_score, match_count

hsa_accession = "NP_000468.1"
window_size = 15
blosum62 = substitution_matrices.load("BLOSUM62")
gap_penalty = -10

hsa_sequence = fetch_sequence(hsa_accession)
blast_records = blast_sequence(hsa_sequence)
alignment = fetch_and_align_sequences(hsa_sequence, blast_records, 6, 0.8)
unique_fragment, conserved_fragment = analyze_alignment(alignment, window_size)

print(f"Unique Fragment: {unique_fragment}")
print(f"Conserved Fragment: {conserved_fragment}")
