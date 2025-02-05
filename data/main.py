from Bio import SeqIO
import pandas as pd
import psutil

# --- Memory Logging Function ---
def log_memory_usage(step):
    mem = psutil.virtual_memory()
    print(f"[{step}] Memory usage: {mem.percent}% used of {mem.total / (1024**3):.2f} GB")

# --- Step 1: Parse FASTA File ---
fasta_file = "data/sequence.fa"
transcripts = []

for record in SeqIO.parse(fasta_file, "fasta"):
    header_parts = record.description.split("|")  # Split header fields
    transcript_id = header_parts[0]  # ENST00000641515.2
    gene_id = header_parts[1]  # ENSG00000186092.7
    gene_name = header_parts[5]  # OR4F5

    # Extract CDS and UTR regions
    utr5 = cds = utr3 = None
    for part in header_parts:
        if part.startswith("UTR5:"):
            utr5 = part.split(":")[1]  # Get 5' UTR positions
        elif part.startswith("CDS:"):
            cds = part.split(":")[1]  # Get CDS positions
        elif part.startswith("UTR3:"):
            utr3 = part.split(":")[1]  # Get 3' UTR positions

    sequence = str(record.seq)  # Full transcript sequence
    transcripts.append([transcript_id, gene_id, gene_name, utr5, cds, utr3, sequence])

# Save parsed transcripts to CSV
transcripts_df = pd.DataFrame(transcripts, columns=["Transcript_ID", "Gene_ID", "Gene_Name", "UTR5", "CDS", "UTR3", "Sequence"])
transcripts_df.to_csv("data/parsed_transcripts.csv", index=False)
print("✅ Transcript sequences saved to parsed_transcripts.csv!")
log_memory_usage("FASTA parsing completed")

# --- Step 2: Chunked Junction Processing ---
junction_file = "data/GTEx_Analysis_v10_STARv2.7.10a_junctions.gct"

# Function to check if a junction maps to an exon
def map_junction_to_exon(junction, cds_start, cds_end):
    return (junction["Start"] >= cds_start) and (junction["End"] <= cds_end)

# Process junctions in chunks
def process_junction_chunks(junction_file, transcripts_df, chunk_size=100000):
    psi_values = []
    splicing_labels = []

    for index, row in transcripts_df.iterrows():
        transcript_id = row["Transcript_ID"]
        cds = row["CDS"]

        if pd.isna(cds) or cds == "":
            print(f"No CDS found for Transcript ID: {transcript_id}")
            psi_values.append(None)
            splicing_labels.append(None)
            continue

        cds_ranges = cds.split(",")
        inclusion_reads, skipping_reads = 0, 0

        print(f"Processing Transcript {transcript_id}...")
        log_memory_usage(f"Before processing {transcript_id}")

        for chunk in pd.read_csv(junction_file, sep="\t", skiprows=2, chunksize=chunk_size):
            chunk = chunk[["Chromosome", "Start", "End", "Strand", "ReadCount"]]  # Keep relevant columns
            print(f"Processing chunk with {len(chunk)} junctions...")

            for range_str in cds_ranges:
                cds_start, cds_end = map(int, range_str.split("-"))
                print(f"Processing CDS range: {cds_start}-{cds_end}")

            # Calculate inclusion reads
            inclusion_reads_chunk = chunk[
                chunk.apply(lambda x: map_junction_to_exon(x, cds_start, cds_end), axis=1)
            ]["ReadCount"].sum()
            inclusion_reads += inclusion_reads_chunk
            print(f"Inclusion reads for this chunk: {inclusion_reads_chunk}")

            # Calculate skipping reads
            skipping_reads_chunk = chunk[
                (chunk["Start"] < cds_start) & (chunk["End"] > cds_end)
            ]["ReadCount"].sum()
            skipping_reads += skipping_reads_chunk
            print(f"Skipping reads for this chunk: {skipping_reads_chunk}")

        psi = inclusion_reads / (inclusion_reads + skipping_reads + 1e-6)
        psi_values.append(psi)
        splicing_labels.append(1 if psi > 0.5 else 0)
        print(f"Total inclusion reads: {inclusion_reads}, Total skipping reads: {skipping_reads}")
        print(f"Calculated PSI: {psi}")

        print(f"Processed Transcript {transcript_id}: PSI = {psi}, Splicing Label = {splicing_labels[-1]}")
        log_memory_usage(f"After processing {transcript_id}")

    return psi_values, splicing_labels

# --- Step 3: Compute PSI and Save Final Dataset ---
print("Processing junction data in chunks...")
psi_values, splicing_labels = process_junction_chunks(junction_file, transcripts_df)

# Add PSI and Splicing_Label to the DataFrame
transcripts_df["PSI"] = psi_values
transcripts_df["Splicing_Label"] = splicing_labels

# Save the final dataset
transcripts_df.to_csv("data/splice_data_with_psi.csv", index=False)
print("✅ Final dataset saved as splice_data_with_psi.csv!")
log_memory_usage("Final dataset saved")
