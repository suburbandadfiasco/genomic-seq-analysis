# Emily Taylor
# SAT4650 Final Project
# Genomic Sequence Analysis

from Bio import Entrez, SeqIO
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import tkinter as tk
from tkinter import messagebox
import os
from collections import Counter
import mysql.connector
import zipfile

# Some good test sequence names:

#Homo sapiens tumor protein p53
#Homo sapiens BRCA1



# NCBI Entrez query for single nucleotide fasta sequence
def fetch_genomic_sequence_fasta():

    global search_term
    search_term = genomic_sequence_entry.get()
    
    # making sure email is valid and has been entered
    user_email = email_entry.get()

    if not validate_email(user_email):
        print("User email in incorrect format. Please enter email in 'example@gmail.com' format")
        return

    Entrez.email = user_email
    
    try:
        # Searching NCBI
        print(f"Searching NCBI for '{search_term}'...")
        handle = Entrez.esearch(db='nucleotide', term=search_term, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            sequence_id = record["IdList"][0]
            print(f"Found sequence ID: {sequence_id}")

            print(f"Fetching sequence {sequence_id} in FASTA format...")
            handle = Entrez.efetch(db='nucleotide', id=sequence_id, rettype="fasta", retmode="text")
            fasta_sequence = handle.read()
            handle.close()

            if validate_fasta(fasta_sequence):
                messagebox.showinfo("Success", "FASTA Sequence retrieved successfully!")
                print(f"FASTA Sequence:\n {fasta_sequence}")
                display_analysis_options(fasta_sequence)
            else:
                messagebox.showerror("Error", "Invalid FASTA format")

        else:
            messagebox.showerror("No Results", "No sequences found for the given query")

        return fasta_sequence
    
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# function to ensure user email is in valid format
def validate_email(user_email):
    email_regex = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    if re.match(email_regex, user_email):
        return True
    else:
        return False

# function to ensure queried nucleotide sequence is valid (title line with carrot and describtion, only ACTG characters otherwise)
def validate_fasta(fasta_string):
    try:
        lines = fasta_string.strip().split("\n")
        if not lines[0].startswith(">"):
            print("Error: FASTA header must start with '>'.")
            return False

        # Track if the current sequence is valid
        for line in lines[1:]:
            if not re.match("^[ATCGatcgN\n]+$", line):
                print(f"Error: Invalid characters found in sequence: {line}")
                return False

        return True
    
    except Exception as e:
        print(f"An error occurred during validation: {e}")
        return False

def display_analysis_options(fasta_sequence):
     # Show buttons for analysis options
    spacer_analysis1 = tk.Frame(root, height = 10)
    spacer_analysis1.pack()
    analysis_label = tk.Label(root, text = "Sequence Analysis Options", font=("Helvetica", 12))
    analysis_label.pack()
    spacer_analysis2 = tk.Frame(root, height = 10)
    spacer_analysis2.pack()

    # Export DNA.fasta file functionality
    export_dna = tk.Button(root, text="Export DNA Sequence", command=lambda:export_fasta(fasta_sequence, "DNA"))
    export_dna.pack()

    # Translate DNA to Protein functionality
    translate_button = tk.Button(root, text="Translate to Protein and Export .FASTA file", command=lambda: translate_seq(fasta_sequence))
    translate_button.pack()

    # Transcribe DNA to RNA functionality
    transcribe_button = tk.Button(root, text="Transcribe to RNA and Export .FASTA file", command=lambda: transcribe_seq(fasta_sequence))
    transcribe_button.pack()

    # Calculate GC (Guanine-Cytosine) Percentage functionality
    gc_button = tk.Button(root, text = "Calculate Guanine-Cytosine content of Sequence", command=lambda: gc_content(fasta_sequence))
    gc_button.pack()

    # Calculate codon frequency functionality
    codon_freq_button = tk.Button(root, text = "Calculate Codon Frequency of Sequence", command=lambda: codon_freq(fasta_sequence))
    codon_freq_button.pack()

    # Database Functionality

    database_frame = tk.Frame(root, height = 10)
    database_frame.pack()
    database_label = tk.Label(root, text = "Database Options", font =("Helvatica", 12))
    database_label.pack()

    # Save Sequence to local database functionality
    save_seq_button = tk.Button(root, text = "Save Sequence to Local Database", command=lambda: save_to_database(fasta_sequence))
    save_seq_button.pack()


def save_to_database(fasta_sequence):

    # connect to database
    db_connection =mysql.connector.connect(
    host = "localhost",
    user = "eetaylor",
    passwd = "xxxxxx",
    database = "genome_analysis")

    # create cursor
    db_cursor = db_connection.cursor()

    # prepare table entry information

    entry_records = list(SeqIO.parse(StringIO(fasta_sequence), "fasta"))
        
    for record in entry_records:
        # Get DNA sequence
        dna_sequence = str(record.seq)
        print(dna_sequence)

        # Transcribe to RNA
        rna_sequence = str(record.seq.transcribe())
        print(rna_sequence)

        # Translate to Protein
        protein_sequence = str(record.seq.translate(to_stop=False))
        print(protein_sequence)

        # Calculate GC percentage
        gc_percentage = gc_content(str(record.seq))
        print(gc_percentage)

        # Calculate Codon Frequency
        codon_frequency = codon_freq(str(record.seq))
  

        # Get first codon frequency as a string (optional, you can change as needed)
        codon_frequency_str = ', '.join([f"{codon}: {freq:.2f}%" for codon, freq in codon_frequency.items()])
        print(codon_frequency_str)


    # Insert into database
        query = """
            INSERT INTO sequences (record_id, description, dna_sequence, rna_sequence, protein_sequence, gc_percentage, codon_frequency)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        values = (
            record.id,
            search_term,  # assuming search_term is passed to describe the sequence
            dna_sequence,
            rna_sequence,
            protein_sequence,
            gc_percentage,
            codon_frequency_str,  # Store codon frequencies as a string
        )

        db_cursor.execute(query, values)

    # Commit and close the connection
    db_connection.commit()
    db_cursor.close()
    db_connection.close()

    print("Sequence data has been saved to the database.")

# Viewing sequences present in database function
def view_database_sequences():
    # connect to database
    db_connection =mysql.connector.connect(
    host = "localhost",
    user = "eetaylor",
    passwd = "ReJw3wd@bq",
    database = "genome_analysis")

    # create cursor
    db_cursor = db_connection.cursor()

    db_cursor.execute("SELECT id, record_id, description FROM sequences")

    result = db_cursor.fetchall()

    print("Sequence(s) present in database: ")

    for entry in result:
        print(entry)

    # close database
    db_cursor.close()
    db_connection.close()

# Extract full sequence report

def sequence_report(sequence_id):
    db_connection =mysql.connector.connect(
    host = "localhost",
    user = "eetaylor",
    passwd = "ReJw3wd@bq",
    database = "genome_analysis")

    # create cursor
    db_cursor = db_connection.cursor()

    query = str(f"SELECT * FROM sequences WHERE id = {sequence_id}")

    db_cursor.execute(query)

    result = db_cursor.fetchall()
    print(result[0][1])
    if not result:
        print(f"No sequence found with ID: {sequence_id}")
        return

    # Extracting all data from the query
    id_ = result[0][1]
    description = result[0][2]
    dna_sequence = result[0][3]
    rna_sequence = result[0][4]
    protein_sequence = result[0][5]
    gc_percentage = result[0][6]
    codon_frequency = result[0][7]
    
    # close database
    db_cursor.close()
    db_connection.close()
    

    # export sequences as .fastas and all other information (gc percentage, codon count, record id, and description) in text doc

    # DNA .fasta file
    dna_fasta_filename = str(f"sequence_{id_}_dna.fasta")
    record = SeqRecord(Seq(dna_sequence), id=str(id_), description=description)
    with open(dna_fasta_filename, "w") as fasta_file:
        SeqIO.write(record, fasta_file, "fasta")
    print(f"DNA sequence exported to {dna_fasta_filename}")

    # RNA .fasta file
    rna_fasta_filename = str(f"sequence_{id_}_rna.fasta")
    record = SeqRecord(Seq(rna_sequence), id=str(id_), description=description)
    with open(rna_fasta_filename, "w") as fasta_file:
        SeqIO.write(record, fasta_file, "fasta")
    print(f"RNA sequence exported to {rna_fasta_filename}")

    # Protein .fasta file
    protein_fasta_filename = str(f"sequence_{id_}_protein.fasta")
    record = SeqRecord(Seq(protein_sequence), id=str(id_), description=description)
    with open(protein_fasta_filename, "w") as fasta_file:
        SeqIO.write(record, fasta_file, "fasta")
    print(f"Protein sequence exported to {protein_fasta_filename}")

    report_filename = str(f"sequence_{id_}_report.txt")
    with open(report_filename, "w") as report_file:
        report_file.write(f"Sequence Report for ID: {id_}\n")
        report_file.write(f"Description: {description}\n")
        report_file.write(f"GC Percentage: {gc_percentage}\n")
        report_file.write("Codon Frequencies:\n")
        for codon in codon_frequency:  # Assuming codon frequencies are stored as a string
            report_file.write(f"  {codon}\n")
    print(f"Report saved to {report_filename}")

    # creation of zip file
    zip_filename = str(f"sequence_{id_}_report.zip")
    with zipfile.ZipFile(zip_filename, "w") as zipf:
        zipf.write(dna_fasta_filename, os.path.basename(dna_fasta_filename))
        zipf.write(rna_fasta_filename, os.path.basename(rna_fasta_filename))
        zipf.write(protein_fasta_filename, os.path.basename(protein_fasta_filename))
        zipf.write(report_filename, os.path.basename(report_filename))
    print(f"All files compressed into {zip_filename}")
    
    
# Export fasta for any molecule type
def export_fasta(fasta_sequence, molecule_type):
    output_filename=str(f"{search_term}_{molecule_type}.fasta")
    try:
        # Ensure the input is in proper FASTA format
        if not fasta_sequence.startswith(">"):
            raise ValueError("Invalid FASTA format: The sequence must start with a header line beginning with '>'.")

        # Write the FASTA sequence to a file
        with open(output_filename, "w") as fasta_file:
            fasta_file.write(fasta_sequence)

        print(f"FASTA sequence successfully exported as '{output_filename}' in {os.getcwd()}")
    
    except Exception as e:
        print(f"An error occurred while exporting FASTA: {e}")


# function to translate DNA nucleotide sequence into protein sequence
def translate_seq(nucleotide_fasta):
    try:
        # Parse the nucleotide FASTA sequence
        nucleotide_records = list(SeqIO.parse(StringIO(nucleotide_fasta), "fasta"))
        
        # Translate each nucleotide sequence into protein
        protein_records = []
        for record in nucleotide_records:
            # Translate the sequence and create a new SeqRecord
            protein_seq = record.seq.translate(to_stop=False)  # Stop at the first stop codon
            protein_record = SeqRecord(
                protein_seq,
                id=record.id,
                description=str(f"{search_term} - translated to protein")
            )
            protein_records.append(protein_record)
        
        # Convert protein records back to FASTA format
        output = StringIO()
        SeqIO.write(protein_records, output, "fasta")
        protein_fasta= output.getvalue()
        print("Protein Sequence:\n", protein_fasta)
        messagebox.showinfo("Translation Complete", "Protein sequence has been translated!")

        export_fasta(protein_fasta, "protein")
        
    except Exception as e:
        messagebox.showerror("Error", f"An error occured during translation: {e}")

# function to transcribe DNA nucleotide sequence into RNA nucleotide sequence
def transcribe_seq(nucleotide_fasta):
    try:
        # Parse the nucleotide FASTA sequence
        nucleotide_records = list(SeqIO.parse(StringIO(nucleotide_fasta), "fasta"))
        
        # Translate each nucleotide sequence into rna
        rna_records = []
        for record in nucleotide_records:
            # Translate the sequence and create a new SeqRecord
            rna_seq = record.seq.transcribe()
            rna_record = SeqRecord(
                rna_seq,
                id=record.id,
                description=str(f"{search_term} - transcribed to RNA"))
            rna_records.append(rna_record)
        
        # Convert rna records back to FASTA format
        output = StringIO()
        SeqIO.write(rna_records, output, "fasta")
        rna_fasta = output.getvalue()
        print("RNA Sequence:\n", rna_fasta)
        messagebox.showinfo("Transcription Complete", "RNA Sequence has been Transcribed!")

        return rna_fasta

        export_fasta(rna_fasta, "RNA")
    
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Function Calculates the Guanine-Cytosine Percentage for a Given Sequence
def gc_content(fasta_sequence):
    try:
        for record in SeqIO.parse(StringIO(fasta_sequence), "fasta"):
            gc_count = record.seq.count("G") + record.seq.count("C")
            gc_percentage = (gc_count /len(record.seq)) *100
            print(f"GC content of {record.id}: {gc_percentage: .2f}%")
            messagebox.showinfo("Calculation Complete", "Guanine-Cytosine Content calculated successfully!")
            return gc_percentage

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred during calculation: {e}")
        return None

# Function Calculates Codon Frequencies for a given sequence
def codon_freq(fasta_sequence):
    try:
        rna_input = transcribe_seq(fasta_sequence)

        sequence_lines = rna_input.strip().splitlines()[1:]
        rna_sequence = ''.join(sequence_lines)

        record = SeqRecord(
            Seq(rna_sequence),
            id="sequence",
            description="Transcribed RNA sequence")

        def get_codons(seq):
            return [seq[i:i+3] for i in range(0, len(seq), 3)]

        codons = get_codons(record.seq)
        codon_counts = Counter(codons)

        total_codons = len(codons)
        codon_frequencies = {codon: (count / total_codons) * 100 for codon, count in codon_counts.items()}

        formatted_frequencies = "\n".join([f"{codon}: {freq:.2f}%" for codon, freq in codon_frequencies.items()])

        messagebox.showinfo("Calculation Complete", f"Codon Frequencies: \n{formatted_frequencies}")
        print(formatted_frequencies)

        return codon_frequencies


    except Exception as e:
        messagebox.showerror("Error", f"An error occurred during calculation: {e}")
        return None

#GUI

# Creating new instance of tkinter
root = tk.Tk()
root.title("Genomic Sequence Query and Analysis Tool")
root.geometry("400x600")

# View sequences from database functionality
view_seq_button = tk.Button(root, text = "View all Sequences in Database", command = lambda: view_database_sequences())
view_seq_button.pack()

# Extract full report on a single sequence from database functionality
extract_full_report_label = tk.Label(root, text = "Extract Full Sequence Report from Database(.zip file)\nEnter sequence id: ")
extract_full_report_label.pack()

extract_full_report_entry = tk.Entry(root)
extract_full_report_entry.pack()

extract_full_report_button = tk.Button(root, text = "Extract Sequence Report", command = lambda: sequence_report(extract_full_report_entry.get()))
extract_full_report_button.pack()


# spacer
spacer_ = tk.Frame(root, height = 10)
spacer_.pack()


# Email input
email_label = tk.Label(root, text = "For sequence search, please enter your email: ")
email_label.pack()
email_entry = tk.Entry(root)
email_entry.pack()

# NCBI query
genomic_sequence_label = tk.Label(root, text = "Enter name of sequence you wish to query: ")
genomic_sequence_label.pack()
genomic_sequence_entry = tk.Entry(root, width = 50)
genomic_sequence_entry.pack()


# Submit button
submit_button = tk.Button(root, text="Submit Query to NCBI", command=fetch_genomic_sequence_fasta)
submit_button.pack()

# Run the GUI
root.mainloop()

