"""
Author:     Domitille Jarrige
Date:       2022-05-06
Title:      Frame shift translation
Purpose:    Translate indel mutant genes to predict the impact of the frame-shift
"""

#______________________________________________________
# Modules import
#______________________________________________________

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

import sys, getopt, os, re
import pandas as pd

#______________________________________________________
# Global variables
#______________________________________________________

# NCBI genetic code table from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
CODON_TABLE = "Bacterial"
MUT = r"[ATGC]+"
COOR = r"^[0-9]+"
FILE = r"^[^.]+"

#______________________________________________________
# Functions
#______________________________________________________


def translate_wildtype(nuc_file,  prot_file):
    """Translate the  CDS from a multifasta file and produce 
    a multifasta file of the predicted proteins
    :param nuc_file: multifasta nucleotides or genbank file
    :param prot_file: multifasta amino-acids file"""
    nuc_format = os.path.splitext(nuc_file)[1][1:]
    if nuc_format in ["fa", "fsa", "fna", "fasta"]:
        dico_nuc = SeqIO.to_dict(SeqIO.parse(nuc_file, format="fasta"))
    
    prot_records = []
       
    for rec in dico_nuc:
        ref_seq = SeqRecord(dico_nuc[rec].seq)
        ref_seq.description = f"Automated translation of WT {rec} with NCBI {CODON_TABLE} codon usage table"
        ref_seq.id= dico_nuc[rec].id+ "_translated_prot"
        ref_seq.name = dico_nuc[rec].name + "_translated_prot"
        ref_seq.seq = ref_seq.seq.translate(table=CODON_TABLE)
        prot_records.append(ref_seq)
        
    SeqIO.write(prot_records, prot_file, "fasta")



def translate_mutants(nuc_file, mut_table, nuc_out_file, prot_out_file, report_file):
    """Translate the  CDS from a multifasta file and produce 
    a multifasta file of the predicted proteins
    :param nuc_file: multifasta nucleotides file
    :param prot_file: multifasta amino-acids file
    :param nuc_out_file: multifasta nucleotides file of the mutant genes.
    :param prot_out_file: multifasta amino-acids file of the mutant proteins.
    :param report_file: tsv report file of the mutations effects.
    """
    dico_fasta = SeqIO.to_dict(SeqIO.parse(nuc_file, format="fasta"))
    mut_df = pd.read_csv(mut_table, sep="\t")
    prot_records = []
    nuc_records = []
    list_strains, list_labels, list_coord, list_mut_id, =  [], [], [], []
    list_prot_effects, list_frame_shift, list_wt_prot_len, list_mut_cds_id = [], [], [], []
    
    # Processing each mutant data
    list_mutants = mut_df.groupby("mutant").count().index
    for mutant in list_mutants:
        tmp_df = mut_df[mut_df["mutant"] == mutant]
        gene = ""
        
        for index, row in tmp_df.iterrows():
            gene_tmp = gene
            gene = row["gene_name"]
            try:
                gene_next = tmp_df.loc[index+1, "gene_name"]
            except KeyError:
                gene_next = ""
                
            if gene_tmp == gene:
                # Keeping previous mutated sequence as is
                # Calculating the coordinates of the mutation taking into account previous mutation
                modified_len = len(dico_fasta[gene].seq) - len(sequence)
                coord = row["relative_coord"] - modified_len
                old_coord = row["relative_coord"]
                print(f"For gene {gene}: modified len {modified_len}. old coord {old_coord} -> new coord: {coord}")
                rel_coord.append(str(coord))

            else:
                # Recovering wildtype gene sequence
                sequence = MutableSeq(dico_fasta[gene].seq)
                # Coordinates of the mutation
                coord = row["relative_coord"]
                print(f"coord: {coord}")
                desc = ""
                mut_identity = []
                mut_cds_id = []
                rel_coord = []
                rel_coord.append(str(coord))
            
            # Recovering mutation information
            mut_type = row["mutation_type"]
            print(mut_type)
            variant = re.findall(MUT, mut_type)[-1]

            # Sequence mutation
            sequence, m_id, cds_id = mutate_sequence(sequence, mut_type, variant, coord, frame=int(row["strand"]))
            mut_identity.append(m_id)
            mut_cds_id.append(cds_id)

            if gene == gene_next:
                desc = desc + f"{mut_type} at {coord} and "
                continue
            else:
                # Mutant record generation for previous gene
                # Nucleotides sequence
                mut_seq = SeqRecord(Seq(sequence))
                mut_seq.description = f"Sequence of {mutant}'s {gene}. Mutation(s): " + desc + f"{mut_type} at {coord}"
                mut_seq.id = dico_fasta[gene].id
                mut_seq.name = dico_fasta[gene].name
                nuc_records.append(mut_seq)
            
                # Amino-acids sequence
                mut_prot = SeqRecord(Seq(sequence))
                mut_prot.description = f"Automated translation of {mutant}'s {gene} with NCBI {CODON_TABLE} " \
                                       f"codon usage table. Mutation(s): " + desc + f"{mut_type} at {coord}"
                mut_prot.id = dico_fasta[gene].id + "_translated_prot"
                mut_prot.name = dico_fasta[gene].name + "_translated_prot"
                mut_prot.seq = mut_prot.seq.translate(table=CODON_TABLE)
                prot_records.append(mut_prot)

                # Compute frame_shift:
                try:
                    prefix = os.path.commonprefix([str(mut_prot.seq),
                                                   str(dico_fasta[gene].seq.translate(table=CODON_TABLE))])
                    new_codon = str(mut_prot.seq).removeprefix(prefix)[0]
                    wt_codon = str(dico_fasta[gene].seq.translate(table=CODON_TABLE)).removeprefix(prefix)[0]
                    if wt_codon != new_codon:
                        prot_effect = f"{(len(prefix)) + 1}{wt_codon} > {new_codon}"
                    else:
                        prot_effect = "synonymous"
                except IndexError:
                    prot_effect = "synonymous"

                new_stop = mut_prot.seq.find("*") + 1
                wt_stop = dico_fasta[gene].seq.translate(table=CODON_TABLE).find("*") + 1
                if new_stop == wt_stop:
                    frame_shift = "none"
                elif new_stop == 0:
                    frame_shift = f"no STOP anymore"
                else:
                    frame_shift = f"new STOP codon at {new_stop}"

                list_strains.append(mutant)
                list_labels.append(gene)
                list_coord.append(", ".join(rel_coord))
                list_mut_id.append(", ".join(mut_identity))
                list_mut_cds_id.append(", ".join(mut_cds_id))
                list_prot_effects.append(prot_effect)
                list_frame_shift.append(frame_shift)
                list_wt_prot_len.append(len(dico_fasta[gene].seq.translate(table=CODON_TABLE)))


    report_df = pd.DataFrame({"Mutant": list_strains,
                              "Label": list_labels,
                              "mut_rel_coord": list_coord,
                              "genomic_mutation_identity": list_mut_id,
                              "cds_mutation_identity": list_mut_cds_id,
                              "first_protein_effect": list_prot_effects,
                              "frame_shift": list_frame_shift,
                              "wt_protein_codon_length":list_wt_prot_len})

    # Output files generation:
    report_df.to_csv(report_file, index=False, sep="\t")
    SeqIO.write(prot_records, prot_out_file, "fasta")
    SeqIO.write(nuc_records, nuc_out_file, "fasta")


def mutate_sequence(sequence, mut_type, variant, coord, frame):
    """Perform the sequence mutation at the given coordinates.
    :return sequence: a string of the mutated sequence.
    :return mut_identity: a string of the mutation identity.
    """
    if frame < 0:
        variant = str(Seq(variant).complement())

    if mut_type.endswith("INS"):
        s1 = sequence[:coord]
        s2 = sequence[coord:]
        sequence = s1 + variant + s2
        cds_mut_identity = variant + ">INS"
        if frame < 0:
            variant = str(Seq(variant).complement())
        mut_identity = variant + ">INS"

    elif mut_type.endswith("DEL"):
        reference = sequence[coord - 1:coord - 1 + len(variant)]
        s1 = sequence[:coord - 1]
        s2 = sequence[coord - 1 + len(variant):]
        sequence = s1 + s2
        cds_mut_identity = reference + ">DEL"
        if frame < 0:
            reference = str(Seq(reference).complement())
        mut_identity = reference + ">DEL"

    else:
        substitution = variant
        wild_type = sequence[coord - 1]
        s1 = sequence[:coord - 1]
        s2 = sequence[coord:]
        sequence = s1 + substitution + s2
        if frame < 0:
            substitution = str(Seq(substitution).reverse_complement())
            wild_type = str(Seq(wild_type).reverse_complement())
        cds_mut_identity = wild_type + ">" + substitution
        mut_identity = wild_type + ">" + substitution


    return sequence, str(mut_identity), str(cds_mut_identity)

#______________________________________________________
# Main program
#______________________________________________________

options="hi:t:o:"
long_options=["help", "ref_fasta=", "mutants_table="]

try:
    opts, args = getopt.getopt(sys.argv[1:], options, long_options)
except getopt.GetoptError:
    print("Usage: frame_shift_translation.py -i <ref_cds.fasta> -t <mutants_table.tsv>")
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-h", "--help"):
        print("Usage: frame_shift_translation.py -i <ref_cds.fasta> -t <mutants_table.tsv>")
        sys.exit()
    elif opt in ("-i", "--mutants_fasta"):
        in_file = arg
    elif opt in ("-t", "--mutants_table"):
        tab_file = arg

        
radical = re.match(FILE, os.path.basename(in_file))[0]
ref_prot_file = radical + "_ref_prot.fasta"

radical = re.match(FILE, os.path.basename(in_file))[0]
output_prot = radical + "_mutants_prot.fasta"
output_nuc = radical + "_mutants_gene.fasta"
output_report = radical + "_report.tsv"

print(f"\nWild-type sequences fasta file: {in_file}, \noutput predicted wild-type protein file: {ref_prot_file}\n")
print(f"Mutants data table file: {tab_file}, \noutput mutants genes file: {output_nuc} \n"
      f"output predicted mutated proteins file: {output_prot}\noutput mutation report table: {output_report}")

translate_wildtype(nuc_file = in_file, prot_file = ref_prot_file)
translate_mutants(nuc_file = in_file,  mut_table = tab_file,
                  nuc_out_file = output_nuc, prot_out_file = output_prot, report_file=output_report)

print("\n*********\n* Done! *\n*********\n")
