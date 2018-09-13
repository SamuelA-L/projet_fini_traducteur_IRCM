import pandas as pd
import numpy as np


# reads .tsv ornament output reference files to return a panda

def file_to_panda_ref(filePath):
    ref_panda = pd.read_csv(filePath, sep="\t", dtype={"Transcript stable ID": str,
                                                       "Exon rank in transcript": int,
                                                       "Chromosome/scaffold name": str,
                                                       "cDNA coding start": float,
                                                       "cDNA coding end": float,
                                                       "Exon region start (bp)": int,
                                                       "Exon region end (bp)": int,
                                                       "5' UTR start": float,
                                                       "5' UTR end": float,
                                                       "3' UTR start": float,
                                                       "3' UTR end": float,
                                                       "CDS start": float,
                                                       "CDS end": float,
                                                       "Strand": float})
    return ref_panda


# writes whole oRNAment panda to an output file with name of the transcript

def panda_to_file(panda, fileName):
    outputFile = lib_dir + fileName + "_processed.tsv"
    panda.to_csv(outputFile, sep="\t", index=False, mode="w", columns=["Transcript stable ID",
                                                                       "Exon rank in transcript",
                                                                       "Chromosome/scaffold name",
                                                                       "cDNA coding start",
                                                                       "cDNA coding end",
                                                                       "Exon region start (bp)",
                                                                       "Exon region end (bp)",
                                                                       "5' UTR start",
                                                                       "5' UTR end",
                                                                       "3' UTR start",
                                                                       "3' UTR end",
                                                                       "CDS start",
                                                                       "CDS end",
                                                                       "Strand",
                                                                       "cDNA exon start",
                                                                       "cDNA exon end"])


# calculate cdna addresses for completely empty transcripts and missing last address

def complete_fruitfly_empty():

    reference_fruitfly = file_to_panda_ref(lib_dir + "sorted_fruitfly.tsv")
    fruitfly_transcripts = file_to_panda_ref(original_dir + "fruitfly_transcripts.tsv")
    addition_start = []
    addition_end = []

    for i in fruitfly_transcripts["Transcript stable ID"]:

        transcript_df = reference_fruitfly.loc[reference_fruitfly["Transcript stable ID"] == i]
        print(i)
        last_row_val = 0
        row_start = 0
        row_end = 0

        for j, row in transcript_df.iterrows():

            if row[1] == 1:
                row_start = 1
                row_end = row_start + (row[6] - row[5])

            else:

                row_start = last_row_val+1
                row_end = row_start + row[6] - row[5]

            last_row_val = row_end
            addition_start.append(row_start)
            addition_end.append(row_end)

    reference_fruitfly["cDNA exon start"] = addition_start
    reference_fruitfly["cDNA exon end"] = addition_end

    reference_fruitfly = reference_fruitfly.replace(0, np.nan,)
    panda_to_file(reference_fruitfly, "exon_reference_fruitfly")

    print(reference_fruitfly)


def complete_human_empty():

    reference_human = file_to_panda_ref(lib_dir + "sorted_human.tsv")
    human_transcripts = file_to_panda_ref(original_dir + "human_transcripts.tsv")
    addition_start = []
    addition_end = []
    print(reference_human.shape)

    for i in human_transcripts["Transcript stable ID"]:

        transcript_df = reference_human.loc[reference_human["Transcript stable ID"] == i]
        print(i)
        last_row_val = 0
        row_start = 0
        row_end = 0

        for j, row in transcript_df.iterrows():

            if row[1] == 1:
                row_start = 1
                row_end = row_start + (row[6] - row[5])

            else:

                row_start = last_row_val + 1
                row_end = row_start + row[6] - row[5]

            last_row_val = row_end
            addition_start.append(row_start)
            addition_end.append(row_end)

    reference_human["cDNA exon start"] = addition_start
    reference_human["cDNA exon end"] = addition_end

    reference_human = reference_human.replace(0, np.nan, )
    panda_to_file(reference_human, "exon_reference_human")

    print(reference_human)


lib_dir = "../lib/"
original_dir = "../original files/"


# fill info files
complete_fruitfly_empty()
complete_human_empty()
