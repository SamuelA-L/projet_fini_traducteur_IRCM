import re
import requests
import sys
import pandas as pd
import subprocess
import os
import math
import numpy as np
import datetime
import ujson
from pandas.io.json import json_normalize
import gzip


# reference arrays for all motifs classified by length

length_4_motifs = ["M344", "M320", "M316", "M292", "M286"]
length_5_motifs = ["M246", "M247", "M333", "M307", "M297", "M275"]
length_6_motifs = ["M354", "M250", "M254", "M261", "M271", "M256", "M298"]
length_7_motifs = ["M176", "M160", "M161", "M162", "M163", "M164", "M165", "M166", "M167", "M168", "M169", "M170",
                   "M171", "M172", "M173", "M174", "M184", "M182", "M181", "M179", "M159", "M175", "M178", "M177",
                   "M158", "M148", "M156", "M128", "M130", "M131", "M132", "M133", "M135", "M137", "M138", "M139",
                   "M140", "M141", "M142", "M143", "M144", "M145", "M146", "M187", "M149", "M150", "M151", "M153",
                   "M154", "M155", "M157", "M189", "M202", "M193", "M227", "M228", "M229", "M230", "M231", "M232",
                   "M234", "M235", "M236", "M238", "M239", "M241", "M242", "M243", "M244", "M262", "M273", "M290",
                   "M291", "M318", "M319", "M348", "M349", "M226", "M225", "M224", "M223", "M194", "M195", "M196",
                   "M198", "M199", "M201", "M127", "M203", "M204", "M206", "M208", "M190", "M209", "M211", "M213",
                   "M214", "M215", "M216", "M217", "M218", "M219", "M220", "M221", "M222", "M210", "M126", "M001",
                   "M124", "M039", "M040", "M041", "M042", "M043", "M044", "M045", "M046", "M047", "M048", "M049",
                   "M050", "M051", "M052", "M053", "M054", "M055", "M056", "M125", "M058", "M059", "M060", "M061",
                   "M062", "M063", "M064", "M065", "M066", "M067", "M037", "M068", "M036", "M034", "M002", "M003",
                   "M004", "M005", "M007", "M008", "M009", "M011", "M012", "M013", "M014", "M015", "M016", "M017",
                   "M018", "M019", "M020", "M021", "M022", "M023", "M024", "M025", "M026", "M027", "M028", "M029",
                   "M030", "M031", "M032", "M035", "M069", "M057", "M147", "M087", "M088", "M089", "M090", "M091",
                   "M092", "M093", "M094", "M095", "M114", "M096", "M097", "M113", "M098", "M099", "M112", "M100",
                   "M101", "M111", "M102", "M103", "M104", "M105", "M109", "M115", "M086", "M083", "M123", "M072",
                   "M121", "M073", "M120", "M084", "M075", "M076", "M108", "M077", "M082", "M078", "M117", "M081",
                   "M079", "M116"]
length_8_motifs = ["M118", "M129", "M110", "M122", "M006", "M272", "M119", "M317", "M010", "M136", "M134", "M070",
                   "M107", "M191", "M071", "M197", "M188", "M200", "M186", "M185", "M183", "M205", "M207", "M180",
                   "M212", "M074", "M085", "M106", "M192", "M080", "M240", "M237", "M033", "M152", "M233", "M038"]
length_9_motifs = ["M353", "M352", "M350", "M274", "M330", "M329"]
length_10_motifs = ["M279", "M331", "M332", "M346", "M351"]
length_11_motifs = ["M269", "M345", "M347", "M296", "M343"]
length_12_motifs = ["M299", "M300"]
length_13_motifs = ["M323", "M260"]
length_14_motifs = ["M328"]


# Translates cDNA addresses to bp addresses with Ensembl's rest API

def rest_api(transcript_start, transcript_name):

    cdna_range_start = int(transcript_start)
    cdna_range_end = int(cdna_range_start)

    # pasted code from REST api web page
    server = "http://rest.ensembl.org"
    ext = "/map/cdna/" + transcript_name + "/" + str(cdna_range_start) + ".." + str(cdna_range_end)

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    ensembl_output = (repr(decoded))

    # split info array to find translation address
    output_array = re.split(":|,|}", ensembl_output)
    genome_address_start = int(output_array[16])
    print("start : ", genome_address_start)

    return genome_address_start


# reads .tsv oRNAment output files and returns a panda dataframe

def file_to_panda(filePath):
    columns = ["address_cDNA", "motif", "score"]  # define column names
    ornament_example = pd.read_csv(filePath, sep="\t", names=columns, header=None)
    return ornament_example


def json_to_panda(filepath):

    with gzip.open(filepath) as ujson_data_file:
        data_D = ujson.load(ujson_data_file)
    pd_dataframe = json_normalize(data_D, ['stats', 'RBP'], ['id', ['stats', 'position']], errors='ignore')
    return pd_dataframe


# reads .tsv exon addresses reference files with specified data type to return a panda dataframe

def ref_file_to_panda(filePath):
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
                                                       "Strand": float})  # desired data type read
    return ref_panda


# reads compendium .csv file to return a panda dataframe

def compendium_to_panda(filePath):

    compendium_panda = pd.read_csv(filePath, sep=",")
    return compendium_panda


# writes specified columns in oRNAment dataframe to an output file with the name of the transcript

def panda_to_file(panda, original_fileName):
    outputFile = output_dir + original_fileName + "_translated.tsv"  # name for output file
    panda.to_csv(outputFile, sep="\t", header=None, index=False,
                 columns=["Chromosome",
                          "stats.position",
                          "translation start",
                          "translation end",
                          "score",
                          "motif",
                          "id",
                          "duplicate line"])       # desired output dataframe columns


# executes a terminal command

def terminal_command(command):

    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()


# returns small dataframe with only the lines containing the info of one transcript

def extract_transcript_info(transcript_name):

    global output_dir

    # look for transcript in fruitfly reference file
    if transcript_name in reference_fruitfly["Transcript stable ID"].values:   # find transcript in fruitfly

        print(transcript_name, "found in fruitfly")

        # assign transcript dataframe to return variable
        transcript_info = reference_fruitfly.loc[reference_fruitfly["Transcript stable ID"] == transcript_name]
        output_dir = "../output/fruitfly/"

    # look for transcript in human reference file
    elif transcript_name in reference_human["Transcript stable ID"].values:    # find transcript in human

        print(transcript_name, "found in human")

        # assign transcript dataframe to return variable
        transcript_info = reference_human.loc[reference_human["Transcript stable ID"] == transcript_name]
        output_dir = "../output/human/"

    # if transcript isn't found in either
    else:

        print("transcript not found : ", transcript_name)

    return transcript_info


# returns motif length for a given motif with reference arrays(only needed if motif length vary)

def extract_motif_length(motif_id):

    if motif_id in length_7_motifs:
        return 7
    elif motif_id in length_8_motifs:
        return 8
    elif motif_id in length_6_motifs:
        return 6
    elif motif_id in length_5_motifs:
        return 5
    elif motif_id in length_4_motifs:
        return 4
    elif motif_id in length_9_motifs:
        return 9
    elif motif_id in length_10_motifs:
        return 10
    elif motif_id in length_11_motifs:
        return 11
    elif motif_id in length_12_motifs:
        return 12
    elif motif_id in length_13_motifs:
        return 13
    elif motif_id in length_14_motifs:
        return 14


# returns cDNA bp start address, end address + info on the transcripts start and end position in an array for 1 address

def extract_transcript_addresses(transcript_info, motif_start, motif_length, last_translated_exon):

    motif_start += 1  # conversion from base 0 to Ensembl's base 1
    motif_end = motif_start + motif_length - 1  # motif end address in Ensembl base
    info = np.empty([2], dtype="<U4")
    start_translation = np.empty([0], dtype=int)
    end_translation = np.empty([0], dtype=int)
    two_lines = np.empty([0], dtype=bool)
    last_exon_in_transcript = transcript_info.tail(1)["Exon rank in transcript"].item()

    # resize to 3 rows for the 3 necessary exons(ex: exon and i+1 for pos strand or exoni and i-1 for neg strand)
    resized_transcript_info = transcript_info.loc[transcript_info["Exon rank in transcript"] <= last_translated_exon+2]

    for i, row in resized_transcript_info.iterrows():

        # variables to identify information in each row
        transcript = row[0]
        exon_rank = row[1]
        chromosome = "chr" + str(row[2])
        cDNA_coding_start = row[3]
        cDNA_coding_end = row[4]
        genome_start = row[5]
        genome_end = row[6]
        five_utr_start = row[7]
        five_utr_end = row[8]
        three_utr_start = row[9]
        three_utr_end = row[10]
        strand = row[13]
        cDNA_start = math.floor(float(row[14]))  # remove .0 float part
        cDNA_end = math.floor(float(row[15]))

        # if start address and end address are in first exon
        if (motif_start in range(cDNA_start, cDNA_end+1)) and (motif_end in range(cDNA_start, cDNA_end+1)):

            if strand == 1:     # translation for positive strands

                start_translation = np.append(start_translation, rest_api(motif_start, transcript))
                end_translation = np.append(end_translation, rest_api(motif_end, transcript))

            elif strand == -1:  # translation for negative strands

                end_translation = np.append(end_translation, rest_api(motif_start, transcript))
                start_translation = np.append(start_translation, rest_api(motif_end, transcript))

            last_translated_exon = exon_rank
            break

        # if start address and end address are on different exons
        else:

            if motif_start in range(cDNA_start, cDNA_end + 1):

                if strand == 1:  # translation for positive strands

                    start_translation = np.append(start_translation, rest_api(motif_start, transcript))
                    end_translation = np.append(end_translation, genome_end)

                elif strand == -1:  # translation for negative strands

                    end_translation = np.append(end_translation, rest_api(motif_start, transcript))
                    start_translation = np.append(start_translation, genome_start)

                last_translated_exon = exon_rank

            # if address is in exon (math function to make int for range)
            if motif_end in range(cDNA_start, cDNA_end + 1):

                if strand == 1:  # translation for positive strands

                    end_translation = np.append(end_translation, rest_api(motif_end, transcript))
                    start_translation = np.append(start_translation, genome_start)

                elif strand == -1:  # translation for negative strands

                    start_translation = np.append(start_translation, rest_api(motif_end, transcript))
                    end_translation = np.append(end_translation, genome_end)

                last_translated_exon = exon_rank

        # for motifs ending out of the last exon
        if exon_rank == last_exon_in_transcript:

            # if array is empty
            if start_translation.size == 0:
                start_translation = np.append(start_translation, transcript_info.tail(1)["Exon region start (bp)"].item())

            if not end_translation.size == 0:
                end_translation = np.append(end_translation, transcript_info.tail(1)["Exon region end (bp)"].item())

    for j in range(len(start_translation)):

        # position info for start
        if cDNA_coding_start <= motif_start <= cDNA_coding_end:
            info[j] += "C"
        if five_utr_start <= start_translation[j] <= five_utr_end:
            info[j] += "5"
        if three_utr_start <= start_translation[j] <= three_utr_end:
            info[j] += "3"

        # position info for end
        if cDNA_coding_start <= motif_end <= cDNA_coding_end:
            info[j] += "C"
        if five_utr_start <= end_translation[j] <= five_utr_end:
            info[j] += "5"
        if three_utr_start <= end_translation[j] <= three_utr_end:
            info[j] += "3"

    # indicate if the translation is on same exon or over 2 different exons
    if len(start_translation) == 1:
        two_lines = np.append(two_lines, False)
    elif len(start_translation) == 2:
        two_lines = np.append(two_lines, True)

    return [chromosome, start_translation, end_translation, info, last_translated_exon, two_lines]


# execute the translation for all the addresses in a transcript oRNAment file and return 4 arrays to create datasets

def translation(dataset, transcript):

    array_length = dataset.shape[0]  # length of oRNAment file
    chromosome_array = np.empty([array_length], dtype="<U32")
    translation_start = np.empty([array_length], dtype=int)
    translation_end = np.empty([array_length], dtype=int)
    info_array = np.empty([array_length], dtype="<U4")  # data type -> string of length 2
    duplicate = np.empty([array_length], dtype=bool)
    last_exon = 1
    dataset = dataset.reset_index(drop=True)
    i = 0

    # extract info dataframe for the transcript being translated
    transcript_info = extract_transcript_info(transcript)

    for k, row in dataset.iterrows():

        same_start_address = row[4]
        cDNA_address_start = row[3]
        motif_length = 7  # choose motif length for translation of ending address

        # if last address is the same as current one
        if same_start_address:

            for j in range(len(last_start)):  # for each address in array (2 lines added if exon splice)

                translation_start[i] = last_start[j]
                translation_end[i] = last_end[j]
                info_array[i] = last_info[j]
                chromosome_array[i] = chromosome
                i += 1

            # to know which lines to duplicate
            duplicate[k] = two_lines

        else:

            # use function to extract addresses and info on a starting cDNA address for a transcript
            extract = extract_transcript_addresses(transcript_info, cDNA_address_start, motif_length, last_exon)

            chromosome = extract[0]
            start = extract[1]
            end = extract[2]
            info = extract[3]
            last_exon = extract[4]
            two_lines = extract[5]

            for j in range(len(start)):  # for each address in array (2 lines added if exon splice)

                translation_start[i] = start[j]
                translation_end[i] = end[j]
                info_array[i] = info[j]
                chromosome_array[i] = chromosome
                last_start = start
                last_end = end
                last_info = info
                i += 1

            # to know which lines to duplicate
            duplicate[k] = two_lines

        if two_lines:
            # increase array sizes to match final length with exon splice
            chromosome_array.resize(len(chromosome_array) + 1)
            translation_start.resize(len(translation_start) + 1)
            translation_end.resize(len(translation_end) + 1)
            info_array.resize(len(info_array) + 1)

    return [chromosome_array, translation_start, translation_end, info_array, duplicate]


# main loop

def main():

    translated_files = (os.listdir(output_human) + os.listdir(output_fruitfly))  # combine output directories

    for file in os.listdir(src_dir):

        # exclude hidden files of the loop
        # and
        # verify if transcript is already translated (presence of file in output dir)
        if file.endswith(".gz") and ((file[:-9] + "_translated.tsv") not in translated_files):

            file_path = str(src_dir + file)  # generate file path
            ornament_panda = json_to_panda(file_path)

            # reorder panda to regroup all identical starting addresses
            sorted_ornament = ornament_panda.sort_values(by=["id", "stats.position", "motif", "score"], axis=0, ascending=True)

            # reorder line index column for translation loop index
            sorted_ornament = sorted_ornament.reset_index(drop=True)

            # new boolean column indicating if last cDNA starting address is the same as the current one
            sorted_ornament["mmAddress"] = sorted_ornament["stats.position"].diff().eq(0)

            # create arrays to fill with translation data
            file_transcripts = sorted_ornament["id"].unique()
            translation_chr = np.empty([0], dtype="<U32")
            translation_start = np.empty([0], dtype=int)
            translation_end = np.empty([0], dtype=int)
            translation_info = np.empty([0], dtype="<U4")
            two_lines = np.empty([0], dtype=bool)

            # try executing translation
            for transcript in file_transcripts:

                try:
                    transcript_data = sorted_ornament.loc[sorted_ornament["id"] == transcript]

                    # if transcript is human
                    if transcript[0:4] == "ENST":
                        transcript = transcript[0:15]  # extract only transcript id(not ending .##)

                    addition = translation(transcript_data, transcript)

                    # add avery transcript translation info to the arrays
                    translation_chr = np.append(translation_chr, addition[0])
                    translation_start = np.append(translation_start, addition[1])
                    translation_end = np.append(translation_end, addition[2])
                    translation_info = np.append(translation_info, addition[3])
                    two_lines = np.append(two_lines, addition[4])

                # if transcript is not found in reference files
                except UnboundLocalError:

                    # write file name in list and move file to "not found" directory
                    not_found_file.write(file+" : " + transcript + "\n")

                    # create fillers for not found values
                    nan_filler = np.full(transcript_data.shape[0], np.nan)  # won't be removed when empty lines are
                    filler = np.empty(transcript_data.shape[0])  # empty data

                    translation_chr = np.append(translation_chr, nan_filler)  # column is used to find empty cells
                    translation_start = np.append(translation_start, filler)
                    translation_end = np.append(translation_end, filler)
                    translation_info = np.append(translation_info, filler)
                    two_lines = np.append(False, two_lines)

            print("trad fini")

            # add series to determine lines that have to be duplicated for translation over 2 different exons
            lines_add = np.count_nonzero(two_lines)  # number of lines to add to dataframe
            sorted_ornament["duplicate line"] = two_lines

            print("lines_add :", lines_add)

            # extract lines to duplicate form sorted ornament
            missing_data = sorted_ornament.loc[sorted_ornament["duplicate line"] == True]

            print("ajout ligne fini")

            # add lines that need to be duplicated to the dataframe
            sorted_ornament = pd.concat([missing_data, sorted_ornament], ignore_index=True)

            print("append fini")

            # reorder lines to place the lines that were just added to the dataframe
            sorted_ornament = sorted_ornament.sort_values(by=["id", "stats.position", "motif", "score"], axis=0, ascending=True)
            print("sort fini")

            # create new dataframe columns with numpy arrays
            sorted_ornament["Chromosome"] = translation_chr
            sorted_ornament["translation start"] = translation_start
            sorted_ornament["translation end"] = translation_end
            sorted_ornament["info"] = translation_info

            print("pret a envoyer en tsv")

            panda_to_file(sorted_ornament, file[:-9])

            print(sorted_ornament.head())
            print("\nlength :  " + str(sorted_ornament.shape[0]) + "\n")
            print(sorted_ornament.info())



# define used directory for file reading and generating
output_dir = "../output/"
lib_dir = "../lib/"
src_dir = "../src/"
original_dir = "../original files"
rbp_classified_dir_human = "../rbp_classified/human/"
rbp_classified_dir_fruitfly = "../rbp_classified/fruitfly/"
output_human = output_dir + "human/"
output_fruitfly = output_dir + "fruitfly/"

# define reference files

# used to find cDNA addresses and recover bp genome addresses
reference_human = ref_file_to_panda(lib_dir + "exon_reference_human_processed.zip")
reference_fruitfly = ref_file_to_panda(lib_dir + "exon_reference_fruitfly_processed.zip")

compendium = compendium_to_panda("../lib/sorted_compendium.csv")

# define file and info for the file keeping track of transcript that are not found in reference files
not_found_file = open("../not found/not found.txt", "a")
time_info = str(datetime.datetime.now())
not_found_file.write(time_info + "\n")

# loop to keep process running until all files are translated or classified as "not found"(not necessary for cluster)
lenSrc = len(os.listdir(src_dir))
lenOutput = len(os.listdir("../output/fruitfly") + os.listdir("../output/human"))

# until there is the same amount of output files as source files
while True:

    main()




