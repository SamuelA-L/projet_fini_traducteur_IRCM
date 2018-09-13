import requests
import sys
import re


def rest_api(transcript_id,cdna_range_start):   # base model to pull traduction from ensembl

    cdna_range_end = cdna_range_start

    server = "http://rest.ensembl.org"
    ext = "/map/cdna/" + transcript_id + "/" + cdna_range_start + ".." + cdna_range_end  # + "?/species=" + species

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    ensembl_output = (repr(decoded))
    print(ensembl_output)

    output_array = []
    output_array = re.split(":|,|}",ensembl_output)

    print("address : ", cdna_range_start, "\n on genome :  ", int(output_array[16]), "\n")
    print(80*"*", "\n")


addreses = [1, 223, 224, 353, 354, 1103, 1104, 2589, 2590, 3264]


# +1 au cDNA address
rest_api("ENST00000155093", "289")
