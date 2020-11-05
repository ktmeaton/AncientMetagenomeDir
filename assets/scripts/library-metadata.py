# -----------------------------------------------------------------------------#
# Environment
# -----------------------------------------------------------------------------#
# conda:
#   conda create -n library-metadata python=3.8 requests=2.24.0
#   conda activate library-metadata
# pip:
#   pip install requests==2.24.0

# -----------------------------------------------------------------------------#
# Packages
# -----------------------------------------------------------------------------#
import argparse  # Command-line argument parsing
import os  # file path searching and creation
import urllib # Catch HTTP Errors
import requests # ENA API

# ------------------------------------------------------------------------------#
# Argument Parsing                                                             #
# ------------------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description=("Add ENA metadata to AncientMetagenomeDir TSV files."),
    add_help=True,
)

parser.add_argument(
    "--input",
    help="Input TSV file of AncientMetagenomeDir sample metadata.",
    action="store",
    dest="inputPath",
    required=True,
)

parser.add_argument(
    "--output",
    help="Output TSV file of AncientMetagenomeDir library metadata.",
    action="store",
    dest="outputPath",
    required=True,
)

args = vars(parser.parse_args())
input_path = args["inputPath"]
output_path = args["outputPath"]

# -----------------------------------------------------------------------------#
# Setup
# -----------------------------------------------------------------------------#

METADATA_COLUMNS = [
                  "Sample_Name",
                  "Library_ID",
                  "Lane",
                  "Colour_Chemistry",
                  "SeqType",
                  "Organism",
                  "Strandedness",
                  "UDG_Treatment",
                  "R1",
                  "R2",
                  "BAM",
                  "Instrument_Platform",
                  "Library_Source",
                  "Library_Strategy",
                  "Read_Count",
        ]
METADATA_HEADER = "\t".join(METADATA_COLUMNS)
metadata_dict = {}


# ENA file will be retrieved with the ENA api
ENA_SAMPLE_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields={}"
# Currently this is just for mental conversion to the nf-core/eager format
ENA_FIELDS_MAP = {
        "run_accession"  : "Library_ID",
        "library_layout" : "SeqType",
        "fastq_ftp" : "R1,R2",
        "instrument_model" : "Instrument_Model",
        "instrument_platform" : "Instrument_Platform",
        "library_source" : "Library_Source",
        "library_strategy" : "Library_Strategy",
        "read_count" : "Read_Count",
        }
ENA_FIELDS_LIST = list(ENA_FIELDS_MAP.keys())
ENA_FIELDS_CSV = ",".join(ENA_FIELDS_LIST)

# Validate I/O
try:
  input_file = open(input_path, "r")
except FileNotFoundError:
  print("Invalid input file: {}".format(input_path))
  exit(1)

output_file_exists = False
if os.path.isfile("output_path"):
    output_file_exists = True
try:
  output_file = open(output_path, "w")
  output_file.close()
except FileNotFoundError:
  print("Invalid output file: {}".format(output_path))
  exit(1)

# -----------------------------------------------------------------------------#
# Processing - Output file does not already exists.
# -----------------------------------------------------------------------------#
if not output_file_exists:
  # Write the new header to the output file
  with open(output_path, "w") as outfile:
      outfile.write(METADATA_HEADER + "\n")
  # File is tsv - perhaps verify first?
  file_data = [line.strip().split("\t") for line in input_file]
  # Header is the first line
  file_header = file_data[0]
  # figure out which columns are accesion, name
  name_index = file_header.index("sample_name")
  acc_index = file_header.index("archive_accession")
  # Depeneding on the tsv, there may not be an organism
  if "singlegenome_species" in file_header:
    organism_index = file_header.index("singlegenome_species")
  else:
    organism_index = None

  # Exclude header from data
  file_data = file_data[1:]

  # Iterate through the rows/samples in the tsv file
  for sample_metadata in file_data:
    # Should sample_name instead be the sra sample accession?
    sample_name = sample_metadata[name_index]
    # Remember, organism is not always present
    if organism_index:
      sample_organism = sample_metadata[organism_index]
    else:
      sample_organism = "N/A"

    # Add sample name to dictionary
    metadata_dict[sample_name] = {}

    sample_acc_csv = sample_metadata[acc_index]
    # Remove whitespaces from acc (possibly after commas)
    sample_acc_csv = sample_acc_csv.replace(" ","")

    # Split up the csv values of acc
    sample_acc_split = sample_acc_csv.split(",")

    # Iterate through the sample accessions
    for sample_acc in sample_acc_split:

      # Here are the default values, in case we can't find them
      run_accession = "N/A"
      run_lane = "N/A"
      color_chemistry = "N/A"
      seq_type = "N/A"
      strandedness = "N/A"
      udg_treatment = "N/A"
      R1 = "N/A"
      R2 = "N/A"
      bam = "N/A"
      instrument_platform = "N/A"
      library_source = "N/A"
      library_strategy = "N/A"
      read_count = "N/A"

      # Get the ENA metadata
      result = requests.get(url = ENA_SAMPLE_URL.format(sample_acc, ENA_FIELDS_CSV))

      # Status code 200 is success
      if result.status_code == 200:
        # The result is multiple lines, exclude header line
        data = result.text.strip().split("\n")[1:]
        split_data = [col.split("\t") for col in data]
        for sample_run in split_data:
            # Parse desired columns from ENA
            run_accession = sample_run[ENA_FIELDS_LIST.index("run_accession")]
            run_lane = "N/A"
            run_model = sample_run[ENA_FIELDS_LIST.index("instrument_model")]
            color_chemistry = (4 if "HiSeq" or "MiSeq" in run_model else 
                               2 if "NextSeq" or "NovaSeq" in run_model else
                               "NA")
            # Parse library layout
            library_layout = sample_run[ENA_FIELDS_LIST.index("library_layout")]
            seq_type = ("PE" if "PAIRED" in library_layout else 
                        "SE" if "SINGLE" in library_layout else 
                        "NA")
            # Fastq ftps will be a semicolon (;) separated string
            fastq_ftp_split = sample_run[ENA_FIELDS_LIST.index("fastq_ftp")].split(";")
            # R1 will always be the first element, R2 is optinal
            R1 = fastq_ftp_split[0]
            # If there's a second item, it's paired end
            if len(fastq_ftp_split) == 2:
                R2 = fastq_ftp_split[1]
            # Parse the library info
            instrument_platform = sample_run[ENA_FIELDS_LIST.index("instrument_platform")]
            library_source = sample_run[ENA_FIELDS_LIST.index("library_source")]
            library_strategy = sample_run[ENA_FIELDS_LIST.index("library_strategy")]
            read_count = sample_run[ENA_FIELDS_LIST.index("read_count")]

      # The final metadata for the library!
      with open(output_path, "a") as outfile:
         outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                                  sample_name,
                                  run_accession,
                                  run_lane,
                                  color_chemistry,
                                  seq_type,
                                  sample_organism,
                                  strandedness,
                                  udg_treatment,
                                  R1,
                                  R2,
                                  bam,
                                  instrument_platform,
                                  library_source,
                                  library_strategy,
                                  read_count) + "\n")
