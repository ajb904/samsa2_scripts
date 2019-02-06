import argparse, re, csv, os
from Bio import SeqIO

def get_org(fasta):
    org = fasta.description.rsplit("[", 1)[1]

    return org


def get_func(fasta):
    # Function is everything between the refseq ID (first word of definition)
    # and the organism (defined by square brackets)
    func = fasta.description.rsplit("[", 1)[0]
    func = func.split(" ", 1)[1]
    func = func.rstrip()

    return func


def clean(annotation, regexp):
    return re.sub(regexp, '', annotation)


def load_refseq_db(refseq_fa):
    # Take refseq SeqIO generator and convert to dictionary
    db = {}

    clean_re = re.compile('[^a-zA-Z0-9-_*. ]')

    seqcounter = 0

    for seq in refseq_fa:
        seqcounter += 1
        if seqcounter % 1000000 == 0:
            print "processed %d sequences" % seqcounter
        org = clean(get_org(seq), clean_re)
        func = clean(get_func(seq), clean_re)
        db[seq.id] = {"org": org, "func": func}

    return db


def parse_m8_file(m8_file, m8_db):
    # Read an m8 file and return as a dictionary with refseq IDs as keys and
    # counts as values.
    sample = get_sample_from_name(m8_file)

    for line in open(m8_file):
        refseqID = line.split("\t")[1]
        if m8_db.has_key(refseqID):
            try:
                m8_db[refseqID][sample] += 1
            except KeyError:
                m8_db[refseqID][sample] = 1
        else:
            m8_db[refseqID] = {"RefSeqID": refseqID, sample: 1}


def annotate_m8(m8_dict, refseq_dict):
    for k in m8_dict.keys():
        try:
            refseq_entry = refseq_dict[k]
            m8_dict[k]["org"] = refseq_entry["org"]
            m8_dict[k]["func"] = refseq_entry["func"]
        except KeyError:
            m8_dict[k]["org"] = "NA"
            m8_dict[k]["func"] = "NA"

    return m8_dict


def get_sample_from_name(filename):
    return os.path.basename(filename).split("_S")[0]


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--db', help='database file to parse')
parser.add_argument('-o', '--output', help='output csv file')
parser.add_argument('-i', '--input_dir', help='directory of input m8 files')

args = parser.parse_args()

refseq_db_fasta = SeqIO.parse(args.db, 'fasta')

refseq_db = load_refseq_db(refseq_db_fasta)

input_files = [f for f in os.listdir(args.input_dir) if f.endswith("RefSeq_annotated")]
input_files = [os.path.join(args.input_dir, f) for f in input_files]
input_samples = [get_sample_from_name(f) for f in input_files]

# test_m8 = 'test_m8/1344_S1_L001.merged.RefSeq_annotated'
test_m8_db = {}
for f, s in zip(input_files, input_samples):
    print "Processing sample: %s" % s
    parse_m8_file(f, test_m8_db)
test_m8_db = annotate_m8(test_m8_db, refseq_db)

n=0
fieldnames = ['RefSeqID'] + input_samples + ['org', 'func']
with open(args.output, 'w') as csvfile:
    writer =csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for k in test_m8_db.keys():
        if n < 20:
            print k, test_m8_db[k]
            n+=1
        writer.writerow(test_m8_db[k])


# outfile = open(args.output, 'w')

# for seq in db_fasta:
#     if org[0].isdigit():
#         print seq.description
#     outfile.write('%s,%s,%s\n' % (seq.id,org,func))
# outfile.close()
