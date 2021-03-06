import argparse, re, csv, os
import cProfile

def parse_refseq_header(header):
    refseqID, desc = tuple(header.split(' ', 1))

    func, org = tuple(desc.rsplit('[', 1))

    #func = re.sub('[^a-zA-Z0-9-_*. ]', '', func)
    #org = re.sub('[^a-zA-Z0-9-_*. ]', '', org)

    func = func.rstrip()
    org = org.rstrip()

    return refseqID, func, org


def load_refseq_db(refseq_fa):
    # Take refseq fasta headers and convert to dictionary
    db = {}

    seqcounter = 0

    for line in open(refseq_fa):
        if line.startswith(">"):
            seqcounter += 1
            if seqcounter % 1000000 == 0:
                print "processed %d sequences" % seqcounter
            refseqID, func, org = parse_refseq_header(line[1:])

            db[refseqID] = {"org": org, "func": func}

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


def clean_m8_row(row):
    clean_re = re.compile('[^a-zA-Z0-9-_*. ]')

    row['func'] = re.sub(clean_re, '', row['func'])
    row['org'] = re.sub(clean_re, '', row['org'])

    return row


def get_sample_from_name(filename):
    return os.path.basename(filename).split("_S")[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--db', help='database file to parse')
    parser.add_argument('-o', '--output', help='output csv file')
    parser.add_argument('-i', '--input_dir', help='directory of input m8 files')

    args = parser.parse_args()


    refseq_db = load_refseq_db(args.db)

    input_files = [f for f in os.listdir(args.input_dir) if f.endswith("RefSeq_annotated")]
    input_files = [os.path.join(args.input_dir, f) for f in input_files]
    input_samples = [get_sample_from_name(f) for f in input_files]

    test_m8_db = {}
    for f, s in zip(input_files, input_samples):
        print "Processing sample: %s" % s
        parse_m8_file(f, test_m8_db)
    test_m8_db = annotate_m8(test_m8_db, refseq_db)

    fieldnames = ['RefSeqID'] + input_samples + ['org', 'func']

    with open(args.output, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for k in test_m8_db.keys():
            row = clean_m8_row(test_m8_db[k])
            writer.writerow(row)

main()
