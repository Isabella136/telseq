from Bio import SeqIO
import pysam
from common import *

def overlapsWith(mge_query_pos, arg_query_pos):
    if mge_query_pos[0] < arg_query_pos[0]:
        if mge_query_pos[1] > arg_query_pos[0]:
            return True
    elif mge_query_pos[1] > arg_query_pos[1]:
        if mge_query_pos[0] < arg_query_pos[1]:
            return True
    else:
        return True
    return False

def overlap_finder(config):
    overlap_threshold = float(config['MISC']['OVERLAP_THRESHOLD'])

    mge_genes_length = dict()
    mge_fasta = config['DATABASE']['MGES']
    for rec in SeqIO.parse(mge_fasta, 'fasta'):
        mge_genes_length[rec.name] = len(rec.seq)

    arg_genes_length = dict()
    arg_fasta = config['DATABASE']['MEGARES']
    for rec in SeqIO.parse(arg_fasta, 'fasta'):
        arg_genes_length[rec.name] = len(rec.seq)

    read_mge_alignments = dict()
    mgeSAM = pysam.AlignmentFile(config['INPUT']['MGES_SAM_FILE'])
    for read in mgeSAM.fetch():
        if read.is_unmapped: continue
        if config['MISC']['USE_SECONDARY_ALIGNMENTS'] not in ['True', 'true'] and read.is_secondary: continue
        if read.query_name not in read_mge_alignments:
            read_mge_alignments[read.query_name] = list()
        read_mge_alignments[read.query_name].append(((read.query_alignment_start, read.query_alignment_end),
                                                     (read.reference_start, read.reference_end),
                                                     read.reference_name))
    
    read_arg_alignments = dict()
    argSAM = pysam.AlignmentFile(config['INPUT']['ARGS_SAM_FILE'])
    for read in argSAM.fetch():
        if read.is_unmapped: continue
        if config['MISC']['USE_SECONDARY_ALIGNMENTS'] not in ['True', 'true'] and read.is_secondary: continue
        if read.query_name not in read_arg_alignments:
            read_arg_alignments[read.query_name] = list()
        read_arg_alignments[read.query_name].append(((read.query_alignment_start, read.query_alignment_end),
                                                     (read.reference_start, read.reference_end),
                                                     read.reference_name))

    overlaps = dict()
    for query_name in read_mge_alignments:
        for mge_alignment in read_mge_alignments[query_name]:
            for arg_alignment in read_arg_alignments[query_name]:
                if overlapsWith(mge_alignment[0], arg_alignment[0]):
                    if (mge_alignment[2], arg_alignment[2]) not in overlaps:
                        overlaps[(mge_alignment[2], arg_alignment[2])] = list()
                    overlaps[(mge_alignment[2], arg_alignment[2])].append((mge_alignment[1], arg_alignment[1]))
    
    overlapping_mges = list()
    for references in overlaps:
        if references[0] in overlapping_mges:
            continue
        arg_length = arg_genes_length[references[1]]
        mge_length = mge_genes_length[references[0]]
        if arg_length <= mge_length:
            overlapping_bases = [False for i in range(arg_length)]
            for reference_positions in overlaps[references]:
                for i in range(reference_positions[1][0], reference_positions[1][1]):
                    overlapping_bases[i] = True
            overlap_count = 0
            for base in overlapping_bases:
                if base: overlap_count += 1
            if overlap_count/arg_length >= overlap_threshold:
                overlapping_mges.append(references[0])
        else:
            overlapping_bases = [False for i in range(mge_length)]
            for reference_positions in overlaps[references]:
                for i in range(reference_positions[0][0], reference_positions[0][1]):
                    overlapping_bases[i] = False
            overlap_count = 0
            for base in overlapping_bases:
                if base: overlap_count += 1
            if overlap_count/mge_length >= overlap_threshold:
                overlapping_mges.append(references[0])

    with open(config['OUTPUT']['OUTPUT_PREFIX'] + config['EXTENSION']['OVERLAP'], 'w') as out_csv:
        csv_write = csv.writer(out_csv)
        csv_write.writerow(overlapping_mges)


def main():
    parser = argparse.ArgumentParser(description='Find MGEs that overlap alignments with ARGs')
    parser.add_argument('-r', help='Reads file', dest='reads_file', required=True)
    parser.add_argument('-a', help='ARGS Alignment file', dest='args_sam', required=True)
    parser.add_argument('-m', help='MGES Alignment file', dest='mges_sam', required=True)
    parser.add_argument('-o', help='Output Prefix', dest='out_prefix', required=True)
    parser.add_argument('-c', help='Config file', dest='config_path', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_path)

    logger = init_logger()

    config['INPUT'] = dict()
    config['INPUT']['ARGS_SAM_FILE'] = args.args_sam
    config['INPUT']['MGES_SAM_FILE'] = args.mges_sam
    config['INPUT']['INPUT_FILE_NAME_EXT'] = os.path.basename(args.reads_file)
    config['INPUT']['INPUT_FILE_NAME_NO_EXT'] = os.path.splitext(config['INPUT']['INPUT_FILE_NAME_EXT'])[0]
    config['INPUT']['INPUT_FILE_PATH'] = os.path.dirname(os.path.abspath(args.reads_file))
    config['INPUT']['INPUT_FILE'] = os.path.join(config['INPUT']['INPUT_FILE_PATH'], config['INPUT']['INPUT_FILE_NAME_EXT'])

    config['OUTPUT'] = dict()
    config['OUTPUT']['OUTPUT_PREFIX'] = args.out_prefix
    config['OUTPUT']['OUT_DIR'] = os.path.dirname(os.path.abspath(config['OUTPUT']['OUTPUT_PREFIX']))

    overlap_finder(config)

if __name__ == "__main__":
    main()