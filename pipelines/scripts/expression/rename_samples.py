import sys

def rename_samples(path_to_matrix, path_to_bams):
    table = {}
    with open(path_to_bams, 'r') as infile:
        for i, line in enumerate(infile):
            if i > 0:
                splitted = line.strip().split(',')
                sample   = splitted[0]
                filename = splitted[1].split('/')[-1].rstrip()
                table[filename] = sample

    with open(path_to_matrix, 'r') as infile:
        with open('counts.matrix.renamed.txt', 'w') as outfile:
            for i, line in enumerate(infile):
                if i == 0:
                    header = line.strip('\n').split('\t')
                    colnames = [header[0]]+[table[i] for i in header[1:]]
                    outfile.write('\t'.join(colnames)+'\n')
                else:
                    outfile.write(line)

if __name__ == '__main__':
    path_to_matrix = sys.argv[1]
    path_to_bams = sys.argv[2]
    rename_samples(path_to_matrix, path_to_bams)
