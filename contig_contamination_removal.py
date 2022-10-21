#

from Bio import SeqIO
import argparse
import sys
import os

def filterBlastn6(blastn, refList, fasta, pident, qcov, path):
    pident = float (pident)
    qcov = float (qcov)
    name1 = blastn.rsplit('.',1)[0].split('/')[-1]
    name2 = fasta.rsplit('.',1)[0].split('/')[-1]
    path=os.path.abspath(path)
    refItems=[]
    promisedCtgs=[]
    questionedCtgs=[]
    contigs=[]
    with open (refList, 'r') as refList:
        for item in refList:
            item = item.strip().upper()
            refItems.append(item)

    with open ('%s/%s.filtered.blastn'%(path, name1), 'w') as filtered:
        with open ('%s/%s.cont.blastn'%(path, name1), 'w') as cont:
            with open (blastn, 'r') as blastn:
                for line in blastn:
                    line = line.strip()
                    Line = line.split('\t')
                    sp = Line[13].split()[0].upper()
                    ctg = Line[0]
                    if ctg not in contigs:
                        contigs.append(ctg)
                    if sp == "PREDICTED:":
                        sp = Line[13].split()[1].upper()
                    if sp in refItems:
                        if ctg not in promisedCtgs:
                            promisedCtgs.append(ctg)
                        filtered.write('%s\n'%line)
                    else:
                        cont.write('%s\n'%line)
                        if ctg in promisedCtgs:
                            pid = float(Line[2])
                            length = float(Line[3])
                            qlen=float(Line[4])
                            if pid >= pident and length/qlen >= qcov:
                                promisedCtgs.pop(ctg)
                                if ctg not in  questionedCtgs:
                                    questionedCtgs.append(ctg)
    with open ('%s/%s.promised.fasta'%(path, name1), 'w') as promised:
        with open ('%s/%s.questioned.fasta'%(path, name1), 'w') as questioned:
            with open ('%s/%s.no.blastHits.fasta'%(path, name1), 'w') as noHits:
                with open (fasta, 'r') as fasta:
                    for seqs in SeqIO.parse(fasta, 'fasta'):
                        ID = seqs.id
                        seq = seqs.seq
                        if ID in promisedCtgs:
                            promised.write('>%s\n'%ID)
                            promised.write('%s\n'%seq)
                        if ID in questionedCtgs:
                            questioned.write('>%s\n'%ID)
                            questioned.write('%s\n'%seq)
                        if ID not in contigs:
                            noHits.write('>%s\n'%ID)
                            noHits.write('%s\n'%seq)

if __name__=="__main__":
    """
    Parse the function and usage.
    """
    parser = argparse.ArgumentParser(description='Filter blastn result using a ref species list')
    parser.add_argument('-v', '--version', action='version', version='1.0')
    parser.add_argument('-b', dest='blastn', help='blastn result with customised format 6', type = str)
    parser.add_argument('-r', dest='refList', help='a list contains all species from the kingdom, such as plantae', type = str)
    parser.add_argument('-f', dest='fasta', help='the query fasta file', type=str)
    parser.add_argument('-p', dest='pid', help='min percentage of identity in blastn. Default: 80.', default = 80, type=float)
    parser.add_argument('-c', dest='qcov', help='min percentage of the query coverage. Default: 50.', default = 50, type=float)
    parser.add_argument('-o', dest='output', help='the output directory', type=str)
    args = parser.parse_args()
    if None not in [args.blastn, args.refList, args.fasta, args.output]:
        filterBlastn6(args.blastn, args.refList, args.fasta, args.pid, args.qcov, args.output)
    else:
        print
        parser.print_help()
        print
        sys.exit(1)
