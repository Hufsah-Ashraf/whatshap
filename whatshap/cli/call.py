"""
Call variants
"""
import logging
import pysam
import pyfaidx
from collections import defaultdict, namedtuple, deque
from contextlib import ExitStack
from whatshap.core import Caller
#from whatshap._call import hash_to_dna
#from whatshap.align import edit_distance
#import cProfile
#import re
#import timeit

logger = logging.getLogger(__name__)


def add_arguments(parser):
	add = parser.add_argument
	add('reference', metavar='FASTA', help='Reference genome FASTA')
	add('bam', metavar='BAM', help='BAM file ')
	add('vcf', metavar='VCF', help='VCF file ')
	add('kmer_size', metavar='KMER', help='kmer_size')
	add('window_size', metavar='WINDOW', help='length of the window in one direction i.e. half the total window size')


def run_call(reference, bam, vcf, kmer_size, window_size):
    fasta = pyfaidx.Fasta(reference, as_raw=True)
    samfile = pysam.Samfile(bam)
    variantslist=list()
    reader= open(vcf, 'r')
    call=0
    for line in reader:
        line_r = line.strip().split("\t")
        varpos= int(line_r[1])
        variantslist.append(varpos)
    variant=0
    encoded_references={}
    k = int(kmer_size)
    window= int(window_size)
    min_support = 2
    chromosome = None
    hashed={}

    for line in samfile:
        firstrefpos=line.reference_start
        break
    for i, bam_alignment in enumerate(samfile):
        if not bam_alignment.is_unmapped:
            if bam_alignment.reference_name != chromosome:
                chromosome = bam_alignment.reference_name
                if chromosome in encoded_references:
                    caller = Caller(encoded_references[chromosome], k)
                else:
                    ref = fasta[chromosome]
                    #caller = Caller(str(ref).encode('UTF-8'), k)
                    encoded_references[chromosome]= str(ref).encode('UTF-8')

                    caller = Caller(encoded_references[chromosome], k)
            if call==0:
                caller.all_variants(variantslist)
                call=1
            else:
                pass
            caller.add_read(bam_alignment.pos, bam_alignment.cigartuples, str(bam_alignment.query).encode('UTF-8'))
            #for ref_pos, ref_kmer, pileup_kmers in caller.process_complete_columns():
                #variantposition=variantslist[variant]
                #varstart= variantposition-window
                #varend= variantposition+window


                #if ref_pos>=varstart and ref_pos<=varend:
                    #pass


                #elif variant<len(variantslist)-1 and ref_pos>= (variantslist[variant+1]-window) and ref_pos<= (variantslist[variant+1]+window):
                        #variant+=1


                #else:

                    #if none_key not in pileup_kmers and (len(pileup_kmers) > 0):


                        #for kmer,count in pileup_kmers.items():
                            #if ref_kmer not in hashed and kmer not in hashed:
                             #   hashed_ref_kmer= caller.hash_to_dna(ref_kmer,k)
                              #  hashed_kmer= caller.hash_to_dna(kmer,k)
                               # hashed[ref_kmer]=hashed_ref_kmer
                                #hashed[kmer]=hashed_kmer
                            #else:
                             #   if ref_kmer not in hashed:
                              #      hashed_ref_kmer= caller.hash_to_dna(ref_kmer,k)
                               #     hashed_kmer= hashed[kmer]
                                #    hashed[ref_kmer]=hashed_ref_kmer
                                #elif kmer not in hashed:
                                 #   hashed_kmer= caller.hash_to_dna(kmer,k)
                                  #  hashed_ref_kmer= hashed[ref_kmer]
                                   # hashed[kmer]=hashed_kmer
                                #else:
                                 #   hashed_kmer=hashed[kmer]
                                  #  hashed_ref_kmer= hashed[ref_kmer]
                            #writer.write (str(ref_pos)+'\t'+str(hashed_ref_kmer)+'\t'+str(hashed_kmer)+'\t'+str(count)+'\n')











def main(args):
	#cProfile.runctx('run_call(**vars(args))', {'run_call': run_call, 'args':args}, {})
	#code_to_test = """
	run_call(**vars(args))
	#"""
	#global argms
	#argms= args
	#elapsed_time = timeit.timeit(stmt= "run_call(**vars(argms))",globals=globals(), number=100)/100
	#print(elapsed_time)
