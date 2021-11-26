# cython: profile=True
from libcpp.string cimport string
cimport cython
#from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
from libcpp.vector cimport vector
from libcpp.deque cimport deque
import collections
cimport cpp
from distutils.core import setup
from cpython.ref cimport PyObject
from cython.operator import address

def enumerate_reference_kmers(string reference, int k):
        cdef int A = ord('A')
        cdef int C = ord('C')
        cdef int G = ord('G')
        cdef int T = ord('T')
        cdef c = 0
        cdef int h = 0
        cdef int mask = (1 << (2*k)) - 1
        cdef int i = 0
        cdef enum_refkmers= collections.deque()
        for i in range(len(reference)):
                c = reference[i]
                if c == A:
                        h = ((h << 2) | 0) & mask
                elif c == C:
                        h = ((h << 2) | 1) & mask
                elif c == G:
                        h = ((h << 2) | 2) & mask
                elif c == T:
                        h = ((h << 2) | 3) & mask
                else:
                        h = ((h << 2) | 0) & mask
                if i >= k-1:
                        #enum_refkmers.append((h, i+1))
                        yield(h, i+1)
        




def enumerate_kmers(int pos, string query_string, int k, vector[vector[int]] cigartuples, str query_name):
        '''
        Generates all kmers in the read and yields pairs (kmer_hash, position),
        where kmer_hash is a binary representation of the kmer and postion is
        the reference position this kmer has been aligned to.
        '''
        cdef h = 0
        cdef int n_kmer=0
        cdef int mask = (1 << (2*k)) - 1
        cdef int cigar_index = 0
        cdef int cigar_op = cigartuples[cigar_index][0]
        cdef int cigar_length = cigartuples[cigar_index][1]
        #cdef int cigar_length = cigartuples[cigar_index][1]
        cdef int BAM_CMATCH = 0     # M
        cdef int BAM_CINS = 1       # I
        cdef int BAM_CDEL = 2       # D
        cdef int BAM_CREF_SKIP = 3  # N
        cdef int BAM_CSOFT_CLIP = 4 # S
        cdef int BAM_CHARD_CLIP = 5 # H
        cdef int BAM_CPAD = 6       # P
        cdef int BAM_CEQUAL = 7     # =
        cdef int BAM_CDIFF = 8      # X
        cdef int BAM_CBACK = 9      # B
        cdef int A = ord('A')
        cdef int C = ord('C')
        cdef int G = ord('G')
        cdef int T = ord('T')
        cdef int i = 0
        cdef int consecutive = 0
        cdef enum_kmers = collections.deque()
        while i < len(query_string):
                #print('query name is ', query_name)
                #print('i ', i, ' consecutive ', consecutive, ' cigar_index ', cigar_index)
                # process cigar entries that don't consume a character from the read
                while True:
                        if (cigar_op == BAM_CDEL) or (cigar_op == BAM_CREF_SKIP):
                                #print(pos, '-->', pos + cigar_length)
                                pos += cigar_length
                                #print("enum_kmer_pos ",pos)
                        elif cigar_op == BAM_CSOFT_CLIP:
                                #i += cigar_length
                                consecutive = 0
                                #print("consecutive set to 0")
                        elif (cigar_length == 0) or (cigar_op == BAM_CHARD_CLIP):
                                pass
                        else:
                                break
                        cigar_index += 1
                        cigar_op, cigar_length = cigartuples[cigar_index]
                        #print("cigar_op ", cigar_op , " cigar_length ", cigar_length)
                if i >= len(query_string):
                        break

                # update hash
        
                if query_string[i] == A:
                        h = ((h << 2) | 0) & mask
                        consecutive += 1
                        #print(" first if, h and consecutive is ", h, " " ,consecutive, "query_string of ", i," is ",query_string[i])
                elif query_string[i] == C:
                        h = ((h << 2) | 1) & mask
                        consecutive += 1
                        #print(" first elseif, h and consecutive is ", h, " " ,consecutive, "query_string of ", i," is ",query_string[i])
                elif query_string[i] == G:
                        h = ((h << 2) | 2) & mask
                        consecutive += 1
                        #print(" second elseif, h and consecutive is ", h, " " ,consecutive, "query_string of ", i," is ",query_string[i])
                elif query_string[i] == T:
                        h = ((h << 2) | 3) & mask
                        consecutive += 1
                        #print(" 3rd else if, h and consecutive is ", h, " " ,consecutive, "query_string of ", i," is ",query_string[i])
                else:
                        n_kmer=1
                        consecutive += 1
                        #print( "N encountered at", i)
                        #print(" last else, h and consecutive is ", h, " " ,consecutive, "query_string of ", i," is ",query_string[i])
                        
                if consecutive >= k:
                        if n_kmer ==0:
                                #print ('enum_kmers yielded the following')
                                #print(h , pos+1)
                                #enum_kmers.append((h,pos+1))
                                yield (h, pos+1)
                        elif n_kmer==1:
                                #print (h, pos+1, "this kmer had N")
                                n_kmer=0
                
                # consume one character of read
                assert cigar_length > 0
                if (cigar_op == BAM_CMATCH) or (cigar_op == BAM_CEQUAL) or (cigar_op == BAM_CDIFF):
                        cigar_length -= 1
                        pos += 1
                elif cigar_op == BAM_CINS:
                        cigar_length -= 1
                else:
                        assert False, 'Unexpected cigar operation'

                i += 1
        


                

def hash_to_dna(int h, int k):
        cdef list l = [None]*k
        cdef int i = 0
        for i in range(k):
                if h&3 == 0:
                        l[k-1-i]='A'   
                elif h&3 == 1:
                        l[k-1-i]='C'   
                elif h&3 == 2:
                        l[k-1-i]='G'
                elif h&3 == 3:
                        l[k-1-i]='T'
                
                h = h >> 2
        return ''.join(l)

class Caller:
        def __init__(self, string reference, int k):
                self.k = k
                # bam record for each active read
                self.bam_records = collections.deque()
                # k-mer generators for each active read
                self.kmer_generators =collections.deque()
                self.kmer_generators_finished = collections.deque()
                # .. corresponding queue with the latest (kmer, position) for each read
                self.current_kmers = collections.deque() 
                self.ref_kmer_generator = enumerate_reference_kmers(reference,k)
                # position that every other operation is relative to
                self.pileup_columns=collections.deque()
                self.ref_kmers = collections.deque()
                kmer, pos = self.ref_kmer_generator.__next__()
                #print ('initial kmer', kmer)
                #print ('initial pos', pos)
                self.pileup_columns.append(collections.defaultdict(int))
                self.ref_kmers.append(int(kmer))
                self.ref_pos = pos

        def add_read(self, bam_alignment):
                #print ("add_read_start")
                self.bam_records.append(bam_alignment)
                self.kmer_generators.append(enumerate_kmers(bam_alignment.pos,str(bam_alignment.query).encode('UTF-8'), self.k, bam_alignment.cigartuples, str(bam_alignment.query_name)))
                self.kmer_generators_finished.append(False)
                #print(' -1 called')
                kmer, pos = self.kmer_generators[-1].__next__()
                #print ('add_read_kmer', kmer)
                #print('add_read_pos',pos)
                self.current_kmers.append((kmer,pos))
                #print('self.current_kmers', self.current_kmers)
                pileup_column, ref_kmer = self.get_column(pos)
                #print ("within add_read get_column is completed successfully")
                pileup_column[kmer] += 1
                #print(self.current_kmers)
                #print ("add_read_end")
        def finish(self):
                pass

        def get_column(self, pos):
                #print ('get_column_start with pos being', pos)
                #print ('get_column_start with ref_pos being', self.ref_pos)
                index = pos - self.ref_pos
                if index >=0:
                        #print ('entered if')
                        #print("index before entering while ", index)
                        #print("size_pileup_columns before entering while ", len(self.pileup_columns))
                        while len(self.pileup_columns) <= index:
                                #print ('entered while')
                                kmer, pos = self.ref_kmer_generator.__next__()
                                #print("get_column_kmer ",kmer)
                                #print("get_column_pos ",pos)
                                self.ref_kmers.append(int(kmer))
                                self.pileup_columns.append(collections.defaultdict(int))
                                #print ('length of pileup_coulmns updated to ', len(self.pileup_columns))
                        #print ('exited while')
                        return self.pileup_columns[index], self.ref_kmers[index]
                else:
                        #print ('entered else')
                        #print ('entered else with length of pileup_coulmns being ', len(self.pileup_columns))
                        kmer, pos = self.ref_kmer_generator.__next__()
                        self.ref_kmers.append(kmer)
                        self.pileup_columns.append(collections.defaultdict(int))
                        #print ('length of pileup_coulmns updated to ', len(self.pileup_columns))
						
                        return collections.defaultdict(int), kmer
        def pop_column(self):
                if len(self.pileup_columns) > 0:
                        result = (self.ref_pos, self.ref_kmers.popleft(), self.pileup_columns.popleft())
                else:
                        kmer, pos = self.ref_kmer_generator.__next__()
                        assert pos == self.ref_pos
                        result = (self.ref_pos, kmer, None)
                self.ref_pos += 1
                return result
                
        def process_complete_columns(self):
                '''
                Perform calling of columns that are complete, i.e. they cannot receive
                more reads because subsequent reads are further to the right.
                '''
                # compute largest position for which k-mer pileup can safely be generated
                #print ('process_complete_columns start')
                target_pos = self.bam_records[-1].pos + self.k - 1
                #print ("bam_record ",self.bam_records[-1].pos, "k ", self.k)
                self.advance_to(target_pos)
                cdef list complete_columns= []
                while self.ref_pos < target_pos:
                        #print ('entered while')
                        complete_columns.append(self.pop_column())
                #print ('process_complete_columns end with length of complete_columns being ', len(complete_columns))
                return complete_columns
        
        def advance_to(self, int target_pos):
                '''
                Add all k-mer from all reads up to target_pos to pileup_columns.
                '''
                #print ('advance_to started with target_pos ',target_pos)
                #print('advance_to', target_pos)
                cdef int i =0
                cdef int kmer=0
                cdef int pos=0
                for i, kmer_generator in enumerate(self.kmer_generators):
                        #print ('length of self.kmer_generators is ', len(self.kmer_generators))
                        #print ('entered for loop with i being ',i)
                        try:
                                
                                kmer, pos = self.current_kmers[i]
                                #print((kmer,pos) , "came out of current kmers of ", i)
                                #print("within try and outside while kmer is ",kmer)
                                #print("within try and outside while pos is ",pos)
                                while pos < target_pos:
                                        #print(' normally called')
                                        kmer, pos = kmer_generator.__next__()
                                        #print("within while kmer is ",kmer)
                                        #print("within while pos is ",pos)
                                        #print('  read', i, 'kmer', hash_to_dna(kmer, self.k), 'pos', pos)
                                        #print ('inside while get_column called')
                                        pileup_column, ref_kmer = self.get_column(pos)
                                        #print ('inside while get_column completed successfully')
                                        pileup_column[kmer] += 1
                                self.current_kmers[i] = (kmer,pos)
                                #print ((kmer,pos), 'was added to current_kmers of ', i)
                        except StopIteration:
                                #print('entered except')
                                self.kmer_generators_finished[i] = True
                                pass
                        #print('self.current_kmers.left ', self.current_kmers[-1], 'self.current_kmers.right ', self.current_kmers[0] )
                while (len(self.kmer_generators) > 0) and self.kmer_generators_finished[0]:
                        #print ('current kmers popping ',self.current_kmers.popleft())
                        self.bam_records.popleft()
                        #print(' kmer_generators popping ',self.kmer_generators.popleft())
                        self.kmer_generators_finished.popleft()
                #print("at advance to end kmer is ", kmer)
                #print("at advance to end pos is ", pos)
                #print('advance_to end')
                
                
                







