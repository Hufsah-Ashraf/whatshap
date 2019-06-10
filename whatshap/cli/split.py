"""
Split reads by haplotype.

Reads FASTQ/BAM file and a list of haplotype assignments (such as generated by 
whatshap haplotag --output-haplotag-list). Outputs one FASTQ/BAM per haplotype.
BAM mode is intended for unmapped BAMs (such as provided by PacBio).
"""
import logging
import sys
import gzip
import pysam
from collections import defaultdict
from subprocess import Popen, PIPE

from contextlib import ExitStack
from whatshap import __version__
from whatshap.timer import StageTimer

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--output-h1', default=None,
		help='Output file to write reads from Haplotype 1 to. Use ending .gz to '
		'create gzipped file.')
	arg('--output-h2', default=None,
		help='Output file to write reads from Haplotype 2 to. Use ending .gz to '
		'create gzipped file.')
	arg('--output-untagged', default=None,
		help='Output file to write untagged reads to. Use ending .gz to '
		'create gzipped file.')
	arg('--add-untagged', default=False, action='store_true',
		help='Add reads without tag to both H1 and H2 output streams.')
	arg('--pigz', default=False, action='store_true',
		help='Use the pigz program for gzipping output.')
	arg('--only-largest-block', default=False, action='store_true',
		help='Only consider reads to be tagged if they belong to the largest '
		'phased block (in terms of read count) on their respective chromosome')
	arg('--read-lengths-histogram', default=None,
		help='Output file to write read lengths histogram to in tab separated format.')
	arg('reads_file', metavar='READS', help='Input FASTQ/BAM file with reads (fastq can be gzipped)')
	arg('list_file', metavar='LIST',
		help='Tab-separated list with (at least) two columns <readname>,<haplotype> (can be gzipped)')


def validate(args, parser):
	if (args.output_h1 is None) and (args.output_h2 is None) and (args.output_untagged is None):
		parser.error('Nothing to be done since neither --output-h1 nor --output-h2 nor --output-untagged are given.')


def open_possibly_gzipped(filename, readwrite='r', pigz=False):
	if filename is None:
		return None
	if readwrite == 'r':
		if filename.endswith('.gz'):
			return gzip.open(filename, 'rt')
		else:
			return open(filename)
	elif readwrite == 'w':
		if filename.endswith('.gz'):
			if pigz:
				g = Popen(['pigz'], stdout=open(filename, 'w'), stdin=PIPE)
				return g.stdin
			else:
				return gzip.open(filename, 'w')
		else:
			return open(filename, 'wb')
	else:
		assert False, 'Invalid open mode'


def read_fastq(filename):
	'''Yields pairs (readname, record) where record is a list of four lines.'''
	f = open_possibly_gzipped(filename)
	n = 0
	while True:
		record = [ f.readline() for _ in range(4) ]
		if record[3] == '':
			break
		assert record[0].startswith('@'), record
		assert record[2].startswith('+'), record
		name = record[0][1:].split()[0]
		yield name, record


def run_split(
		reads_file,
		list_file,
		output_h1=None,
		output_h2=None,
		output_untagged=None,
		add_untagged=False,
		pigz=False,
		only_largest_block=False,
		read_lengths_histogram=None,
	):

	timers = StageTimer()
	timers.start('overall')

	with ExitStack() as stack:

		HAPLOTYPE_TO_INT = {'none': 0, 'H1': 1, 'H2': 2}
		largest_block_map = None
		if only_largest_block:
			# do one first pass and determine largest blocks
			with open_possibly_gzipped(list_file) as f:
				block_sizes = defaultdict(int)
				logger.info('Reading %s to determine block sizes', list_file)
				for line in f:
					if line.startswith('#'):
						continue
					fields = line.split()
					assert len(fields) == 4, 'Error parsing input file "{}"'.format(list_file)
					readname, haplotype_name, phaseset, chromosome = fields
					if HAPLOTYPE_TO_INT[haplotype_name] != 0:
						block_sizes[(chromosome,phaseset)] += 1
			largest_block_map = dict()
			chromosomes = set(chromosome for chromosome,phaseset in block_sizes.keys())
			logger.info('Largest blocks per chromosome:')
			for chromosome in sorted(chromosomes):
				blocks = sorted(((size, phaseset) for (_chromosome,phaseset), size in block_sizes.items() if _chromosome == chromosome), reverse=True)
				largest_block_size, largest_block = blocks[0]
				largest_block_map[chromosome] = largest_block
				logger.info('   %s: phaseset "%s" with %d tagged reads', chromosome, largest_block, largest_block_size)

		# mapping of read names to haplotypes, 0 means untagged
		haplotype = defaultdict(int)
		tags_removed = 0
		assigned_to_haplotype = 0
		with open_possibly_gzipped(list_file) as f:
			logger.info('Reading %s', list_file)
			for line in f:
				if line.startswith('#'):
					continue
				fields = line.split()
				assert len(fields) == 4, 'Error parsing input file "{}"'.format(list_file)
				readname, haplotype_name, phaseset, chromosome = fields
				haplotype_int = HAPLOTYPE_TO_INT[haplotype_name]
				if only_largest_block:
					if (haplotype_int != 0) and (largest_block_map[chromosome] != phaseset):
						tags_removed += 1
						haplotype_int = 0
				haplotype[readname] = haplotype_int
				if haplotype_int != 0:
					assigned_to_haplotype += 1
		logger.info('... read %d records:', len(haplotype))
		logger.info('...   %d reads are assigned to a haplotype.', assigned_to_haplotype)
		if only_largest_block:
			logger.info('...   %d reads have been unassigned because they are not from the largest block.', tags_removed)

		output_h1_file = None
		output_h2_file = None
		output_untagged_file = None

		read_lengths_histogram_dict = dict()

		n_reads = 0
		# TODO: Avoid code duplication in the two code blocks
		# TODO: Detect file type in a smarter way, not based on ending
		if reads_file.endswith('.bam'):
			bamreader = pysam.AlignmentFile(reads_file, mode='rb', check_sq=False)
			if output_h1 is not None:
				output_h1_file = pysam.AlignmentFile(output_h1, 'wb', template=bamreader)
			if output_h2 is not None:
				output_h2_file = pysam.AlignmentFile(output_h2, 'wb', template=bamreader)
			if output_untagged is not None:
				output_untagged_file = pysam.AlignmentFile(output_untagged, 'wb', template=bamreader)

			for record in bamreader:
				n_reads += 1
				h = haplotype[record.query_name]

				if read_lengths_histogram is not None:
					read_length = len(record.query_sequence)
					if read_length not in read_lengths_histogram_dict:
						read_lengths_histogram_dict[read_length] = [0,0,0]
					read_lengths_histogram_dict[read_length][h] += 1

				if output_h1_file is not None:
					if (h==1) or (h==0 and add_untagged):
						output_h1_file.write(record)

				if output_h2_file is not None:
					if (h==2) or (h==0 and add_untagged):
						output_h2_file.write(record)

				if output_untagged_file is not None:
					if h==0:
						output_untagged_file.write(record)

		else:
			output_h1_file = open_possibly_gzipped(output_h1, 'w', pigz)
			output_h2_file = open_possibly_gzipped(output_h2, 'w', pigz)
			output_untagged_file = open_possibly_gzipped(output_untagged, 'w', pigz)

			for name, record in read_fastq(reads_file):
				n_reads += 1
				h = haplotype[name]

				if read_lengths_histogram is not None:
					read_length = len(record[1].strip())
					if read_length not in read_lengths_histogram_dict:
						read_lengths_histogram_dict[read_length] = [0,0,0]
					read_lengths_histogram_dict[read_length][h] += 1

				if output_h1_file is not None:
					if (h==1) or (h==0 and add_untagged):
						for line in record:
							output_h1_file.write(line.encode('utf-8'))

				if output_h2_file is not None:
					if (h==2) or (h==0 and add_untagged):
						for line in record:
							output_h2_file.write(line.encode('utf-8'))

				if output_untagged_file is not None:
					if h==0:
						for line in record:
							output_untagged_file.write(line.encode('utf-8'))

		if read_lengths_histogram is not None:
			with open(read_lengths_histogram, 'wt') as t:
				print('#length', 'count-untagged', 'count-h1', 'count-h2', sep='\t', file=t)
				for length in sorted(read_lengths_histogram_dict.keys()):
					print(length, *(read_lengths_histogram_dict[length]), sep='\t', file=t)
				t.close()

		if output_h1_file is not None:
			output_h1_file.close()
		if output_h2_file is not None:
			output_h2_file.close()
		if output_untagged_file is not None:
			output_untagged_file.close()

	logger.info('\n== SUMMARY ==')
	logger.info('Total reads processed:              %12d', n_reads)


def main(args):
	run_split(**vars(args))
