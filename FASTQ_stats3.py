### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.25 ###


__usage__ = """
		python3 FASTQ_stats3.py
		--in_file <FULL_PATH_TO_FASTQ_FILE> |	--in_dir <FULL_PATH_TO_DIRECTORY>
		
		bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
		"""

import sys, gzip, glob

# --- end of imports --- #

def calculate_n50( total_length ):
	"""! @brief calculate N50 based on list of given read lengths """
	
	total = sum( total_length )
	half = total / 2.0
	sorted_lengths = sorted( total_length )[::-1]
	counter = 0
	i = 0
	for i in range( len( sorted_lengths ) ):
		counter += sorted_lengths[ i ]
		if counter >= half:
			return sorted_lengths[ i ]
	sys.exit( "ERROR: no read lengths detected." )


def analyze_FASTQ( filename ):
	"""! @brief analysis of FASTQ file """
	
	gzip_state = False
	try:
		file_extension = filename.split('.')[-1]
		if file_extension in [ "gz", "gzip", "GZ", "GZIP" ]:
			gzip_state = True
	except:
		pass
	
	if not gzip_state:
		with open( filename, "r" ) as f:
			total_length = []
			total_GC = []
			line = f.readline()	#header
			while line:
				seq = str( f.readline().strip().upper() )
				total_length.append( len( seq ) )
				total_GC.append( seq.count('C') )
				total_GC.append( seq.count('G') )
				f.readline()	#useless line
				f.readline()	#quality line
				line = f.readline()
			print( "number of reads: " + str( len( total_length ) ) )
			total_len = sum( total_length )
			total_gc = sum( total_GC )
			n50 = calculate_n50( total_length )
			print( filename )
			print( "total number of nucleotides:\t" + str( total_len ) )
			print( "average read length:\t" + str( total_len / float( len( total_length ) ) ) )
			print( "GC content:\t" + str( total_gc / float( total_len ) ) )
			print( "N50: " + str( n50 ) )
		
	else:
		with gzip.open( filename, "rb" ) as f:
			total_length = []
			total_GC = []
			line = f.readline()	#header
			while line:
				seq = str( f.readline().strip().upper() )
				total_length.append( len( seq ) )
				total_GC.append( seq.count('C') )
				total_GC.append( seq.count('G') )
				f.readline()	#useless line
				f.readline()	#quality line
				line = f.readline()
			print( "number of reads: " + str( len( total_length ) ) )
			total_len = sum( total_length )
			total_gc = sum( total_GC )
			n50 = calculate_n50( total_length )
			print( filename )
			print( "total number of nucleotides:\t" + str( total_len ) )
			print( "average read length:\t" + str( total_len / float( len( total_length ) ) ) )
			print( "GC content:\t" + str( total_gc / float( total_len ) ) )
			print( "N50: " + str( n50 ) )


def main( arguments ):
	"""! @brief runs everything """
	
	if '--in_file' in arguments:
		input_file = arguments[ arguments.index( '--in_file' )+1 ]
		analyze_FASTQ( input_file )
	else:
		directory = arguments[ arguments.index( '--in_dir' )+1 ]
		if directory[-1] != '/':
			directory += "/"
		input_files = []
		extensions = [ ".fq", ".fastq", ".fq.gzip", ".fastq.gzip", ".fq.gz", ".fastq.gz", ".FQ", ".FASTQ", ".FQ.GZIP", ".FASTQ.GZIP", ".FQ.GZ", ".FASTQ.GZ" ]
		for extension in extensions:
			input_files += glob.glob( directory + '*' + extension )
		for filename in input_files:
			try:
				analyze_FASTQ( filename )
			except:
				print( "ERROR while processing " + filename )


if __name__ == '__main__':
	
	if '--in_file' in sys.argv or '--in_dir' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
