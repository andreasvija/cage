"""
tss_id
gene_name
chr
peak_start
peak_end
strand
->
seqname - chr
source - fantom5
feature - promotor
start - peak_start
end - peak_end
score - .
strand - strand
frame - .
attribute:
gene_id gene_name; // meta-feature
promotor_id tss_id; // feature
"""

inputfilename = 'FANTOM5_promoter_annotations.tsv'
outputfilename = 'FANTOM5_promoter_annotations.gtf'

with open(inputfilename, 'r') as infile:
	with open(outputfilename, 'w') as outfile:

		infile.readline()
		for inline in infile:

			incols = inline.strip().split('\t')
			outline = ''

			outline += incols[2] + '\t' # seqname - chr
			outline += 'fantom5' + '\t' # source - fantom5
			outline += 'promotor' + '\t' # feature - promotor
			outline += incols[3] + '\t' # start - peak_start
			outline += incols[4] + '\t' # end - peak_end
			outline += '+' + '\t' # score - .
			outline += incols[5] + '\t' # strand - strand
			outline += '.' + '\t' # frame - .

			outline += 'gene_id' + ' "' + incols[0] + '";' # promotor_id tss_id; // feature
			outline += 'geeni_id' + ' "' + incols[1] + '";\t' # gene_id gene_name; // meta-feature
			outline += '\n'

			outfile.write(outline)
