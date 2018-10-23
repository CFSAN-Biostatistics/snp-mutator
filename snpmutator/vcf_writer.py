"""
VCF file writer.
"""

from __future__ import print_function
from __future__ import absolute_import

from collections import defaultdict
import datetime
import sys

from snpmutator import __version__


class VcfWriter(object):
    """Class to write VCF files for the mutated replicates.
    """

    def __init__(self, original_seq, file_path):
        """Initialize the VCF writer.

        Parameters
        ----------
        original_seq : str
            Original sequence string with no mutations.
        file_path : str
            Path to output vcf file.
        """
        self.reference = original_seq
        self.file_path = file_path
        self.replicate_names = list()
        self.known_alts = defaultdict(list) # key=pos value=list of alts
        self.alt_dict = dict()              # key=(pos, replicate) value=alt num

    def store_replicate_mutations(self, replicate_name, new_indexed_seq, subs_positions, insertion_positions, deletion_positions):
        """Store all of the mutations for a specified replicate.

        Parameters
        ----------
        replicate_name : str
            Name of sythetic replicate.
        new_indexed_seq : list of str
            List indexed by original position containing strings of bases.  Each
            string can be 0 - 2 bases long, where a zero length string indicates
            a deletion at the position, a string containing a single base indicates
            a snp at the position, and a string containing two bases
            indicates an insertion prior to the original base at the position.
        subs_positions : list of int
            List of positions where substitions are introduced.
        insertion_positions : list of int
            List of positions where insertions are introduced.
        deletion_positions : list of int
            List of positions where bases are deleted.
        """
        self.replicate_names.append(replicate_name)

        # Store the snps
        for pos in subs_positions:
            alt = new_indexed_seq[pos]
            if alt not in self.known_alts[pos]:
                self.known_alts[pos].append(alt)
            alt_num = 1 + self.known_alts[pos].index(alt) # ALT zero is reference match
            self.alt_dict[(pos, replicate_name)] = alt_num

        # Store the insertions
        for pos in insertion_positions:
            alt = new_indexed_seq[pos]
            if alt not in self.known_alts[pos]:
                self.known_alts[pos].append(alt)
            alt_num = 1 + self.known_alts[pos].index(alt) # ALT zero is reference match
            self.alt_dict[(pos, replicate_name)] = alt_num

        # Store the deletions
        for pos in deletion_positions:
            alt = "*"
            if alt not in self.known_alts[pos]:
                self.known_alts[pos].append(alt)
            alt_num = 1 + self.known_alts[pos].index(alt) # ALT zero is reference match
            self.alt_dict[(pos, replicate_name)] = alt_num

    def write(self, reference_name, chrom):
        """Write the stored mutations to the VCF output file.

        Parameters
        ----------
        reference_name : str
            Name of reference to write to the header.
        chrom : str
            Name of the chrom to be written in the VCF CHROM column. It will be the same for
            all rows because there is only one contig.
        """
        print("Writing VCF file.", file=sys.stderr)
        with open(self.file_path, "w") as f:
            # Write the VCF file header sections.
            f.write('##fileformat=VCFv4.2\n')
            f.write(datetime.datetime.strftime(datetime.datetime.now(), '##fileDate=%Y%m%d\n'))
            f.write('##source=snpmutator %s\n' % __version__)
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##reference=%s\n' % reference_name)
            header = '#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT %s\n' % ' '.join(self.replicate_names)
            header = header.replace(' ', '\t')
            f.write(header)

            # Write the variant calls
            for pos in sorted(self.known_alts):
                ref = self.reference[pos]
                alt_str = ','.join(self.known_alts[pos])
                fields = [chrom, str(pos + 1), '.', ref, alt_str, '.', '.', '.', "GT"]
                # default to ref (0) if alt_dict has no variant for this (pos, replicate_name)
                genotypes = [str(self.alt_dict.get((pos, replicate_name), 0)) for replicate_name in self.replicate_names]
                fields.extend(genotypes)
                f.write('\t'.join(fields))
                f.write("\n")
