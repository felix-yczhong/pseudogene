
import numpy as np

from smaca.bam import Bam

from pseudogene.utils import BASE_TO_INDEX_MAPPING, BASE_MATCH, BASES

class PseudoBam(Bam):
    # 0-based coordinate
    def get_cov(self, ranges):
        c, start, stop = ranges
        cov_details = np.zeros((stop - start, len(BASES) + 1), dtype=np.int32)

        for read in self.samfile.fetch(c, start, stop, until_eof=True):
            if not read.is_duplicate and not read.is_qcfail:
                for seq_pos, genomic_pos in read.get_aligned_pairs(matches_only=False, with_seq=False):
                    if genomic_pos is not None and seq_pos is not None and start <= genomic_pos < stop:
                        seq = read.get_forward_sequence()
                        if read.is_reverse:
                            aligned_read_nucleotide = BASE_MATCH[seq[::-1][seq_pos]]
                        else:
                            aligned_read_nucleotide = seq[seq_pos]
                        cov_details[genomic_pos - start][BASE_TO_INDEX_MAPPING[aligned_read_nucleotide]] += 1
        cov_details[:, len(BASES)] = np.sum(cov_details, axis=1)
        return cov_details