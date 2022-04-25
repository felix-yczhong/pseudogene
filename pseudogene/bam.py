import numpy as np

from smaca.bam import Bam
from smaca.utils import get_total_depth

class PseudoBam(Bam):
    # 0-based coordinate
    def get_alignment_cov_util(self, pos_info):
        chr, pos = pos_info
        pos = pos - 1
        return self.samfile.count_coverage(chr, pos, pos + 1)

    def get_cov_ranges_(self, true_gene_ranges, pseudo_gene_ranges):
        # use np.array views to avoid copies
        cov = np.zeros(len(true_gene_ranges) + len(pseudo_gene_ranges))
        for i, (chr, start, end) in enumerate(true_gene_ranges):
            # x = self.samfile.count_coverage(chr, start, end).sum(axis=1) / 4
            _, cov[i] = get_total_depth(self.samfile, chr, start, end)
        for j, (chr, start, end) in enumerate(pseudo_gene_ranges):
            _, cov[len(true_gene_ranges) + j] = get_total_depth(self.samfile, chr, start, end)
        return cov