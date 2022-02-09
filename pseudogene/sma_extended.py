import tempfile

import numpy as np
import pandas as pd
import click
from joblib import Parallel, delayed

from smaca.sma import SmaCalculator
from smaca.bam import Bam
import smaca.constants as C
from pseudogene.utils import *

class SmaCalculatorExtended(SmaCalculator):
    def read_config(self, config, ref, gene):
        self.ref = ref
        self.gene = gene

        self.vcf_path_template = config["temp_files"]["vcf_path"]
        self.nirvana_output_path_template = config["temp_files"]["nirvana_path"]
        self.summary_template = config["temp_files"]["summary"]
        self.details_template = config["temp_files"]["details"]

        self.fasta_loc = config["tools"]["fasta_loc"][ref]
        self.nirvana_path = config["tools"]["nirvana_path"]

        self.true_gene = config[gene]["true_gene"]
        self.pseudo_genes = config[gene]["pseudo_genes"]

        seq_align_map = pd.read_csv(get_project_root() / config["seq_align_map"][ref][gene].format(gene=gene), sep=',', comment='#', header=[0, 1, 2])
        self.seq_align_map = seq_align_map.sort_values(by=("True Gene", self.true_gene, "Start"), axis=0, ascending=True, ignore_index=True)
        force = pd.read_csv(get_project_root() / config["force"][self.ref][gene].format(gene=gene), sep=',', comment='#', header=[0, 1, 2])
        self.force = force.sort_values(by=("True Gene", self.true_gene, "Start"), axis=0, ascending=True, ignore_index=True)

        section = config[gene][ref]
        self.true_gene_region = section["true_gene_region"]
        self.pseudo_gene_regions = section["pseudo_gene_regions"]
        self.NM_number = section["NM_number"]
        self.CNV_ratio = section["CNV_ratio_range"]

        true_gene_region = convert_coordinate_one_to_zero(self.true_gene_region)
        self.true_gene_region = {f"{self.true_gene}_{pos}": genomic_range for pos, genomic_range in enumerate([true_gene_region])}
        pseudo_gene_region = [convert_coordinate_one_to_zero(pseudo) for pseudo in self.pseudo_gene_regions]
        self.pseudo_gene_regions = {f"{self.pseudo_genes[0]}_{pos}": genomic_range for pos, genomic_range in enumerate([pseudo_gene_region[0]])}
        self.gene_regions = {**self.true_gene_region, **self.pseudo_gene_regions}

    def read_align_map(self):
        self.num_seq = len(self.seq_align_map)
        self.num_true_gene = len(self.seq_align_map["True Gene"].groupby(axis=1, level=0))
        self.num_pseudo_gene = len(self.seq_align_map["Pseudo Gene"].groupby(axis=1, level=0))
        self.num_gene = self.num_true_gene + self.num_pseudo_gene

        true_genomic_ranges = self.seq_align_map['True Gene'][self.true_gene][['Chrom', 'Start', 'End']].values.tolist()
        self.true_gene_pos = {f"{self.true_gene}_{pos}": genomic_range for pos, genomic_range in enumerate(true_genomic_ranges)}
        pseudo_genomic_ranges = self.seq_align_map['Pseudo Gene'][self.pseudo_genes[0]][['Chrom', 'Start', 'End']].values.tolist()
        self.pseudo_gene_pos = {f"{self.pseudo_genes[0]}_{pos}": genomic_range for pos, genomic_range in enumerate(pseudo_genomic_ranges)}

    def __init__(self, bam_list, ref, config, gene, n_jobs=1):
        self.read_config(config, ref, gene)
        self.read_align_map()

        self.bam_list = np.array(bam_list)
        self.n_bam = len(self.bam_list)
        # number of reads that align to SMN1 at position x
        self.D1_ij = np.zeros((self.n_bam, self.num_seq))
        # number of reads that align to SMN2 at position x
        self.D2_ij = np.zeros((self.n_bam, self.num_seq))
        # total number of reads aligned to the SMN1 region at position j
        # and the analogous SMN2 region
        self.r_ij = np.zeros((self.n_bam, self.num_seq))
        # average coverage for HK gene k
        self.H_ik = np.zeros((self.n_bam, len(C.POSITIONS[ref]["GENES"])))
        # average coverage for the SMNx gene region
        self.c_ix = np.zeros((self.n_bam, self.num_gene))
        # scaled coverage for SMN1 and SMN2 for HK gene k
        self.z_ik = np.zeros((self.n_bam, len(C.POSITIONS[ref]["GENES"])))
        # scaled proportion of SMN reads that align to SMN1 in exon 7
        self.pi_ij = np.zeros((self.n_bam, self.num_seq))
        # averaged scaled coverages for SMN1 and SMN2 in the N subjects
        self.zmean_k = np.zeros((len(C.POSITIONS[ref]["GENES"]), ))
        # weighted average of the coverage of SMN1 to our K housekeeping genes
        self.theta_i = np.zeros((self.n_bam, ))
        # std of coverage in housekeeping genes in sample i
        self.std_i = np.zeros((self.n_bam, ))
        # std of coverage in housekeeping gene k for each sample
        self.std_k = np.zeros((len(C.POSITIONS[ref]["GENES"]), ))
        # get consensus sequence on dup. markers
        self.dup_id = np.empty((self.n_bam, len(C.POSITIONS[ref]["DUP_MARK"])),
                                dtype="S100")

        bam_list = np.array(bam_list)
        n_bam = len(bam_list)

        if n_jobs == 1:
            with click.progressbar(length=self.n_bam,
                                    label='processing BAM files') as bar:
                for i in bar:
                    self.compute(bam_list[i], i, self.D1_ij, self.D2_ij,
                                    self.H_ik, self.c_ix, self.dup_id, ref, self)
        else:
            from pathlib import Path

            tmp_dir = tempfile.mkdtemp()

            D1_ij_fname_memmap = Path(tmp_dir).joinpath("D1_ij_memmap")
            D1_ij_memmap = np.memmap(D1_ij_fname_memmap.as_posix(),
                                        dtype=np.float,
                                        shape=(n_bam,
                                            self.num_seq),
                                        mode='w+')

            D2_ij_fname_memmap = Path(tmp_dir).joinpath("D2_ij_memmap")
            D2_ij_memmap = np.memmap(D2_ij_fname_memmap,
                                        dtype=np.float,
                                        shape=(n_bam,
                                            self.num_seq),
                                        mode='w+')

            H_ik_fname_memmap = Path(tmp_dir).joinpath("H_ik_memmap")
            H_ik_memmap = np.memmap(H_ik_fname_memmap,
                                    dtype=np.float,
                                    shape=(n_bam,
                                            len(C.POSITIONS[ref]["GENES"])),
                                    mode='w+')

            c_ix_fname_memmap = Path(tmp_dir).joinpath("c_ix_memmap")
            c_ix_memmap = np.memmap(c_ix_fname_memmap,
                                    dtype=np.float,
                                    shape=(n_bam,
                                            self.num_gene),
                                    mode='w+')

            dup_id_fname_memmap = Path(tmp_dir).joinpath("dup_id_memmap")
            dup_id_memmap = np.memmap(
                dup_id_fname_memmap,
                dtype="S100",
                shape=(n_bam, len(C.POSITIONS[ref]["DUP_MARK"])),
                mode='w+')

            Parallel(n_jobs=n_jobs)(delayed(
                self.compute)(bam_list[idx], idx, D1_ij_memmap, D2_ij_memmap,
                                H_ik_memmap, c_ix_memmap, dup_id_memmap, ref, self)
                                    for idx in range(self.n_bam))

            self.D1_ij[:] = D1_ij_memmap[:]
            self.D2_ij[:] = D2_ij_memmap[:]
            self.H_ik[:] = H_ik_memmap[:]
            self.c_ix[:] = c_ix_memmap[:]
            self.dup_id[:] = dup_id_memmap[:]

        self.r_ij = self.D1_ij + self.D2_ij
        #TODO: consider using only "the bests" HK
        self.z_ik = self.c_ix.sum(axis=1).reshape((self.n_bam, 1)) / self.H_ik
        self.std_k = np.std(self.H_ik, axis=0)
        self.std_i = np.std(self.H_ik, axis=1)
        self.zmean_k = self.z_ik.sum(axis=0) / self.n_bam
        self.theta_i = (self.z_ik / self.zmean_k).sum(axis=1) / len(
            C.POSITIONS[ref]["GENES"])
        self.pi_ij = self.theta_i.reshape(
            (self.n_bam, 1)) * (self.D1_ij / self.r_ij)

    @staticmethod
    def compute(bam_file, i, D1_ij, D2_ij, H_ik, c_ix, dup_id, ref, self):
        b = Bam(bam_file)
        D1_ij[i] = b.get_cov_ranges(self.true_gene_pos)
        D2_ij[i] = b.get_cov_ranges(self.pseudo_gene_pos)
        H_ik[i] = b.get_cov_ranges(C.POSITIONS[ref]["GENES"])
        c_ix[i] = b.get_cov_ranges(self.gene_regions)
        for d in range(len(C.POSITIONS[ref]["DUP_MARK"])):
            dup_id[i][d] = b.get_consensus_sequence(
                list(C.POSITIONS[ref]["DUP_MARK"].values())[d])