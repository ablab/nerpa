from typing import Dict, NamedTuple, List
from src.data_types import (
    BGC_Module_Prediction,
    NRP_Monomer,
    BGC_Variant,
    LogProb,
    NRP_Variant,
    NRP_Fragment)
from dataclasses import dataclass
from io import StringIO
import csv


class AlignmentStep(NamedTuple):
    bgc_module_pos: int  # rank is a dubious name but I don't want to confuse them with actual indexes
    nrp_monomer_pos: int
    score: LogProb
    action: str

    def to_dict(self,
                bgc_predictions: List[BGC_Module_Prediction],
                nrp_fragment: NRP_Fragment) -> Dict[str, str]:  # fragment is needed only to retrieve rban indexes
        NA = '---'
        bgc_module = bgc_predictions[self.bgc_module_pos] if self.bgc_module_pos else None
        nrp_monomer = nrp_fragment.monomers[self.nrp_monomer_pos] if self.nrp_monomer_pos else None
        return {'Gene': bgc_module.gene_id if bgc_module else NA,
                'A-domain_idx': bgc_module.module_idx if bgc_module else NA,
                'Modifications': ','.join(mod.name for mod in bgc_module.modifications)
            if bgc_module and bgc_module.modifications else NA,
                'NRP_residue': nrp_monomer.residue if nrp_monomer else NA,
                'NRP_chirality': nrp_monomer.chirality.name if nrp_monomer else NA,
                'NRP_modifications': ','.join(ptm.name for ptm in nrp_monomer.ptms)
            if nrp_monomer and nrp_monomer.ptms else NA,
                'rBAN_name': nrp_monomer.rban_name if nrp_monomer else NA,
                'rBAN_idx': nrp_fragment.rban_indexes[self.nrp_monomer_pos] if nrp_monomer else NA,
                'Alignment_step': self.action,
                'Score': self.score}


Alignment = List[AlignmentStep]

def alignment_score(alignment: Alignment) -> LogProb:
    return sum(alignment_step.score for alignment_step in alignment)


@dataclass
class Match:
    bgc_variant: BGC_Variant
    nrp_variant: NRP_Variant
    alignments: List[Alignment]  # alignments of each fragment
    normalised_score: float

    def raw_score(self) -> LogProb:
        return sum(map(alignment_score, self.alignments))

    def __str__(self):
        out = StringIO('')
        out.write('\n'.join([f'Genome={self.bgc_variant.genome_id}',
                             f'BGC={self.bgc_variant.bgc_id}',
                             f'NRP={self.nrp_variant.nrp_id}',
                             f'NormalisedScore={self.normalised_score}',
                             f'Score={self.raw_score()}']))
        out.write('\n')

        for i, alignment in enumerate(self.alignments):
            if len(self.alignments) > 1:
                out.write(f'Fragment_#{i}\n')
            rows = [alignment_step.to_dict(self.bgc_variant.tentative_assembly_line,
                                           self.nrp_variant.fragments[i])
                    for alignment_step in alignment]
            csv_writer = csv.DictWriter(out, fieldnames=rows[0].keys(), delimiter='\t')
            csv_writer.writeheader()
            csv_writer.writerows(rows)