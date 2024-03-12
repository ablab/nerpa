from __future__ import annotations
from typing import (
    Dict,
    NamedTuple,
    List,
    Tuple,
    Union
)
from src.data_types import (
    BGC_Module,
    NRP_Monomer,
    BGC_Variant,
    LogProb,
    NRP_Variant,
    NRP_Fragment
)
from dataclasses import dataclass
from prettytable import PrettyTable
from io import StringIO


class AlignmentStep(NamedTuple):
    bgc_module: Union[BGC_Module, None]
    nrp_monomer: Union[NRP_Monomer, None]
    score: Union[LogProb, Tuple[LogProb, LogProb, LogProb]]
    action: str

    def get_score(self) -> LogProb:
        return sum(self.score) if self.action == 'match' else self.score

    def to_dict(self) -> Dict[str, str]:
        NA = '---'
        return {'Gene': self.bgc_module.gene_id if self.bgc_module else NA,
                'A-domain_idx': self.bgc_module.module_idx if self.bgc_module else NA,
                'Top_scoring_residue': max(self.bgc_module.residue_score.items(),
                                            key=lambda p: p[1])[0] if self.bgc_module else NA,
                'Modifying_domains': ','.join(mod.name for mod in self.bgc_module.modifications)
            if self.bgc_module and self.bgc_module.modifications else NA,
                'NRP_residue': self.nrp_monomer.residue if self.nrp_monomer else NA,
                'NRP_chirality': self.nrp_monomer.chirality.name if self.nrp_monomer else NA,
                'NRP_modifications': ','.join(mod.name for mod in self.nrp_monomer.modifications)
            if self.nrp_monomer and self.nrp_monomer.modifications else NA,
                'rBAN_name': self.nrp_monomer.rban_name if self.nrp_monomer else NA,
                'rBAN_idx': self.nrp_monomer.rban_idx if self.nrp_monomer else NA,
                'Alignment_step': self.action,
                'Score': round(self.get_score(), 3),
                'ResidueScore': round(self.score[0], 3) if self.action == 'match' else round(self.score, 3),
                'MethylationScore': round(self.score[1], 3) if self.action == 'match' else NA,
                'ChiralityScore': round(self.score[2], 3) if self.action == 'match' else NA}

    @classmethod
    def from_yaml_dict(cls, data: dict) -> AlignmentStep:
        return cls(bgc_module=BGC_Module.from_yaml_dict(data['bgc_module']) if data['bgc_module'] else None,
                   nrp_monomer=NRP_Monomer.from_yaml_dict(data['nrp_monomer']) if data['nrp_monomer'] else None,
                   score=data['score'],
                   action=data['action'])


Alignment = List[AlignmentStep]

def alignment_score(alignment: Alignment) -> LogProb:
    return sum(alignment_step.get_score() for alignment_step in alignment)


def show_alignment(alignment: Alignment) -> str:
    rows = [alignment_step.to_dict()
            for alignment_step in alignment]

    t = PrettyTable(rows[0].keys(), align='l', border=False)
    t.add_rows(row.values() for row in rows)
    return str(t)

@dataclass
class Match:
    bgc_variant: BGC_Variant
    nrp_variant: NRP_Variant
    alignments: List[Alignment]  # alignments of each fragment
    normalized_score: float

    def raw_score(self) -> LogProb:
        return sum(map(alignment_score, self.alignments))

    def __str__(self):
        out = StringIO()
        out.write('\n'.join([f'Genome={self.bgc_variant.genome_id}',
                             f'BGC={self.bgc_variant.bgc_id}',
                             f'BGC_variant={self.bgc_variant.variant_idx}',
                             f'NRP={self.nrp_variant.nrp_id}',
                             f'NRP_variant={self.nrp_variant.variant_idx}',
                             f'NormalisedScore={self.normalized_score}',
                             f'Score={self.raw_score()}']))
        out.write('\n')

        for i, alignment in enumerate(self.alignments):
            if len(self.alignments) > 1:
                out.write(f'Fragment_#{i}\n')
            out.write(show_alignment(alignment) + '\n')

        return out.getvalue()

    def to_dict_light(self) -> dict:  # because full Match is too big to write
        return {'Genome': self.bgc_variant.genome_id,
                'BGC': self.bgc_variant.bgc_id,
                'BGC_variant_idx': self.bgc_variant.variant_idx,
                'NRP': self.nrp_variant.nrp_id,
                'NRP_variant_idx': self.nrp_variant.variant_idx,
                'NormalisedScore': self.normalized_score,
                'Score': self.raw_score(),
                'Alignments': [[alignment_step.to_dict()
                                for alignment_step in alignment]
                               for alignment in self.alignments]}

    @classmethod
    def from_yaml_dict(cls, data: dict) -> Match:
        return cls(bgc_variant=BGC_Variant.from_yaml_dict(data['bgc_variant']),
                   nrp_variant=NRP_Variant.from_yaml_dict(data['nrp_variant']),
                   alignments=[[AlignmentStep.from_yaml_dict(alignment_step_data)
                                for alignment_step_data in alignment_data]
                               for alignment_data in data['alignments']],
                   normalized_score=data['normalised_score'])