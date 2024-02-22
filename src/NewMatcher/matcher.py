from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Tuple)
from src.data_types import (
    BGC_Module_Modification,
    BGC_Module_Prediction,
    Chirality,
    LogProb,
    NRP_Monomer,
    NRP_Monomer_Modification,
    BGC_Variant,
    NRP_Variant,
    NRP_Fragment,
)
from src.NewMatcher.dp_helper import DP_Helper
from src.NewMatcher.dp_config import DP_Config, load_config
from src.NewMatcher.alignment_types import Alignment, alignment_score, Match
from src.NewMatcher.dp import get_alignment

from collections import defaultdict
from pathlib import Path


def get_alignment_fragment(assembly_line: List[BGC_Module_Prediction],
                           nrp_fragment: NRP_Fragment,
                           dp_helper: DP_Helper) -> Alignment:

    nrp_cyclic_shifts = [nrp_fragment.monomers] if not nrp_fragment.is_cyclic \
            else [nrp_fragment.monomers[i:] + nrp_fragment.monomers[:i]
                  for i in range(len(nrp_fragment.monomers))]

    alignments = [get_alignment(assembly_line, nrp_cyclic_shift, dp_helper)
                  for nrp_cyclic_shift in nrp_cyclic_shifts]
    return max(alignments, key=alignment_score)


def null_hypothesis_score(nrp_fragment: NRP_Fragment,
                          dp_helper: DP_Helper) -> LogProb:
    return sum(dp_helper.null_hypothesis_score(mon)
               for mon in nrp_fragment.monomers)


def get_match(bgc_variant: BGC_Variant,
              nrp_variant: NRP_Variant,
              dp_helper: DP_Helper) -> Match:
    fragment_alignments = [get_alignment_fragment(bgc_variant.tentative_assembly_line,
                                                  bgc_fragment,
                                                  dp_helper)
                           for bgc_fragment in nrp_variant.fragments]
    final_score = sum(alignment_score(alignment) - null_hypothesis_score(fragment, dp_helper)
                      for fragment, alignment in zip(nrp_variant.fragments, fragment_alignments))
    return Match(bgc_variant,
                 nrp_variant,
                 fragment_alignments,
                 final_score)


def get_matches(bgc_variants: List[BGC_Variant],
                nrp_variants: List[NRP_Variant],
                dp_config: DP_Config) -> List[Match]:
    dp_helper = DP_Helper(dp_config)
    return sorted([get_match(bgc_variant, nrp_variant, dp_helper)
                   for bgc_variant in bgc_variants
                   for nrp_variant in nrp_variants],
                  key=lambda match: match.normalised_score, reverse=True)

def main():
    nrp_fragments = [NRP_Fragment(monomers=[NRP_Monomer(residue='val',
                                                        modifications=(),#(NRP_Monomer_Modification.METHYLATION,),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Val'),
                                            NRP_Monomer(residue='leu',
                                                        modifications=(),#(NRP_Monomer_Modification.UNKNOWN,),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Leu')],
                                  is_cyclic=False,
                                  rban_indexes=[0,1]),
                     NRP_Fragment(monomers=[NRP_Monomer(residue='leu',
                                                        modifications=(),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Leu'),
                                            NRP_Monomer(residue='leu',
                                                        modifications=(),#(NRP_Monomer_Modification.UNKNOWN, NRP_Monomer_Modification.METHYLATION),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Leu'),
                                            NRP_Monomer(residue='val',
                                                        modifications=(),#(NRP_Monomer_Modification.METHYLATION,),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Val')],
                                  is_cyclic=True,
                                  rban_indexes=[3, 7, 2])]

    nrp_variant = NRP_Variant(fragments=nrp_fragments, nrp_id='dragon_sneeze')

    bgc_preds = [BGC_Module_Prediction(residue_score=defaultdict(lambda: -3, {'val': -1}),
                                       modifications=(),
                                       iterative_module=False,
                                       iterative_gene=False,
                                       gene_id='gene1',
                                       module_idx=0),
                 BGC_Module_Prediction(residue_score=defaultdict(lambda: -3, {'leu': -1}),
                                       modifications=(), #(BGC_Module_Modification.EPIMERIZATION, BGC_Module_Modification.METHYLATION,),
                                       iterative_module=False,
                                       iterative_gene=False,
                                       gene_id='gene1',
                                       module_idx=1)]
    bgc_variant = BGC_Variant(tentative_assembly_line=bgc_preds,
                              genome_id='genome#189',
                              bgc_id='bgc#39')


    dp_config = load_config(Path(__file__).parent / 'dp_config.yaml')
    with (Path(__file__).parent / Path('test_matches_output.txt')).open('w') as out:
        for match in get_matches([bgc_variant], [nrp_variant], dp_config):
            out.write(str(match) + '\n\n')

if __name__ == "__main__":
    main()