from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Tuple)
from src.data_types import (
    BGC_Module_Prediction,
    Chirality,
    LogProb,
    NRP_Monomer,
    BGC_Variant,
    NRP_Variant,
    NRP_Fragment,
)
from src.NewMatcher.dp_helper import DP_Helper
from src.NewMatcher.alignment_types import Alignment, alignment_score, Match
from src.NewMatcher.dp import get_alignment


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
                  key=lambda match: match.score, reverse=True)