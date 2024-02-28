from typing import (
    Any,
    Dict,
    Iterable,
    List,
    NamedTuple,
    Tuple,
    TypeVar,
    Set,
    Union
)

from src.data_types import (
    BGC_Module_Prediction,
    BGC_Variant,
    Chirality,
    LogProb,
    NRP_Monomer,
    NRP_Variant,
    NRP_Fragment,
    MonomerResidue
)
from src.NewMatcher.dp import get_alignment
from src.NewMatcher.dp_helper import DP_Helper
from src.NewMatcher.dp_config import DP_Config, load_config
from src.NewMatcher.matcher import get_matches

from src.NewMatcher.alignment_types import AlignmentStep, Alignment, show_alignment
from itertools import groupby, pairwise, chain
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


T = TypeVar('T')
def split_iterable(xs: Iterable[T], sep: T) -> Iterable[Iterable[T]]:
    return (g
            for k, g in groupby(xs, key=lambda x: x != sep)
            if k)


TestAlignment = Tuple[List[Union[MonomerResidue, None]],
List[Union[MonomerResidue, None]]]

def get_iterative(xs: list) -> Set[int]:
    iterative = set()
    for i, (fst, snd) in enumerate(pairwise(xs)):
        if snd == '*':
            iterative.add(i)
    return iterative


def parse_contig(modules: List[str],
                 gene_id: str,
                 iterative_gene: bool) -> List[BGC_Module_Prediction]:
    iterative = get_iterative(modules)
    modules = filter(lambda x: x != '*', modules)
    return [BGC_Module_Prediction(residue_score=defaultdict(lambda: -10, {MonomerResidue(residue): 10}),
                                  modifications=(),
                                  gene_id=gene_id,
                                  module_idx=i,
                                  iterative_module=i in iterative,
                                  iterative_gene=iterative_gene)
            for i, residue in enumerate(modules)]


def parse_contigs(pred_str: str) -> List[BGC_Module_Prediction]:
    genes = eval(pred_str[len('BGC: '):])
    iterative = get_iterative(genes)
    genes = filter(lambda x: x != '*', genes)
    return list(chain(*(parse_contig(contig, str(i), i in iterative)
                        for i, contig in enumerate(genes))))


def parse_nrp_monomers(mons_str: str) -> List[NRP_Monomer]:
    return [NRP_Monomer(residue=MonomerResidue(residue),
                        modifications=(),
                        chirality=Chirality.UNKNOWN,
                        rban_name=residue,
                        rban_idx=i)
            for i, residue in enumerate(eval(mons_str[len('NRP: '):]))]


def parse_alignment(alignment_str: str) -> TestAlignment:
    def parse_monomer(monomer: str) -> Union[MonomerResidue, None]:
        return MonomerResidue(monomer) if monomer != '-' else None
    bgc, nrp = eval(alignment_str[len('Alignment: '):])
    return list(map(parse_monomer, bgc)), list(map(parse_monomer, nrp))


@dataclass
class AlignmentTest:
    predictions: List[BGC_Module_Prediction]
    nrp_monomers: List[NRP_Monomer]
    correct_alignment: TestAlignment
    num_mods: int
    test_description: str

    def __init__(self, test_str_lines):
        lines = iter(test_str_lines)
        self.test_description = next(lines)
        self.num_mods = int(next(lines).split()[-1])
        self.predictions = parse_contigs(next(lines))
        self.nrp_monomers = parse_nrp_monomers(next(lines))
        self.correct_alignment = parse_alignment(next(lines))


def load_tests() -> Iterable[AlignmentTest]:
    with open(Path(__file__).parent / Path("tests_prototype.txt"), 'r') as inp:
        tests_str = split_iterable(filter(lambda line: not line.startswith('###'),
                                          inp.readlines()),
                                   '\n')

    return map(AlignmentTest, tests_str)


def parse_result_alignment(alignment: Alignment) -> TestAlignment:
    def best_pred(bgc_module: BGC_Module_Prediction):
        return max(bgc_module.residue_score.keys(),
                   key=lambda res: bgc_module.residue_score[res])

    bgc = [best_pred(trans.bgc_module) if trans.bgc_module else None
           for trans in alignment
           if trans.action not in ('iterate_module', 'iterate_gene')]
    nrp = [trans.nrp_monomer.residue if trans.nrp_monomer else None
           for trans in alignment
           if trans.action not in ('iterate_module', 'iterate_gene')]
    return bgc, nrp


def check(alignment: Alignment, correct_alignment: TestAlignment) -> bool:
    return parse_result_alignment(alignment) == correct_alignment


def run_tests(tests: Iterable[AlignmentTest]) -> bool:
    dp_helper = DP_Helper(load_config(Path(__file__).parent / 'dp_config.yaml'))
    for i, test in enumerate(tests):
        if test.num_mods > 0:
            continue
        alignment = get_alignment(assembly_line=test.predictions,
                                  nrp_monomers=test.nrp_monomers,
                                  dp_helper=dp_helper)

        if not check(alignment, test.correct_alignment):
            print(f'Test {test.test_description} failed!')
            print('Alignment:')
            print(show_alignment(alignment))
            print('Correct alignment:')
            print(', '.join(res if res else '-' for res in test.correct_alignment[0]))
            print(', '.join(res if res else '-' for res in test.correct_alignment[1]))
            return False
        else:
            print(f'Test {test.test_description} passed')

    return True


def test_matcher():
    nrp_fragments = [NRP_Fragment(monomers=[NRP_Monomer(residue='val',
                                                        modifications=(),#(NRP_Monomer_Modification.METHYLATION,),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Val',
                                                        rban_idx=0),
                                            NRP_Monomer(residue='leu',
                                                        modifications=(),#(NRP_Monomer_Modification.UNKNOWN,),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Leu',
                                                        rban_idx=1)],
                                  is_cyclic=False),
                     NRP_Fragment(monomers=[NRP_Monomer(residue='leu',
                                                        modifications=(),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Leu',
                                                        rban_idx=0),
                                            NRP_Monomer(residue='leu',
                                                        modifications=(),#(NRP_Monomer_Modification.UNKNOWN, NRP_Monomer_Modification.METHYLATION),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Leu',
                                                        rban_idx=1),
                                            NRP_Monomer(residue='val',
                                                        modifications=(),#(NRP_Monomer_Modification.METHYLATION,),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Val',
                                                        rban_idx=2)],
                                  is_cyclic=True)]

    nrp_variant = NRP_Variant(fragments=nrp_fragments, nrp_id='dragon_sneeze')

    bgc_preds = [BGC_Module_Prediction(residue_score=defaultdict(lambda: -3, {'val': -1}),
                                       modifications=(),
                                       iterative_module=False,
                                       iterative_gene=False,
                                       gene_id='gene1',
                                       module_idx=0),
                 BGC_Module_Prediction(residue_score=defaultdict(lambda: -3, {'leu': -1}),
                                       modifications=(), #(BGC_Module_Modification.EPIMERIZATION, BGC_Module_Modification.METHYLATION,),
                                       iterative_module=True,
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


def main():
    run_tests(load_tests())


if __name__ == '__main__':
    main()
