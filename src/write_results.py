from typing import (
    Callable,
    Iterable,
    List,
    TypeVar,
    Tuple
)
from src.NewMatcher.alignment_types import Match
from src.data_types import BGC_Variant, NRP_Variant
from src.nerpa_pipeline.handle_rban import GraphRecord

from pathlib import Path
from io import StringIO
import csv
import yaml
from itertools import groupby

T = TypeVar('T')
U = TypeVar('U')

def sort_groupby(items: Iterable[T],
                 key: Callable[[T], U]) -> Iterable[Tuple[U, Iterable[T]]]:
    return groupby(sorted(items, key=key, reverse=True), key=key)


def write_yaml(data, out_file: Path,
               compress: bool = False):
    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    if not compress:
        yaml.Dumper.ignore_aliases = lambda *args: True

    with open(out_file, 'w') as out:
        yaml.dump(data, out,
                  default_flow_style=None, sort_keys=False)


def build_report(matches: List[Match]) -> str:
    result = StringIO()
    csv_writer = csv.DictWriter(result, fieldnames=('Score', 'NRP_ID', 'BGC_ID'), delimiter='\t')
    csv_writer.writeheader()
    csv_writer.writerows({'Score': match.normalized_score,
                          'NRP_ID': match.nrp_variant.nrp_id,
                          'BGC_ID': match.bgc_variant.bgc_id}
                         for match in matches)
    return result.getvalue()


def write_records_per_id(matches: List[T],
                         output_dir: Path,
                         get_id: Callable[[T], str]):
    for id_, id_matches in sort_groupby(matches, get_id):  # python sort is stable so groups will be sorted by score
        (output_dir / Path(id_)).write_text('\n\n'.join(map(str, matches)))


def write_results(bgc_variants: List[BGC_Variant],
                  nrp_variants: List[NRP_Variant],
                  matches: List[Match],
                  rban_graphs: List[GraphRecord],
                  output_dir: Path):
    (output_dir / Path('report.tsv')).write_text(build_report(matches))

    write_yaml(rban_graphs, output_dir / Path('rban_graphs.yaml'))

    (output_dir / Path('NRP_variants')).mkdir()
    for nrp_id, nrp_id_variants in sort_groupby(nrp_variants, key=lambda nrp_variant: nrp_variant.nrp_id):
        write_yaml(nrp_variants, output_dir / Path(f'NRP_variants/{nrp_id}.yaml'))

    (output_dir / Path('matches_details')).mkdir()
    write_yaml(matches, output_dir / Path('matches_details/matches.yaml'),
               compress=True)

    (output_dir / Path('matches_details/per_BGC')).mkdir()
    write_records_per_id(matches, output_dir / Path('matches_details/per_BGC'),
                         get_id=lambda match: match.bgc_variant.bgc_id)

    (output_dir / Path('matches_details/per_NRP')).mkdir()
    write_records_per_id(matches, output_dir / Path('matches_details/per_NRP'),
                         get_id=lambda match: match.nrp_variant.nrp_id)


