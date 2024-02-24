from typing import List, Literal
from src.NewMatcher.alignment_types import Match
from pathlib import Path
from io import StringIO
import csv
import yaml
from itertools import groupby
from dataclasses import asdict


def write_yaml(data, out_file: Path):
    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    yaml.Dumper.ignore_aliases = lambda *args: True
    with open(out_file, 'w') as out:
        yaml.dump(data, out,
                  default_flow_style=None, sort_keys=False)


def build_report(matches: List[Match]) -> str:
    result = StringIO()
    csv_writer = csv.DictWriter(result, fieldnames=('Score', 'NRP_ID', 'BGC_ID'), delimiter='\t')
    csv_writer.writeheader()
    csv_writer.writerows({'Score': match.normalised_score,
                          'NRP_ID': match.nrp_variant.nrp_id,
                          'BGC_ID': match.bgc_variant.bgc_id}
                         for match in matches)
    return result.getvalue()


def write_matches_per_id(matches: List[Match], output_dir: Path,
                         id_=Literal['BGC', 'NRP']):
    get_id = (lambda m: m.bgc_variant.bgc_id) if id_ == 'BGC' else \
        (lambda m: m.nrp_variant.nrp_id)
    matches_per_id = groupby(sorted(matches, key=get_id, reverse=True),
                              key=get_id)  # python sort is stable so groups will be sorted by score

    for id_, id_matches in matches_per_id:
        (output_dir / Path(id_)).write_text('\n\n'.join(map(str, matches)))


def write_results(matches: List[Match], output_dir: Path):
    (output_dir / Path('report.tsv')).write_text(build_report(matches))

    (output_dir / Path('details')).mkdir()
    write_yaml(matches, output_dir / Path('matches.yaml'))

    (output_dir / Path('details/per_BGC')).mkdir()
    write_matches_per_id(matches, output_dir / Path('details/per_BGC'), id_=Literal['BGC'])

    (output_dir / Path('details/per_NRP')).mkdir()
    write_matches_per_id(matches, output_dir / Path('details/per_BGC'), id_=Literal['NRP'])
