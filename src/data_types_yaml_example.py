from src.data_types import *


def main():
    nrp_fragments = [NRP_Fragment(monomers=[NRP_Monomer(residue='arg',
                                                        modifications=(NRP_Monomer_Modification.METHYLATION,),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Arg',
                                                        rban_idx=0),
                                            NRP_Monomer(residue='leu',
                                                        modifications=(NRP_Monomer_Modification.UNKNOWN,),
                                                        chirality=Chirality.L,
                                                        rban_name='Leu',
                                                        rban_idx=1)],
                                  is_cyclic=False),
                     NRP_Fragment(monomers=[NRP_Monomer(residue='ile',
                                                        modifications=(),
                                                        chirality=Chirality.UNKNOWN,
                                                        rban_name='Ile',
                                                        rban_idx=0),
                                            NRP_Monomer(residue='pro',
                                                        modifications=(NRP_Monomer_Modification.UNKNOWN, NRP_Monomer_Modification.METHYLATION),
                                                        chirality=Chirality.L,
                                                        rban_name='Pro',
                                                        rban_idx=1),
                                            NRP_Monomer(residue='val',
                                                        modifications=(NRP_Monomer_Modification.METHYLATION,),
                                                        chirality=Chirality.D,
                                                        rban_name='Val',
                                                        rban_idx=2)],
                                  is_cyclic=True)]
    nrp_variants = [NRP_Variant(fragments=[nrp_fragments[0]], nrp_id='terrific_nrp#42'),
                    NRP_Variant(fragments=[nrp_fragments[1]], nrp_id='terrific_nrp#48')]
    bgc_preds = [BGC_Module(residue_score={'arg': 0.8, 'leu': 0.7},
                            modifications=(),
                            iterative_module=False,
                            iterative_gene=False,
                            gene_id='gene1',
                            module_idx=0),
                 BGC_Module(residue_score={'trp': 0.8, 'phe': 0.7},
                            modifications=(BGC_Module_Modification.EPIMERIZATION, BGC_Module_Modification.METHYLATION,),
                            iterative_module=True,
                            iterative_gene=False,
                            gene_id='gene1',
                            module_idx=1)]
    bgc_variant = BGC_Variant(tentative_assembly_line=bgc_preds,
                              genome_id='genome#189',
                              bgc_id='bgc#39')

    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    yaml.Dumper.ignore_aliases = lambda *args: True
    with open('test_bgc_predictions_output.yaml', 'w') as out:
        yaml.dump(asdict(bgc_variant), out,
                  default_flow_style=None, sort_keys=False)
    with open('test_nrp_variants_output.yaml', 'w') as out:
        yaml.dump(list(map(asdict, nrp_variants)), out,
                  default_flow_style=None, sort_keys=False)


if __name__ == "__main__":
    main()

