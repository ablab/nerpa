import handle_helper


def get_PK_hybrids(dirname):
    PK_hybrids_list = []
    domains = handle_helper.get_domains_list(dirname)
    orf_ori = handle_helper.get_orf_orientation(dirname)

    # other options: 'PKS_AT'
    PK_DOMAIN = 'PKS_KS'
    A_DOMAIN = 'AMP-binding'

    for orfds in domains:
        step = -1 if orfds and orf_ori[orfds[0][0]] == '-' else 1
        coding_domains_gen = filter(lambda x: x[2] == A_DOMAIN or x[2] == PK_DOMAIN, orfds[::step])
        prev_d = None
        for d in coding_domains_gen:
            if prev_d and prev_d[2] == A_DOMAIN and d[2] == PK_DOMAIN:
                PK_hybrids_list.append(prev_d[0] + "_" + prev_d[1].split('_')[-1])
            prev_d = d

    return PK_hybrids_list


def test(dirname):
    orf_ori = handle_helper.get_orf_orientation(dirname)
    for orf in handle_helper.get_domains_list(dirname):
        print(orf[0][0])
        step = -1 if orf_ori[orf[0][0]] == '-' else 1
        for d in orf[::step]:
            print('\t', d[1:])
    print(get_PK_hybrids(dirname))
    print()


if __name__ == '__main__':
    dirname = '/Data/Projects/CAB/Nerpa/nerpa/resources/test_dataset/antismash/BGC0000383'
    test(dirname)
    # NOT SUPPORTED: PK module on a different orf: BGC0000963, BGC0001058
