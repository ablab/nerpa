#include <cassert>
#include <iostream>
#include "Aminoacids.h"

namespace aminoacid {
    aminoacid::Aminoacids::Aminoacid aminoacid::Aminoacids::get_aminoacid(std::string aminoacid_name) {
        for (int i = 0; i < AMINOACID_CNT; ++i) {
            if (AMINOACID_NAMES[i] == aminoacid_name) {
                return static_cast<Aminoacid>(i);
            }
        }

        std::cerr << aminoacid_name << "\n";
        assert(false);
    }

    aminoacid::Aminoacids::Aminoacid aminoacid::Aminoacids::get_aminoacid_from_formula(std::string fotmula) {
        for (int i = 0; i < AMINOACID_CNT; ++i) {
            if (FORMULS[i] == fotmula) {
                return static_cast<Aminoacid>(i);
            }
        }

        assert(false);
    }

    bool aminoacid::Aminoacids::same(aminoacid::Aminoacids::Aminoacid a, aminoacid::Aminoacids::Aminoacid b) {
        if (FORMULS[a] == FORMULS[b]) return true;
        return false;
    }

    //uda=ureidoalanine
    //dhp=dehydroalanine
    //phg=phenylglycine
    //asn=asparagine
    //dab=diaminobutyric acid
    //hpg=hydroxyphenylglycine
    //me-asp=methylaspartic acid
    //asp=aspartic acid
    //kyn=kynurenine
    //me-glu=methylglutamic acid
    //orn=ornithine
    //sar=sarcosine
    //pip=pipecolinic acid
    //dpg=dihydroxyphenylglycine

    //cha=β-cyclohexyl-l-alanine
    //nle=l-norleucine
    //nva=l-norvaline
    //abu=l-2-aminobutyric acid

    //dhpg=3,5-hydroxyL-phenylglycine
    //bmt= (4R)-4-[(E)-2-butenyl]-4-methyl-L-threonine
    //dhb= 2,-3-dihydroxybenzolglycine synthetase

    //aad= 2-amino-adipic acid
    //dht= dehydro-threonine
    //vol=valinol
    //bht=beta-hydroxy-tyrosine
    //sal=salicylic acid
    //tcl=(4 S )-5,5,5-trichloro-leucine
    //aeo=2-amino-9,10-epoxy-8-oxodecanoic acid
    //ala_d= D-alanine

    //hse=homoserine
    //hty=homotyrosine biosynthetic cassette
    //apa=aminobutyric acid
    //hyv=a-hydroxy-isovaleric acid

    //cit=citrulline
    //dht=dehydro-threonine
    //ahp=3-amino-6-hydroxy-2-piperidone

    //b_ala= beta-alanine
    //allothr=D-Allothreonine

    const std::string Aminoacids::AMINOACID_NAMES[Aminoacids::AMINOACID_CNT] = {"trp", "ser", "gly", "uda", "thr",
                                                                                "dhp", "gln", "dab", "arg", "lys",
                                                                                "ala_d", "phe", "val", "cha", "dhpg",
                                                                                "phg", "his", "aeo", "bmt", "hse",
                                                                                "met", "ala", "tcl", "sal", "allothr",
                                                                                "b-ala", "dhb", "ile", "end", "leu",
                                                                                "gua", "hty", "glu", "bht", "hpg",
                                                                                "apa", "pro", "tyr", "hyv", "asn",
                                                                                "cit", "vol", "cys", "asp", "dht",
                                                                                "ahp", "orn", "apc", "abu", "aad",
                                                                                "pip", "dpg"};


    const std::string Aminoacids::FORMULS[Aminoacids::AMINOACID_CNT] = {"C11H10N2O", "C3H5NO2", "C2H3NO", "C4H7N3O2", "C4H7NO2",
                                                                        "C3H3NO", "C5H8N2O2", "C4H8N2O", "C6H12N4O", "C6H12N2O",
                                                                        "C3H5NO", "C9H9NO", "C5H9NO", "C9N15NO", "C8H7NO2",
                                                                        "C8H7NO", "C6H7N3O", "C10H15NO3", "C9H15NO2", "C4H7NO2",
                                                                        "C5H9NOS", "C3H5NO", "C6H8Cl3NO", "C7H4O2", "C4H7NO2",
                                                                        "C3H5NO", "--", "C6H11NO", "-", "C6H11NO",
                                                                        "-", "--", "C5H7NO3", "C9H9NO3", "C8H7NO2",
                                                                        "C4H7NO", "C5H7NO", "C9H9NO2", "--", "C4H6N2O2",
                                                                        "C6H11N3O2", "C5N11N", "C3H5NOS", "C4H5NO3", "--",
                                                                        "C5H8N2O", "C5H10N2O", "-", "C4H7NO", "C6H9NO3",
                                                                        "C6H9NO", "C8H7NO3"};
}