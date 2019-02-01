#include <cassert>
#include <iostream>
#include "Aminoacid.h"

namespace aminoacid {
    aminoacid::Aminoacid::AminoacidId aminoacid::Aminoacid::get_aminoacid(std::string aminoacid_name) {
        for (int i = 0; i < AMINOACID_CNT; ++i) {
            if (AMINOACID_NAMES[i] == aminoacid_name) {
                return static_cast<AminoacidId>(i);
            }
        }

        std::cerr << aminoacid_name << "\n";
        assert(false);
    }

    aminoacid::Aminoacid::AminoacidId aminoacid::Aminoacid::get_aminoacid_from_formula(std::string formula) {
        for (int i = 0; i < AMINOACID_CNT; ++i) {
            if (FORMULS[i] == formula) {
                return static_cast<AminoacidId>(i);
            }
        }

        //std::cerr << formula << "\n";
        //assert(false);

        return static_cast<AminoacidId>(AMINOACID_CNT - 1);
    }

    bool aminoacid::Aminoacid::same(aminoacid::Aminoacid::AminoacidId a, aminoacid::Aminoacid::AminoacidId b) {
        assert(a >= 0 && a <= AMINOACID_CNT);
        assert(b >= 0 && b <= AMINOACID_CNT);
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
    //hty=Homotyrosine
    //dhb=Threonine

    //TODO finish formuls, sure have mis
    const std::string Aminoacid::AMINOACID_NAMES[Aminoacid::AMINOACID_CNT] = {"trp", "ser", "gly", "uda", "thr",
                                                                                "dhp", "gln", "dab", "arg", "lys",
                                                                                "ala-d", "phe", "val", "cha", "dhpg",
                                                                                "phg", "his", "aeo", "bmt", "hse",
                                                                                "met", "ala", "tcl", "sal", "allothr",
                                                                                "b-ala", "dhb", "ile", "end", "leu",
                                                                                "gua", "hty", "glu", "bht", "hpg",
                                                                                "apa", "pro", "tyr", "hyv", "asn",
                                                                                "cit", "vol", "cys", "asp", "dht",
                                                                                "ahp", "orn", "apc", "abu", "aad",
                                                                                "pip", "dpg", "3-me-glu", "OH-asn",
                                                                                "none"};


    const std::string Aminoacid::FORMULS[Aminoacid::AMINOACID_CNT] = {"C11H12N2O2", "C3H7NO3", "C2H5NO2", "C4H9N3O3", "C4H9NO3",
                                                                        "C3H5NO2", "C5H10N2O3", "C4H10N2O2", "C6H14N4O2", "C6H14N2O2",
                                                                        "C3H7NO2", "C9H11NO2", "C5H11NO2", "C9H17NO2", "C8H9NO4",
                                                                        "C8H9NO2", "C6H9N3O2", "C10H17NO4", "C9H17NO3", "C4H9NO3",
                                                                        "C5H11NO2S", "C3H7NO2", "C6H10Cl3NO2", "C7H6O3", "C4H9NO3",
                                                                        "C3H7NO2", "C7H6O4", "C6H13NO2", "C6H12N4O2", "C6H13NO2",
                                                                        "-", "C10H13NO3", "C5H9NO4", "C9H11NO4", "C8H9NO3",
                                                                        "C4H9NO2", "C5H9NO2", "C9H11NO3", "-", "C4H8N2O3",
                                                                        "C6H13N3O3", "C5H13NO", "C3H7NO2S", "C4H7NO4", "C9H9NO3",
                                                                        "C9H16N2O5", "C5H12N2O2", "-", "C4H9NO2", "C6H11NO4",
                                                                        "C6H11NO2", "C8H9NO4", "C6H11NO4", "C4H8N2O4",
                                                                        "*"};
}