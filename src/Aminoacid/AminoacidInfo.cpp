//
// Created by olga on 21.03.19.
//

#include <fstream>
#include <sstream>
#include <vector>
#include "AminoacidInfo.h"

namespace aminoacid {
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
    std::vector<std::string> AminoacidInfo::AMINOACID_NAMES = {"trp", "ser", "gly", "uda", "thr",
                                                                                      "dhp", "gln", "dab", "arg", "lys",
                                                                                      "ala-d", "phe", "val", "cha", "dhpg",
                                                                                      "phg", "his", "aeo", "bmt", "hse",
                                                                                      "met", "ala", "tcl", "sal", "allothr",
                                                                                      "b-ala", "dhb", "ile", "end", "leu",
                                                                                      "gua", "hty", "glu", "bht", "hpg",
                                                                                      "apa", "pro", "tyr", "hyv", "asn",
                                                                                      "cit", "vol", "cys", "asp", "dht",
                                                                                      "ahp", "orn", "apc", "abu", "aad",
                                                                                      "pip", "dpg", "none"};

    std::vector<std::string> AminoacidInfo::NAME_ID = {"trp", "ser", "gly", "uda", "thr",
                                                               "dhp", "gln", "dab", "arg", "lys",
                                                               "ala-d", "phe", "val", "cha", "dhpg",
                                                               "phg", "his", "aeo", "bmt", "hse",
                                                               "met", "ala", "tcl", "sal", "allothr",
                                                               "b-ala", "dhb", "ile", "end", "leu",
                                                               "gua", "hty", "glu", "bht", "hpg",
                                                               "apa", "pro", "tyr", "hyv", "asn",
                                                               "cit", "vol", "cys", "asp", "dht",
                                                               "ahp", "orn", "apc", "abu", "aad",
                                                               "pip", "dpg", "none"};


    std::vector<Formula> AminoacidInfo::FORMULS = {Formula("C11H12N2O2"), Formula("C3H7NO3"), Formula("C2H5NO2"), Formula("C4H9N3O3"), Formula("C4H9NO3"),
                                                                          Formula("C3H5NO2"), Formula("C5H10N2O3"), Formula("C4H10N2O2"), Formula("C6H14N4O2"), Formula("C6H14N2O2"),
                                                                          Formula("C3H7NO2"), Formula("C9H11NO2"), Formula("C5H11NO2"), Formula("C9H17NO2"), Formula("C8H9NO4"),
                                                                          Formula("C8H9NO2"), Formula("C6H9N3O2"), Formula("C10H17NO4"), Formula("C9H17NO3"), Formula("C4H9NO3"),
                                                                          Formula("C5H11NO2S"), Formula("C3H7NO2"), Formula("C6H10Cl3NO2"), Formula("C7H6O3"), Formula("C4H9NO3"),
                                                                          Formula("C3H7NO2"), Formula("C7H6O4"), Formula("C6H13NO2"), Formula("C6H12N4O2"), Formula("C6H13NO2"),
                                                                          Formula("C179"), Formula("C10H13NO3"), Formula("C5H9NO4"), Formula("C9H11NO4"), Formula("C8H9NO3"),
                                                                          Formula("C4H9NO2"), Formula("C5H9NO2"), Formula("C9H11NO3"), Formula("H179"), Formula("C4H8N2O3"),
                                                                          Formula("C6H13N3O3"), Formula("C5H13NO"), Formula("C3H7NO2S"), Formula("C4H7NO4"), Formula("C9H9NO3"),
                                                                          Formula("C9H16N2O5"), Formula("C5H12N2O2"), Formula("S179"), Formula("C4H9NO2"), Formula("C6H11NO4"),
                                                                          Formula("C6H11NO2"), Formula("C8H9NO4"), Formula()};

    std::vector<double> AminoacidInfo::LogP = {};

    int AminoacidInfo::AMINOACID_CNT = int(AMINOACID_NAMES.size());
    const double AminoacidInfo::DEFAULT_LOG_P = -6.64;

    void AminoacidInfo::init(std::string filename, std::string predictor) {
        std::ifstream in(filename);
        std::string buffer;
        getline(in, buffer);
        std::stringstream ss(buffer);
        std::vector<std::string> header;
        std::string s;
        int nameid = 0;
        int formulaid = 0;
        int logpid = 0;
        int i = 0;
        while (getline(ss, s, '\t')) {
            if (s.size() > 0 && (s.back() == '\r' || s.back() == '\b')) {
                s.pop_back();
            }

            header.push_back(s);
            if (s == "Formula") {
                formulaid = i;
            }
            if (s == predictor) {
                nameid = i;
            }
            if (s == "NRPSPred2_LogP") {
                logpid = i;
            }
            i += 1;
        }

        FORMULS.resize(0);
        AMINOACID_NAMES.resize(0);
        NAME_ID.resize(0);
        LogP.resize(0);

        while (getline(in, buffer)) {
            std::stringstream ss1(buffer);
            std::vector<std::string> parts;
            std::string tmp;
            while (getline(ss1, tmp, '\t')) {
                if (tmp.size() > 0 && (tmp.back() == '\r' || tmp.back() == '\b')) {
                    tmp.pop_back();
                }
                parts.push_back(tmp);
            }
            if (parts[nameid] != "-") {
                NAME_ID.push_back(parts[0]);
                AMINOACID_NAMES.push_back(parts[nameid]);
                if (parts[formulaid] != "-") {
                    FORMULS.push_back(Formula(parts[formulaid]));
                } else {
                    FORMULS.push_back(Formula());
                }
                if (parts[logpid] != "-") {
                    std::stringstream ss(parts[logpid]);
                    double lgp;
                    ss >> lgp;
                    LogP.push_back(lgp);
                } else {
                    LogP.push_back(DEFAULT_LOG_P);
                }
            }
        }

        FORMULS.push_back(Formula());
        AMINOACID_NAMES.push_back("none");
        NAME_ID.push_back("none");
        LogP.push_back(DEFAULT_LOG_P);

        AMINOACID_CNT = AMINOACID_NAMES.size();
        in.close();
    }

    int AminoacidInfo::getIdByNameId(std::string name) {
        for (int i = 0; i < AMINOACID_CNT; ++i) {
            if (NAME_ID[i] == name) {
                return i;
            }
        }
        return AMINOACID_CNT;
    }
}
