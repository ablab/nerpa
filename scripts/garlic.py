import os

genomeList = [
"GCF_001654515.1_ASM165451v1_genomic",
"GCF_000819665.1_ASM81966v1_genomic",
"GCF_001590605.1_ASM159060v1_genomic",
"GCF_000012445.1_ASM1244v1_genomic",
"GCF_000250655.1_ASM25065v1_genomic",
"GCF_000336465.1_ASM33646v1_genomic",
"GCF_000009225.2_ASM922v1_genomic",
"GCF_000026105.1_ASM2610v1_genomic",
"GCF_002007785.1_ASM200778v1_genomic",
"GCF_002117165.1_ASM211716v1_genomic",
"GCF_000283695.1_ASM28369v1_genomic",
"GCF_000828695.1_ASM82869v1_genomic",
"GCF_001612705.1_ASM161270v1_genomic",
"GCF_002208745.1_ASM220874v1_genomic",
"GCF_001553895.1_ASM155389v1_genomic",
"GCF_002157265.1_ASM215726v1_genomic",
"GCF_000012265.1_ASM1226v1_genomic",
"GCF_001579805.1_ASM157980v1_genomic",
"GCF_001708485.1_ASM170848v1_genomic",
"GCF_002202055.1_ASM220205v1_genomic",
"GCF_001874385.1_ASM187438v1_genomic",
"GCF_002192235.1_ASM219223v1_genomic",
"GCF_001889385.1_ASM188938v1_genomic",
"GCF_002005345.1_ASM200534v1_genomic",
"GCF_000015785.1_ASM1578v1_genomic",
"GCF_000195515.1_ASM19551v1_genomic",
"GCF_000833005.1_ASM83300v1_genomic",
"GCF_001593765.1_ASM159376v1_genomic",
"GCF_002057535.1_ASM205753v1_genomic",
"GCF_001593395.2_ASM159339v2_genomic",
"GCF_001647715.1_ASM164771v1_genomic",
"GCF_001705195.1_ASM170519v1_genomic",
"GCF_001808235.1_ASM180823v1_genomic",
"GCF_002205715.1_ASM220571v1_genomic",
"GCF_000408885.1_ASM40888v1_genomic",
"GCF_000341875.1_ASM34187v1_genomic",
"GCF_000827045.1_ASM82704v1_genomic",
"GCF_000973585.1_ASM97358v1_genomic",
"GCF_001586105.1_ASM158610v1_genomic",
"GCF_001685645.1_ASM168564v1_genomic",
"GCF_001687745.1_ASM168774v1_genomic",
"GCF_000455565.1_ucmb5033v1_genomic",
"GCF_000494835.1_ASM49483v1_genomic",
"GCF_000769555.1_ASM76955v1_genomic",
"GCF_000772165.1_ASM77216v1_genomic",
"GCF_000987825.1_ASM98782v1_genomic",
"GCF_001465815.1_ASM146581v1_genomic",
"GCF_001534785.1_ASM153478v1_genomic",
"GCF_001723585.1_ASM172358v1_genomic",
"GCF_001752685.1_ASM175268v1_genomic",
"GCF_002055965.1_ASM205596v1_genomic",
"GCF_002173495.1_ASM217349v1_genomic",
"GCF_000338735.1_ASM33873v1_genomic",
"GCF_000397205.1_ASM39720v1_genomic",
"GCF_000455585.1_ucmb5113v1_genomic",
"GCF_000583065.1_ASM58306v1_genomic",
"GCF_000772205.1_ASM77220v1_genomic",
"GCF_000782835.1_ASM78283v1_genomic",
"GCF_000789275.1_ASM78927v1_genomic",
"GCF_000789295.1_ASM78929v1_genomic",
"GCF_000953615.1_BS49Ch_genomic",
"GCF_000973605.1_ASM97360v1_genomic",
"GCF_001023595.1_ASM102359v1_genomic",
"GCF_001483885.1_ASM148388v1_genomic",
"GCF_001602135.1_ASM160213v1_genomic",
"GCF_001636055.1_ASM163605v1_genomic",
"GCF_001704095.1_ASM170409v1_genomic",
"GCF_001857985.1_ASM185798v1_genomic",
"GCF_001889285.1_ASM188928v1_genomic",
"GCF_002006545.1_ASM200654v1_genomic",
"GCF_002072735.1_ASM207273v1_genomic",
"GCF_002096095.1_ASM209609v1_genomic",
"GCF_002113805.1_ASM211380v1_genomic",
"GCF_002201875.1_ASM220187v1_genomic",
"GCF_002209305.1_ASM220930v1_genomic",
"GCF_000876525.1_ASM87652v1_genomic",
"GCF_001597265.1_ASM159726v1_genomic",
"GCF_002202015.1_ASM220201v1_genomic",
"GCF_000815145.1_ASM81514v1_genomic",
"GCF_002068155.1_ASM206815v1_genomic",
"GCF_002072695.1_ASM207269v1_genomic",
"GCF_000165925.1_ASM16592v1_genomic",
"GCF_000186745.1_ASM18674v1_genomic",
"GCF_000204275.1_ASM20427v1_genomic",
"GCF_000319475.1_ASM31947v1_genomic",
"GCF_000830075.1_ASM83007v1_genomic",
"GCF_000953355.1_XNC2_genomic",
"GCF_001536925.1_BMC_PRJEB12378v1_genomic",
"GCF_000009045.1_ASM904v1_genomic",
"GCF_000493375.1_BAPNAU_WGS_genomic",
"GCF_000685725.1_ASM68572v1_genomic",
"GCF_000699465.1_ASM69946v1_genomic",
"GCF_000835145.1_ASM83514v1_genomic",
"GCF_000877815.1_ASM87781v1_genomic",
"GCF_000952895.1_BS34ACh_genomic",
"GCF_000971925.1_ASM97192v1_genomic",
"GCF_000973485.1_ASM97348v1_genomic",
"GCF_000988345.1_ASM98834v1_genomic",
"GCF_001015095.1_ASM101509v1_genomic",
"GCF_001541905.1_ASM154190v1_genomic",
"GCF_001703495.1_ASM170349v1_genomic",
"GCF_001854345.1_ASM185434v1_genomic",
"GCF_001999205.1_ASM199920v1_genomic",
"GCF_002105595.1_ASM210559v1_genomic",
"GCF_002142595.1_ASM214259v1_genomic",
"GCF_000008425.1_ASM842v1_genomic",
"GCF_000146565.1_ASM14656v1_genomic",
"GCF_000221645.1_ASM22164v1_genomic",
"GCF_000227485.1_ASM22748v1_genomic",
"GCF_000293765.1_ASM29376v1_genomic",
"GCF_000344745.1_ASM34474v1_genomic",
"GCF_000497485.1_ASM49748v1_genomic",
"GCF_000597985.1_ASM59798v1_genomic",
"GCF_000737405.1_ASM73740v1_genomic",
"GCF_000772125.1_ASM77212v1_genomic",
"GCF_000800255.1_ASM80025v1_genomic",
"GCF_000978495.1_ASM97849v1_genomic",
"GCF_001037985.1_ASM103798v1_genomic",
"GCF_001565875.1_ASM156587v1_genomic",
"GCF_001593785.1_ASM159378v1_genomic",
"GCF_001594225.1_ASM159422v1_genomic",
"GCF_001660525.1_ASM166052v1_genomic",
"GCF_001719185.1_ASM171918v1_genomic",
"GCF_001720505.1_ASM172050v1_genomic",
"GCF_002163815.1_ASM216381v1_genomic",
"GCF_002173615.1_ASM217361v1_genomic",
"GCF_002173635.1_ASM217363v1_genomic",
"GCF_002173695.1_ASM217369v1_genomic",
"GCF_002173715.1_ASM217371v1_genomic",
"GCF_002201955.1_ASM220195v1_genomic",
"GCF_002202035.1_ASM220203v1_genomic",
"GCF_000739105.1_ASM73910v1_genomic",
"GCF_001889625.1_ASM188962v1_genomic",
"GCF_000203835.1_ASM20383v1_genomic",
"GCF_000508265.1_ASM50826v1_genomic",
"GCF_001922005.1_ASM192200v1_genomic",
"GCF_000242855.2_ASM24285v2_genomic",
"GCF_000196735.1_ASM19673v1_genomic",
"GCF_000227465.1_ASM22746v1_genomic",
"GCF_001596755.1_ASM159675v1_genomic",
"GCF_001719225.1_ASM171922v1_genomic",
"GCF_000011645.1_ASM1164v1_genomic",
"GCF_000017885.4_ASM1788v4_genomic",
"GCF_000164985.3_ASM16498v2_genomic",
"GCF_000237325.1_ASM23732v1_genomic",
"GCF_000259365.1_ASM25936v1_genomic",
"GCF_000523045.1_ASM52304v1_genomic",
"GCF_000590455.1_ASM59045v1_genomic",
"GCF_000699525.1_ASM69952v1_genomic",
"GCF_000706705.1_ASM70670v1_genomic",
"GCF_000827065.1_ASM82706v1_genomic",
"GCF_001431145.1_ASM143114v1_genomic",
"GCF_001534745.1_ASM153474v1_genomic",
"GCF_001548215.1_ASM154821v1_genomic",
"GCF_001578165.1_ASM157816v1_genomic",
"GCF_001596055.1_ASM159605v1_genomic",
"GCF_001604995.1_ASM160499v1_genomic",
"GCF_001704195.1_ASM170419v1_genomic",
"GCF_001719045.1_ASM171904v1_genomic",
"GCF_001746575.1_ASM174657v1_genomic",
"GCF_001890405.1_ASM189040v1_genomic",
"GCF_001896025.1_ASM189602v1_genomic",
"GCF_001922145.1_ASM192214v1_genomic",
"GCF_002074075.1_ASM207407v1_genomic",
"GCF_002074115.1_ASM207411v1_genomic",
"GCF_002074135.1_ASM207413v1_genomic",
"GCF_002173675.1_ASM217367v1_genomic",
"GCF_002174255.1_ASM217425v1_genomic",
"GCF_002201995.1_ASM220199v1_genomic",
"GCF_000010165.1_ASM1016v1_genomic",
"GCF_000020025.1_ASM2002v1_genomic",
"GCF_000146875.3_ASM14687v2_genomic",
"GCF_000209795.2_ASM20979v2_genomic",
"GCF_000321395.1_ASM32139v1_genomic",
"GCF_000359525.1_ASM35952v1_genomic",
"GCF_000498975.2_ASM49897v2_genomic",
"GCF_000517305.1_ASM51730v1_genomic",
"GCF_000800825.1_ASM80082v1_genomic",
"GCF_000829415.1_ASM82941v1_genomic",
"GCF_000968195.1_ASM96819v1_genomic",
"GCF_000972685.1_ASM97268v1_genomic",
"GCF_001007005.1_Pseudomonas_sp._strain_CCOS_191_genomic",
"GCF_001191605.1_ASM119160v1_genomic",
"GCF_001272655.2_ASM127265v2_genomic",
"GCF_001307275.1_ASM130727v1_genomic",
"GCF_001431765.1_ASM143176v1_genomic",
"GCF_001431785.1_ASM143178v1_genomic",
"GCF_001442805.1_ASM144280v1_genomic",
"GCF_001511775.1_YEL.embl.gz_genomic",
"GCF_001515585.1_ASM151558v2_genomic",
"GCF_001577385.1_ASM157738v1_genomic",
"GCF_001578205.1_ASM157820v1_genomic",
"GCF_001595725.1_ASM159572v1_genomic",
"GCF_001655595.1_ASM165559v1_genomic",
"GCF_001700735.1_ASM170073v1_genomic",
"GCF_001704975.1_ASM170497v1_genomic",
"GCF_001726125.1_ASM172612v1_genomic",
"GCF_001747445.1_ASM174744v1_genomic",
"GCF_001857925.1_ASM185792v1_genomic",
"GCF_001874405.1_ASM187440v2_genomic",
"GCF_001874425.2_ASM187442v2_genomic",
"GCF_001895885.1_ASM189588v1_genomic",
"GCF_001902555.1_ASM190255v1_genomic",
"GCF_001908475.1_ASM190847v1_genomic",
"GCF_001938665.1_ASM193866v1_genomic",
"GCF_001938705.1_ASM193870v1_genomic",
"GCF_002021815.1_ASM202181v1_genomic",
"GCF_002074095.1_ASM207409v1_genomic",
"GCF_002077215.1_ASM207721v1_genomic",
"GCF_002215075.1_ASM221507v1_genomic",
"GCF_000511405.1_ASM51140v1_genomic",
"GCF_000219535.2_ASM21953v3_genomic",
"GCF_002082155.1_ASM208215v1_genomic",
"GCF_000258535.2_ASM25853v2_genomic",
"GCF_000218915.1_ASM21891v1_genomic",
"GCF_000400635.2_ASM40063v2_genomic",
"GCF_002162355.1_ASM216235v1_genomic",
"GCF_000949425.1_ASM94942v1_genomic",
"GCF_000452705.1_ASM45270v3_genomic",
"GCF_000219105.1_ASM21910v1_genomic",
"GCF_000520015.2_ASM52001v2_genomic",
"GCF_001547935.1_ASM154793v1_genomic",
"GCF_001482725.1_ASM148272v1_genomic",
"GCF_001767395.1_ASM176739v1_genomic",
"GCF_000271665.2_ASM27166v2_genomic",
"GCF_000988565.1_ASM98856v1_genomic",
"GCF_000012245.1_ASM1224v1_genomic",
"GCF_000014985.1_ASM1498v1_genomic",
"GCF_000027225.1_ASM2722v1_genomic",
"GCF_000196475.1_ASM19647v1_genomic",
"GCF_000221045.1_ASM22104v1_genomic",
"GCF_000245355.1_ASM24535v1_genomic",
"GCF_000255295.1_ASM25529v1_genomic",
"GCF_000340845.1_ASM34084v1_genomic",
"GCF_000513215.1_DB11_genomic",
"GCF_000747565.1_ASM74756v1_genomic",
"GCF_000756615.1_ASM75661v1_genomic",
"GCF_000829075.1_ASM82907v1_genomic",
"GCF_000988395.1_ASM98839v1_genomic",
"GCF_000988485.1_ASM98848v1_genomic",
"GCF_001281365.1_ASM128136v1_genomic",
"GCF_001584145.1_ASM158414v1_genomic",
"GCF_001642805.2_ASM164280v2_genomic",
"GCF_001936215.1_ASM193621v1_genomic",
"GCF_002196515.1_ASM219651v1_genomic",
"GCF_000943515.2_ASM94351v2_genomic",
"GCF_000250675.2_ASM25067v3_genomic",
"GCF_002209125.1_ASM220912v1_genomic",
"GCF_000364805.1_ASM36480v1_genomic",
"GCF_001704275.1_ASM170427v1_genomic",
"GCF_000092385.1_ASM9238v1_genomic",
"GCF_000176115.2_ASM17611v2_genomic",
"GCF_000177195.2_ASM17719v2_genomic",
"GCF_000269985.1_ASM26998v1_genomic",
"GCF_000732925.1_ASM73292v1_genomic",
"GCF_000954115.1_ASM95411v1_genomic",
"GCF_001406115.1_Streptomyces_venezuelae_ATCC_15439_genomic",
"GCF_001443625.1_ASM144362v1_genomic",
"GCF_001579845.1_ASM157984v1_genomic",
"GCF_001620005.1_ASM162000v1_genomic",
"GCF_001620305.1_ASM162030v1_genomic",
"GCF_001975025.1_ASM197502v1_genomic",
"GCF_001984015.1_ASM198401v1_genomic"
]

molList = [
"BOW79-M",
"PVS77-L",
"antimarin2012_46572",
"mibig_798",
"mibig_682",
"HLY99-V",
"HLY96-S",
"QYK92-E",
"QYN99-A",
"QCP30-L",
"QCW69-M",
"QYN98-Z",
"BCN63-E",
"LRY16-C",
"antimarin2012_27348",
"antimarin2012_47286",
"LQY02-O",
"MWQ86-D",
"MWQ88-F",
"antimarin2012_27351",
"QMK48-P",
"QMN08-S",
"PJQ94-Y",
"streptomedb.3052",
"streptomedb.2916",
"PJQ77-V",
"QMN09-T",
"CLL31-F",
"CLP93-T",
"CLL30-E",
"mibig_92",
"QRY79-A",
"CMD95-U",
"antimarin2012_35354",
"CLL36-K",
"antimarin2012_47853",
"antimarin2012_27355",
"BOV12-I",
"BOV11-H",
"BOV10-G",
"mibig_466",
"CLL34-I",
"antimarin2012_27353",
"antimarin2012_27354",
"streptomedb.3001",
"streptomedb.3018",
"streptomedb.2987",
"HLY92-O",
"HLY91-N",
"streptomedb.3019",
"BCN67-I",
"BCN65-G",
"antimarin2012_10211",
"BOV09-M",
"PVY02-P",
"JJN76-R",
"QSG56-M",
"JXY52-E",
"KYS61-U",
"KYS58-Y",
"KYS53-T",
"KYS54-U",
"BOV06-J",
"MRX90-A",
"BOV01-E",
"streptomedb.3906",
"streptomedb.3905",
"streptomedb.3904",
"JJM23-U",
"mibig_724",
"mibig_154",
"antimarin2012_10213",
"antimarin2012_10212",
"MWQ85-C",
"BCN69-K",
"MWQ89-G",
"MWQ87-E",
"JXY56-I",
"streptomedb.3903",
"NZJ30-F",
"NZJ35-K",
"NZJ26-I",
"NZJ23-F",
"antimarin2012_27002",
"HKY27-R",
"antimarin2012_32086",
"antimarin2012_25889",
"antimarin2012_33077",
"streptomedb.3025",
"streptomedb.3014",
"mibig_723",
"streptomedb.3000",
"streptomedb.3026",
"PJQ96-A",
"streptomedb.3024",
"KYS56-W",
"KYS59-Z",
"JNJ55-S",
"antimarin2012_47332",
"QYO04-Z",
"QYO01-W",
"QYO00-V",
"QYO03-Y",
"QYO09-E",
"antimarin2012_44712",
"KCH49-V",
"antimarin2012_49214",
"QYK70-W",
"antimarin2012_54195",
"antimarin2012_45491",
"PGN93-N",
"antimarin2012_45497",
"JOW50-H",
"antimarin2012_47554",
"JOW91-U",
"antimarin2012_25137",
"mibig_686",
"POG24-O",
"PNZ76-P",
"POG25-P",
"antimarin2012_27347",
"QCP28-Q",
"antimarin2012_58939",
"QCP29-R",
"antimarin2012_48367",
"antimarin2012_48368",
"antimarin2012_27350",
"JTO33-Z",
"antimarin2012_27349",
"antimarin2012_44297",
"QMN74-J",
"QKX14-B",
"BCQ66-W",
"antimarin2012_25136",
"QKX11-Y",
"QMN75-K",
"antimarin2012_44118",
"antimarin2012_44120",
"GYF27-J",
"QKX13-A",
"QKX12-Z",
"antimarin2012_38808",
"CKC12-A",
"MBT08-F",
"QFP97-F",
"CKC11-Z",
"CKC14-C",
"CKC17-F",
"mibig_905",
"GQY87-O",
"antimarin2012_37242",
"antimarin2012_37236",
"antimarin2012_37246",
"antimarin2012_37240",
"antimarin2012_37238",
"antimarin2012_37244",
"antimarin2012_24106",
"antimarin2012_37234",
"antimarin2012_37243",
"antimarin2012_37237",
"GQY89-Q",
"antimarin2012_37241",
"antimarin2012_37245",
"antimarin2012_37239",
"antimarin2012_37235",
"GQY90-K",
"antimarin2012_24097",
"antimarin2012_24100",
"antimarin2012_28820",
"OKQ78-S",
"OKQ80-N",
"OKQ83-Q",
"OKQ82-P",
"QYK71-X",
"antimarin2012_24096",
"antimarin2012_28815",
"antimarin2012_28817",
"NRT47-J",
"NRT46-I",
"JMS03-N",
"mibig_597",
"antimarin2012_14887",
"QWK05-S",
"antimarin2012_50943",
"LCX55-H",
"QWK04-R",
"antimarin2012_37581",
"antimarin2012_51000",
"antimarin2012_48252",
"mibig_601",
"antimarin2012_35777",
"LQY20-S",
"antimarin2012_35880",
"antimarin2012_37580",
"antimarin2012_20281",
"antimarin2012_52500",
"CKS73-V",
"KZY54-F",
"antimarin2012_54635",
"antimarin2012_44988",
"antimarin2012_45007",
"antimarin2012_35881",
"QVX16-C",
"antimarin2012_37584",
"CVY76-B",
"antimarin2012_18345",
"antimarin2012_49375",
"OJL16-S",
"QVX17-D",
"QWK35-B",
"antimarin2012_35879",
"antimarin2012_35475",
"antimarin2012_48029",
"antimarin2012_54646",
"QYL64-C",
"antimarin2012_32179",
"antimarin2012_37582",
"QMD52-D",
"mibig_614",
"antimarin2012_33986",
"antimarin2012_14888",
"QVX15-B",
"antimarin2012_37587",
"antimarin2012_51455",
"antimarin2012_37583",
"antimarin2012_7778",
"JRD64-S",
"LFH04-K",
"MSZ65-N",
"antimarin2012_37825",
"CRO40-M",
"JXY58-K",
"PVC12-M",
"antimarin2012_56362",
"antimarin2012_47768",
"antimarin2012_4425",
"PVK29-K",
"NRX07-R",
"antimarin2012_48308",
"antimarin2012_54056",
"PKS82-K",
"antimarin2012_14748",
"antimarin2012_19563",
"NRX09-T",
"GYC74-G",
"antimarin2012_51851",
"PKS83-L",
"QFC32-V",
"HSZ94-S",
"antimarin2012_51001",
"antimarin2012_52453",
"antimarin2012_52120",
"antimarin2012_52502",
"CHC01-B",
"CHB53-N"
]



garlicPATH = "/home/dereplicator/kolga/soft/garlic.jar"
predictionJsonPATH = "/home/dereplicator/kolga/data/bacteria_complete/json/"
molJsonPATH = "/home/dereplicator/kolga/data/bacteria_complete/PNP/pnpdatabase/json/"

finfo = open('/home/dereplicator/kolga/data/bacteria_complete/PNP/pnpdatabase/library.info')
infos = finfo.readlines()
finfo.close()

#for i in range(len(infos)):
#    info = ((infos[i].split(' ')[0]).split('/')[-1]).split('.')[0]
#    if (info == "streptomedb"):
#        info = info + "." + ((infos[i].split(' ')[0]).split('/')[-1]).split('.')[1]
#    molList.append(info)

for i in range(len(genomeList)):
    predPath = predictionJsonPATH + genomeList[i] + ".fna.json"
    for j in range(len(molList)):
        molPath = molJsonPATH + molList[j] + "/json/molecule.json"
        os.system("mkdir " + genomeList[i] + "_" + molList[j])
        os.system("java -jar " + garlicPATH + " -q " + predPath + " -s " + molPath + " -o " + genomeList[i] + "_" + molList[j])
