# ChemoKernel
C++ Library to compute graph kernels on molecular graphs. Mainly my PhD work. 

## packages required :
* libx11-dev
* cmake
* liblapack-dev


## To Compile : 

- cd GraphKernels
- mkdir build
- cd build
- cmake ../
- make
- cd ../../
- mkdir build
- cd build
- cmake ../Molecules 
- make

## References
This library is mainly an implementations of these papers:

* Benoit Gaüzère, Pierre-Anthony Grenier, Luc Brun and Didier
  Villemin. Treelet kernel incorporating cyclic, stereo and inter
  pattern information in chemoinformatics. In Pattern Recognition,
  2014.
* Benoit Gaüzère, Luc Brun and Didier Villemin. Two new graphs kernels
  in chemoinformatics. In Pattern Recognition Letters, 2012.
* Benoit Gaüzère, Luc Brun, Didier Villemin. Graph kernel encoding
  substituents' relative positioning . In International Conference on
  Pattern Recognition (ICPR), 2014. 
* Benoit Gaüzère, Luc Brun, Didier Villemin. Relevant Cycle Hypergraph
  Representation for Molecules . In 9th IAPR-TC15 Graph Based
  Representations in Pattern Recognition , 2013. 
* Benoit Gaüzère, Luc Brun, Didier Villemin and Myriam Brun. Graph
  Kernels Based on Relevant Patterns and Cycle Information for
  Chemoinformatics. In International Conference on Pattern Recognition
  (ICPR), 2012. 
*  Benoit Gaüzère, Luc Brun, Didier Villemin and Myriam Brun. Graph
   Kernels: Crossing Information from Different Patterns Using Graph
   Edit Distance. In Structural, Syntactic, and Statistical Pattern
   Recognition (S+SSPR), 2012.
  
It allows also to compute others graph kernels.

## Disclaimer

This is a dirty PhD code. It may bug without any warning. This is
published for experiments reproduction purpose. 

## The Library

This library allows to compute treelet kernel in its various forms. It
take datasets as ct (ChemDraw Connection Table) format, as available [here](https://brunl01.users.greyc.fr/CHEMISTRY/)

To use a kernel, you may be interested in the following bin files:
* moleculesTreeletData : Extract treelet distribution from a dataset
file. Example :
`$>./bin/moleculesTreeletData ~/work/Datasets/Acyclic/ dataset_bps.ds -o data`

This will create 4 files :
 * `data_treelets_types` : Cumulative sum of number of treelets for each kind of pattern 
 * `data_treelets_codes` : Canonical codes associated to each treelet
 * `data_y` : The value to predict associated to each graph of a dataset
 * `data_treelets_distrib` : The number of occurences of each treelet for each graph

* `./bin/moleculesComputeGramMatrice ~/work/Datasets/Acyclic/ dataset_bps.ds -gram gram -k 6`
  Compute a Gram Matrix using treelet kernel and store it in "gram"
  file. This file is a plain text files encoding a NxN Matrix with
  kernel values. You can then use it with any SVM or kernel machine.
* `./bin/moleculesComputeGramMatrices ~/work/Datasets/Acyclic/  dataset_bps.ds -k 0 -K 0 -o data`
  Same as previous but a Gram Matrix is provided for each
  treelet. This is required as an input for Multiple Kernel Learning framework.
* `./bin/moleculesRegression10percent  /home/bgauzere/work/Datasets/Acyclic/ dataset_bps.ds -k 6 -s 0.5 -K
  0 -a 0.01 `
  Compute a 10-fold cross validation for a regression problem with
  various parameters. This set of parameters is optimal for Acyclic
  dataset. Check the output at the end of this file. 
* `bin/moleculesClassification`
  Same but for classification problems

If any problem, email at benoit.gauzere at gmail . com


## Appendix 

  $>bin/moleculesRegression10percent /home/bgauzere/work/Datasets/Acyclic/ dataset_bps.ds -k 6 -s 0.5 -K 0 -a 0.01 

  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:139030
  Molecule 0 : -24.3281752461511 (-23.7)	 --> 	0.628175246151084°
  Molecule 1 : 37.8118531342269 (40)	 --> 	2.18814686577306°
  Molecule 2 : 159.425402411766 (154)	 --> 	5.42540241176633°
  Molecule 3 : 87.9121640991155 (88)	 --> 	0.0878359008844711°
  Molecule 4 : 177.212293826742 (173.7)	 --> 	3.5122938267416°
  Molecule 5 : 80.7560878008165 (81.2)	 --> 	0.443912199183515°
  Molecule 6 : 132.011347942363 (132.5)	 --> 	0.488652057637353°
  Molecule 7 : 128.368044308597 (133.6)	 --> 	5.23195569140347°
  Molecule 8 : 208.94993853767 (211)	 --> 	2.05006146232986°
  Molecule 9 : 105.718453042034 (106.5)	 --> 	0.781546957966228°
  Molecule 10 : 151.377458842552 (146)	 --> 	5.37745884255247°
  Molecule 11 : 152.653291504656 (159)	 --> 	6.34670849534413°
  Molecule 12 : 123.034542308033 (122.2)	 --> 	0.834542308032951°
  Molecule 13 : 149.673273592903 (147)	 --> 	2.67327359290283°
  Molecule 14 : 231.68588195903 (226)	 --> 	5.68588195903047°
  Molecule 15 : 148.318062258627 (139)	 --> 	9.31806225862738°
  Molecule 16 : 157.579312184428 (166.5)	 --> 	8.92068781557211°
  Molecule 17 : 147.956229274723 (162)	 --> 	14.0437707252772°
  Molecule 18 : 209.985021965291 (215)	 --> 	5.0149780347088°
  Mean Prediction Time : 9655
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:407354
  Molecule 19 : 17.1188578598173 (14)	 --> 	3.11885785981733°
  Molecule 20 : 35.2252384547334 (34.6)	 --> 	0.625238454733349°
  Molecule 21 : 157.022467088595 (156)	 --> 	1.02246708859485°
  Molecule 22 : 77.0963601287867 (83)	 --> 	5.90363987121334°
  Molecule 23 : 163.635676718576 (165.5)	 --> 	1.86432328142396°
  Molecule 24 : 87.7791721172915 (93)	 --> 	5.22082788270848°
  Molecule 25 : 116.233815067419 (123.5)	 --> 	7.26618493258125°
  Molecule 26 : 116.233532639876 (120.4)	 --> 	4.16646736012407°
  Molecule 27 : 127.951629999569 (125)	 --> 	2.95162999956902°
  Molecule 28 : 104.718515304628 (101)	 --> 	3.71851530462773°
  Molecule 29 : 139.233449352486 (145)	 --> 	5.76655064751421°
  Molecule 30 : 156.155226405016 (163.5)	 --> 	7.34477359498425°
  Molecule 31 : 122.706998234051 (112)	 --> 	10.7069982340511°
  Molecule 32 : 158.799563106421 (165)	 --> 	6.20043689357914°
  Molecule 33 : 208.868435951635 (215)	 --> 	6.13156404836514°
  Molecule 34 : 163.373760179689 (162)	 --> 	1.37376017968941°
  Molecule 35 : 146.595695380612 (146)	 --> 	0.595695380612455°
  Molecule 36 : 161.08848842402 (162)	 --> 	0.911511575980484°
  Molecule 37 : 210.647236174079 (216)	 --> 	5.35276382592144°
  Mean Prediction Time : 9160
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:651766
  Molecule 38 : 49.6999106195394 (37.3)	 --> 	12.3999106195394°
  Molecule 39 : 30.4092781768384 (32)	 --> 	1.59072182316156°
  Molecule 40 : 156.567779015236 (166)	 --> 	9.43222098476417°
  Molecule 41 : 103.566914175643 (104.5)	 --> 	0.933085824356539°
  Molecule 42 : 174.76648719166 (181)	 --> 	6.23351280833981°
  Molecule 43 : 77.1586380731574 (69)	 --> 	8.15863807315738°
  Molecule 44 : 115.025633943629 (120.3)	 --> 	5.27436605637082°
  Molecule 45 : 119.549886212784 (120)	 --> 	0.450113787215912°
  Molecule 46 : 116.812598609236 (118)	 --> 	1.18740139076445°
  Molecule 47 : 95.283455569108 (99.3)	 --> 	4.01654443089195°
  Molecule 48 : 130.425227178775 (141)	 --> 	10.5747728212247°
  Molecule 49 : 206.931260945931 (229.5)	 --> 	22.5687390540693°
  Molecule 50 : 96.6454730566092 (106)	 --> 	9.35452694339079°
  Molecule 51 : 189.675960264429 (188.9)	 --> 	0.775960264428733°
  Molecule 52 : 181.111475241989 (201)	 --> 	19.8885247580111°
  Molecule 53 : 174.780249491224 (173)	 --> 	1.78024949122408°
  Molecule 54 : 177.916802818836 (165)	 --> 	12.9168028188356°
  Molecule 55 : 174.08849821132 (173.2)	 --> 	0.888498211320382°
  Molecule 56 : 232.292198423554 (240)	 --> 	7.70780157644555°
  Mean Prediction Time : 8797
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:882401
  Molecule 57 : 107.705220947439 (109.7)	 --> 	1.99477905256141°
  Molecule 58 : 56.4297349251341 (63)	 --> 	6.57026507486592°
  Molecule 59 : 169.387392061573 (183)	 --> 	13.6126079384275°
  Molecule 60 : 101.948311973125 (102)	 --> 	0.0516880268754392°
  Molecule 61 : 102.301008779956 (99.5)	 --> 	2.80100877995638°
  Molecule 62 : 86.2004502512818 (86.3)	 --> 	0.0995497487182035°
  Molecule 63 : 152.396904593896 (145)	 --> 	7.39690459389556°
  Molecule 64 : 144.411111298084 (137)	 --> 	7.41111129808434°
  Molecule 65 : 118.637607088911 (117.1)	 --> 	1.53760708891086°
  Molecule 66 : 111.237976678891 (90)	 --> 	21.2379766788913°
  Molecule 67 : 171.830681246313 (171)	 --> 	0.830681246312793°
  Molecule 68 : 142.208269295778 (142)	 --> 	0.208269295778166°
  Molecule 69 : 111.565147361701 (114.5)	 --> 	2.93485263829855°
  Molecule 70 : 166.766973708094 (170)	 --> 	3.23302629190636°
  Molecule 71 : 201.368111579662 (205)	 --> 	3.63188842033816°
  Molecule 72 : 159.112797831761 (159.5)	 --> 	0.387202168239156°
  Molecule 73 : 155.223217744827 (159)	 --> 	3.77678225517292°
  Molecule 74 : 188.471882801062 (186.8)	 --> 	1.67188280106205°
  Mean Prediction Time : 8659
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:1107560
  Molecule 75 : 5.03101639149838 (10.8)	 --> 	5.76898360850162°
  Molecule 76 : 57.4005341117461 (53.5)	 --> 	3.90053411174607°
  Molecule 77 : 71.9074928978082 (70.3)	 --> 	1.60749289780817°
  Molecule 78 : 100.490338650428 (92)	 --> 	8.49033865042831°
  Molecule 79 : 94.0742217883446 (92.3)	 --> 	1.77422178834459°
  Molecule 80 : 83.5748561648457 (82)	 --> 	1.57485616484567°
  Molecule 81 : 148.279948810023 (144.2)	 --> 	4.07994881002324°
  Molecule 82 : 190.953653287028 (195.8)	 --> 	4.84634671297209°
  Molecule 83 : 107.703367889279 (107)	 --> 	0.703367889279306°
  Molecule 84 : 138.707530012616 (137)	 --> 	1.70753001261554°
  Molecule 85 : 170.599146108412 (166)	 --> 	4.59914610841167°
  Molecule 86 : 133.791161890264 (125)	 --> 	8.79116189026365°
  Molecule 87 : 153.652778718362 (151)	 --> 	2.65277871836196°
  Molecule 88 : 182.873182580313 (178)	 --> 	4.8731825803134°
  Molecule 89 : 146.445203518336 (151.5)	 --> 	5.0547964816642°
  Molecule 90 : 156.272636167033 (159.5)	 --> 	3.22736383296737°
  Molecule 91 : 211.988725399332 (195)	 --> 	16.9887253993316°
  Molecule 92 : 175.688411914943 (173)	 --> 	2.68841191494269°
  Mean Prediction Time : 8563
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:1332302
  Molecule 93 : 38.941662029324 (39)	 --> 	0.0583379706759715°
  Molecule 94 : 67.5453677824893 (64.4)	 --> 	3.14536778248925°
  Molecule 95 : 62.2976697386275 (63.6)	 --> 	1.30233026137248°
  Molecule 96 : 109.626781930352 (107.4)	 --> 	2.22678193035182°
  Molecule 97 : 84.2295543238165 (90.1)	 --> 	5.87044567618346°
  Molecule 98 : 102.588580309409 (103)	 --> 	0.411419690590904°
  Molecule 99 : 140.319044522254 (142.8)	 --> 	2.48095547774614°
  Molecule 100 : 177.718209097898 (177.2)	 --> 	0.518209097897568°
  Molecule 101 : 102.648225789402 (102.5)	 --> 	0.14822578940219°
  Molecule 102 : 115.387764945301 (114)	 --> 	1.38776494530066°
  Molecule 103 : 151.115089695928 (155)	 --> 	3.88491030407249°
  Molecule 104 : 131.040930825716 (132)	 --> 	0.959069174283513°
  Molecule 105 : 119.002430575814 (128.5)	 --> 	9.49756942418563°
  Molecule 106 : 159.564452733536 (148.5)	 --> 	11.0644527335355°
  Molecule 107 : 159.29990716573 (165.5)	 --> 	6.20009283427046°
  Molecule 108 : 163.833297257609 (155.5)	 --> 	8.33329725760916°
  Molecule 109 : 208.776634695136 (218)	 --> 	9.22336530486393°
  Molecule 110 : 179.229192732816 (187)	 --> 	7.77080726718401°
  Mean Prediction Time : 8480
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:1554796
  Molecule 111 : 51.8968378417414 (42)	 --> 	9.89683784174142°
  Molecule 112 : 86.9802976060384 (84.7)	 --> 	2.28029760603836°
  Molecule 113 : 52.1855698470543 (52.5)	 --> 	0.314430152945746°
  Molecule 114 : 123.234348536773 (123.2)	 --> 	0.0343485367725265°
  Molecule 115 : 78.165711877495 (80.2)	 --> 	2.03428812250499°
  Molecule 116 : 98.1102462588063 (103.5)	 --> 	5.38975374119369°
  Molecule 117 : 136.368503834857 (132)	 --> 	4.36850383485717°
  Molecule 118 : 185.266553760578 (181)	 --> 	4.26655376057801°
  Molecule 119 : 105.531804810588 (112)	 --> 	6.46819518941189°
  Molecule 120 : 127.486118916511 (126)	 --> 	1.48611891651119°
  Molecule 121 : 150.670137798131 (145)	 --> 	5.67013779813055°
  Molecule 122 : 129.855693552027 (130.5)	 --> 	0.644306447973435°
  Molecule 123 : 110.277848021086 (109.5)	 --> 	0.777848021085873°
  Molecule 124 : 168.104780265245 (165)	 --> 	3.10478026524507°
  Molecule 125 : 150.52337317439 (157)	 --> 	6.47662682561028°
  Molecule 126 : 141.986773101938 (141)	 --> 	0.986773101937899°
  Molecule 127 : 246.777639196288 (250)	 --> 	3.22236080371198°
  Molecule 128 : 173.25259953402 (174)	 --> 	0.747400465980348°
  Mean Prediction Time : 8440
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:1777887
  Molecule 129 : 75.5427115805266 (66.6)	 --> 	8.94271158052659°
  Molecule 130 : 101.099149026897 (95.5)	 --> 	5.5991490268967°
  Molecule 131 : 54.662927581659 (59)	 --> 	4.33707241834104°
  Molecule 132 : 113.828417649685 (112.5)	 --> 	1.32841764968538°
  Molecule 133 : 77.8693062280687 (82)	 --> 	4.13069377193128°
  Molecule 134 : 96.6017751925423 (96)	 --> 	0.601775192542348°
  Molecule 135 : 133.436667638644 (134.2)	 --> 	0.763332361355538°
  Molecule 136 : 186.072789643 (185.9)	 --> 	0.172789643000158°
  Molecule 137 : 93.0012038861267 (97.4)	 --> 	4.39879611387329°
  Molecule 138 : 131.522001456978 (124)	 --> 	7.52200145697796°
  Molecule 139 : 155.685925398326 (159)	 --> 	3.31407460167406°
  Molecule 140 : 117.54380834091 (125)	 --> 	7.45619165908954°
  Molecule 141 : 124.115667000344 (126)	 --> 	1.88433299965583°
  Molecule 142 : 180.216142969628 (177)	 --> 	3.21614296962804°
  Molecule 143 : 130.363797521203 (139)	 --> 	8.63620247879729°
  Molecule 144 : 126.729485223868 (126)	 --> 	0.729485223868267°
  Molecule 145 : 212.673843502615 (235)	 --> 	22.3261564973848°
  Molecule 146 : 191.036584101807 (188.5)	 --> 	2.53658410180677°
  Mean Prediction Time : 8420
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:2002690
  Molecule 147 : 131.771664565311 (135)	 --> 	3.22833543468897°
  Molecule 148 : 89.9335633782254 (92)	 --> 	2.06643662177457°
  Molecule 149 : 59.7913152334319 (59.5)	 --> 	0.291315233431902°
  Molecule 150 : 115.184632601919 (118.5)	 --> 	3.31536739808077°
  Molecule 151 : 90.9088105735942 (91.2)	 --> 	0.291189426405836°
  Molecule 152 : 113.23978842908 (112)	 --> 	1.23978842907971°
  Molecule 153 : 143.066453526927 (137)	 --> 	6.0664535269274°
  Molecule 154 : 172.358353860901 (175.7)	 --> 	3.34164613909877°
  Molecule 155 : 101.639396807234 (91.5)	 --> 	10.1393968072339°
  Molecule 156 : 137.538112876012 (140.5)	 --> 	2.96188712398751°
  Molecule 157 : 141.854642403291 (138)	 --> 	3.85464240329142°
  Molecule 158 : 116.580220532949 (122)	 --> 	5.41977946705137°
  Molecule 159 : 146.01558045644 (147)	 --> 	0.984419543560477°
  Molecule 160 : 165.692901408809 (167)	 --> 	1.30709859119108°
  Molecule 161 : 158.563981747604 (163)	 --> 	4.43601825239631°
  Molecule 162 : 147.962404540204 (164)	 --> 	16.0375954597962°
  Molecule 163 : 191.0966324151 (186.5)	 --> 	4.5966324151002°
  Molecule 164 : 182.324933202676 (199)	 --> 	16.6750667973235°
  Mean Prediction Time : 8395
  Spectrum Kernel : Intersection KernelGram Matrix computed
  Learning Time:2226576
  Molecule 165 : 144.111585362524 (148.5)	 --> 	4.38841463747562°
  Molecule 166 : 86.3777698939996 (84.4)	 --> 	1.97776989399962°
  Molecule 167 : 62.8321864701092 (55.2)	 --> 	7.63218647010916°
  Molecule 168 : 102.731389204222 (101.5)	 --> 	1.23138920422204°
  Molecule 169 : 91.7242795866036 (91.5)	 --> 	0.224279586603615°
  Molecule 170 : 106.888999823352 (104)	 --> 	2.88899982335201°
  Molecule 171 : 141.342264167235 (139)	 --> 	2.34226416723547°
  Molecule 172 : 198.541406185377 (186)	 --> 	12.5414061853771°
  Molecule 173 : 94.5260562022613 (87.6)	 --> 	6.92605620226129°
  Molecule 174 : 148.262019269227 (157.5)	 --> 	9.23798073077347°
  Molecule 175 : 145.147244873906 (142)	 --> 	3.14724487390643°
  Molecule 176 : 133.546473813741 (121)	 --> 	12.5464738137407°
  Molecule 177 : 155.649958864346 (158)	 --> 	2.3500411356537°
  Molecule 178 : 190.121291467466 (195)	 --> 	4.87870853253386°
  Molecule 179 : 162.469799984499 (153.5)	 --> 	8.96979998449919°
  Molecule 180 : 161.293419941354 (163)	 --> 	1.70658005864649°
  Molecule 181 : 145.313801441483 (156.5)	 --> 	11.1861985585172°
  Molecule 182 : 228.02931998909 (228)	 --> 	0.0293199890904248°
  Mean Prediction Time : 8472

	Average error : 4.69482116282133
	Standard Deviation : 6.45314284880126
	R : 0.990933085827914


