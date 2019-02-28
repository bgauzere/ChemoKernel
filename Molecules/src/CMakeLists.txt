set(EXECUTABLE_OUTPUT_PATH ../bin)

# Inclusion of the Molecules library headers
include_directories(../include)

# Inclusion of the headers of Pandore, Gtb, GraphKernels,
# CImg, LibSVM, and RIlib

include_directories(../../GraphKernels/include/pandore)
include_directories(../../GraphKernels/include/gtb)
include_directories(../../GraphKernels/include/GraphKernels)
include_directories(../../GraphKernels/include/cimg)
include_directories(../../GraphKernels/include/libsvm)

# The directory in which are located the libraries

# link_directories(/usr/local/lib)
link_directories(../../GraphKernels/build/src/pandore)
link_directories(../../GraphKernels/build/src/gtb)
link_directories(../../GraphKernels/build/src/GraphKernels)
link_directories(../../GraphKernels/build/src/libsvm)
#link_directories(~/dev/tbb40_src/build/linux_intel64_gcc_cc4.6.1_libc2.13_kernel3.0.0_release)
link_directories(/usr/usr/local/lib)
#link_directories(../../GraphKernels/build/src/cimg)
file(GLOB_RECURSE headers
  ../include/*.h
)

# Creation of the library 'Molecules'

add_library(Molecules STATIC
	${headers}
	GraphFileKernel.cpp
	MoleculeGraph.cpp
	MoleculeGraphUnlabeledSpectrum.cpp
	MoleculeGraphLabeledSpectrum.cpp
	MoleculeGraphCycle.cpp
	MoleculeGraphContractedCycleHypergraph.cpp
	MoleculeGraphAugmentedCycle.cpp
	MoleculesDataset.cpp
	MoleculeGraphEditDistanceMCS.cpp
	MoleculeGraphEditDistanceV2.cpp
	MoleculeGraphEditDistance.cpp
	GraphletCountKernel.cpp
	gnuplot_i.c
	utils.cpp
	string_utils.cpp
	TreeletKernel.cpp
	TreeletGaussianKernel.cpp
	TreeletCompleteGaussianKernel.cpp
	TreeletInnerProductKernel.cpp
	TreeletInnerGaussianKernel.cpp
	TreeletPolynomialKernel.cpp
	TanimotoKernel.cpp
	MinMaxKernel.cpp
	TreeletBinaryKernel.cpp
	TreeletIntersectionKernel.cpp
	TreeletRandomKernel.cpp
	TreeletJaccardKernel.cpp
	MoleculeSpectrumKernel.cpp
	WeightedRegression.cpp
	KPrototypeEditDistance.cpp
	TreeletEnumerator.cpp
	TreeletEnumeratorAugmentedCycles.cpp
	CycleKernel.cpp
	AugmentedCycleKernel.cpp
	TreeletCycleKernel.cpp
	CombinedKernel.cpp
	HorvathKernel.cpp
	InterTreeletKernel.cpp
	IntersectionKernel.cpp
	GaussianKernel.cpp
	InnerProductKernel.cpp
	TreeletEditDistance.cpp
	VectorFileKernel.cpp
	ContractedCycleKernel.cpp
	StereoName.cpp
	SubTreeKernel.cpp
	)

# Creation of the executables

add_executable(moleculesRegression moleculesRegression.cpp)

target_link_libraries(moleculesRegression
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesWeightedRegression moleculesWeightedRegression.cpp)
	
target_link_libraries(moleculesWeightedRegression
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegression10percent moleculesRegression10percent.cpp)
	
target_link_libraries(moleculesRegression10percent
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)
add_executable(moleculesRegression10percent_i moleculesRegression10percent_i.cpp)
	
target_link_libraries(moleculesRegression10percent_i
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

# add_executable(moleculesRegression10percent_tbb moleculesRegression10percent_tbb.cpp)
	
# target_link_libraries(moleculesRegression10percent_tbb
#   Molecules
#   GraphKernels
#   svm
#   pandore
#   gtb
# )

add_executable(moleculesClassification moleculesClassification.cpp)
	
target_link_libraries(moleculesClassification
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)


add_executable(moleculesCalculNbAtoms moleculesCalculNbAtoms.cpp)
	
target_link_libraries(moleculesCalculNbAtoms
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesCalculNbCycles moleculesCalculNbCycles.cpp)
	
target_link_libraries(moleculesCalculNbCycles
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)


add_executable(moleculesCorrelation moleculesCorrelation.cpp)
	
target_link_libraries(moleculesCorrelation
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesTestInversion moleculesTestInversion.cpp)
	
target_link_libraries(moleculesTestInversion
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)


add_executable(moleculesRegressionCluster moleculesRegressionCluster.cpp)
	
target_link_libraries(moleculesRegressionCluster
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegressionFastAdd moleculesRegressionFastAdd.cpp)
	
target_link_libraries(moleculesRegressionFastAdd
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegressionSlowAdd moleculesRegressionSlowAdd.cpp)
	
target_link_libraries(moleculesRegressionSlowAdd
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegressionSlowAdd2 moleculesRegressionSlowAdd2.cpp)
	
target_link_libraries(moleculesRegressionSlowAdd2
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegressionDirectory moleculesRegressionDirectory.cpp)
	
target_link_libraries(moleculesRegressionDirectory
  pandore
  gtb
  GraphKernels
  Molecules
  svm
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegressionFile moleculesRegressionFile.cpp)
	
target_link_libraries(moleculesRegressionFile
  pandore
  gtb
  GraphKernels
  Molecules
  svm
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesClassificationFile moleculesClassificationFile.cpp)
	
target_link_libraries(moleculesClassificationFile
  pandore
  gtb
  GraphKernels
  Molecules
  svm
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesClassificationSlowAdd moleculesClassificationSlowAdd.cpp)
	
target_link_libraries(moleculesClassificationSlowAdd
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesKPCA moleculesKPCA.cpp)
	
target_link_libraries(moleculesKPCA
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesClassificationKpp moleculesClassificationKpp.cpp)
	
target_link_libraries(moleculesClassificationKpp
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesClassificationKPCA moleculesClassificationKPCA.cpp)
	
target_link_libraries(moleculesClassificationKPCA
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)


add_executable(moleculesTreeletData moleculesTreeletData.cpp)
	
target_link_libraries(moleculesTreeletData
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesComputeGramMatrices moleculesComputeGramMatrices.cpp)
	
target_link_libraries(moleculesComputeGramMatrices
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesComputeGramMatrice moleculesComputeGramMatrice.cpp)
	
target_link_libraries(moleculesComputeGramMatrice
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(testLabeledSpectrum testLabeledSpectrum.cpp)
	
target_link_libraries(testLabeledSpectrum
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)


add_executable(test_string_utils test_string_utils.cpp)
	
target_link_libraries(test_string_utils
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(testTreeletEnumerator testTreeletEnumerator.cpp)
	
target_link_libraries(testTreeletEnumerator
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(splitDataset splitDataset.cpp)
	
target_link_libraries(splitDataset
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesDatasetAnalysis moleculesDatasetAnalysis.cpp)
	
target_link_libraries(moleculesDatasetAnalysis
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(enumerateTreelets enumerateTreelets.cpp)
	
target_link_libraries(enumerateTreelets
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(testEditDistance testEditDistance.cpp)
	
target_link_libraries(testEditDistance
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesVectorFileClassification moleculesVectorFileClassification.cpp)
	
target_link_libraries(moleculesVectorFileClassification
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesVectorFileRegression moleculesVectorFileRegression.cpp)
	
target_link_libraries(moleculesVectorFileRegression
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(testHypergraph testHypergraph.cpp)

target_link_libraries(testHypergraph
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegressionTrainAndTestSet moleculesRegressionTrainAndTestSet.cpp)

target_link_libraries(moleculesRegressionTrainAndTestSet
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesRegressionFileTrainAndTestSet moleculesRegressionFileTrainAndTestSet.cpp)

target_link_libraries(moleculesRegressionFileTrainAndTestSet
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

add_executable(moleculesEditDistance moleculesEditDistance.cpp)

target_link_libraries(moleculesEditDistance
  Molecules
  GraphKernels
  svm
  pandore
  gtb
  ${LINK_CORE_LIBRARIES}
)

# Installation of the executables

INSTALL(TARGETS moleculesRegression DESTINATION bin)
INSTALL(TARGETS moleculesWeightedRegression DESTINATION bin)
INSTALL(TARGETS moleculesRegression10percent DESTINATION bin)
INSTALL(TARGETS moleculesRegression10percent_i DESTINATION bin)
# INSTALL(TARGETS moleculesRegression10percent_tbb DESTINATION bin)
INSTALL(TARGETS moleculesClassification DESTINATION bin)
INSTALL(TARGETS moleculesCalculNbAtoms DESTINATION bin)
INSTALL(TARGETS moleculesCalculNbCycles DESTINATION bin)
INSTALL(TARGETS moleculesCorrelation DESTINATION bin)
INSTALL(TARGETS moleculesRegressionCluster DESTINATION bin)
INSTALL(TARGETS moleculesRegressionFastAdd DESTINATION bin)
INSTALL(TARGETS moleculesRegressionSlowAdd DESTINATION bin)
INSTALL(TARGETS moleculesRegressionSlowAdd2 DESTINATION bin)
INSTALL(TARGETS moleculesRegressionDirectory DESTINATION bin)
INSTALL(TARGETS moleculesRegressionFile DESTINATION bin)
INSTALL(TARGETS moleculesClassificationFile DESTINATION bin)
INSTALL(TARGETS moleculesClassificationSlowAdd DESTINATION bin)
INSTALL(TARGETS moleculesKPCA DESTINATION bin)
INSTALL(TARGETS moleculesClassificationKpp DESTINATION bin)
INSTALL(TARGETS moleculesClassificationKPCA DESTINATION bin)
INSTALL(TARGETS moleculesTreeletData DESTINATION bin)
INSTALL(TARGETS moleculesComputeGramMatrices DESTINATION bin)
INSTALL(TARGETS moleculesComputeGramMatrice DESTINATION bin)
INSTALL(TARGETS splitDataset DESTINATION bin)
INSTALL(TARGETS moleculesDatasetAnalysis DESTINATION bin)
#Tests
INSTALL(TARGETS moleculesClassificationKPCA DESTINATION bin)
INSTALL(TARGETS moleculesTestInversion DESTINATION bin)
INSTALL(TARGETS test_string_utils DESTINATION bin)
INSTALL(TARGETS testTreeletEnumerator DESTINATION bin)
INSTALL(TARGETS enumerateTreelets DESTINATION bin)
INSTALL(TARGETS testEditDistance DESTINATION bin)
INSTALL(TARGETS moleculesVectorFileClassification DESTINATION bin)
INSTALL(TARGETS moleculesVectorFileRegression DESTINATION bin)

INSTALL(TARGETS moleculesRegressionTrainAndTestSet DESTINATION bin)
INSTALL(TARGETS moleculesRegressionFileTrainAndTestSet DESTINATION bin)
INSTALL(TARGETS moleculesEditDistance DESTINATION bin)

