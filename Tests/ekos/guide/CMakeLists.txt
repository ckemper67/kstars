SET( StarCorrespondenceTests_SRCS teststarcorrespondence.cpp  )

IF (INDI_FOUND)
  INCLUDE_DIRECTORIES(${INDI_INCLUDE_DIR})
  ADD_EXECUTABLE( testguidestars testguidestars.cpp )
  TARGET_LINK_LIBRARIES( testguidestars ${TEST_LIBRARIES})
  ADD_TEST( NAME GuideStarsTest COMMAND testguidestars )
  SET_TESTS_PROPERTIES( GuideStarsTest PROPERTIES LABELS "stable")
ENDIF ()

ADD_EXECUTABLE( teststarcorrespondence teststarcorrespondence.cpp )
TARGET_LINK_LIBRARIES( teststarcorrespondence ${TEST_LIBRARIES})
ADD_TEST( NAME StarCorrespondenceTest COMMAND teststarcorrespondence )
SET_TESTS_PROPERTIES( StarCorrespondenceTest PROPERTIES LABELS "stable")

ADD_EXECUTABLE( testcalibrationprocess testcalibrationprocess.cpp )
TARGET_LINK_LIBRARIES( testcalibrationprocess ${TEST_LIBRARIES})
ADD_TEST( NAME CalibrationProcessTest COMMAND testcalibrationprocess )
SET_TESTS_PROPERTIES( CalibrationProcessTest PROPERTIES LABELS "stable")

FIND_PACKAGE(Eigen3 REQUIRED)
ADD_EXECUTABLE( testguidealgorithms testguidealgorithms.cpp )
TARGET_INCLUDE_DIRECTORIES( testguidealgorithms PRIVATE
    ${CMAKE_BINARY_DIR}/kstars
    ${EIGEN3_INCLUDE_DIR}
)
TARGET_LINK_LIBRARIES( testguidealgorithms ${TEST_LIBRARIES})
ADD_TEST( NAME GuideAlgorithmsTest COMMAND testguidealgorithms )
SET_TESTS_PROPERTIES( GuideAlgorithmsTest PROPERTIES LABELS "stable")

ADD_EXECUTABLE( benchmarkguidealgorithms benchmarkguidealgorithms.cpp )
TARGET_INCLUDE_DIRECTORIES( benchmarkguidealgorithms PRIVATE
    ${CMAKE_BINARY_DIR}/kstars
    ${EIGEN3_INCLUDE_DIR}
)
TARGET_LINK_LIBRARIES( benchmarkguidealgorithms ${TEST_LIBRARIES})
