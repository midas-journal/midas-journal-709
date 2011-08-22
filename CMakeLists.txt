cmake_minimum_required(VERSION 2.6)

if(COMMAND cmake_policy)
cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

PROJECT(vtkRANSACPlane)
ENABLE_TESTING()

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

############
ADD_EXECUTABLE(vtkRANSACPlaneExample vtkRANSACPlaneExample.cxx
vtkRANSACPlane.cxx)
TARGET_LINK_LIBRARIES(vtkRANSACPlaneExample vtkHybrid )

###########
ADD_EXECUTABLE(vtkRANSACPlaneTest vtkRANSACPlaneTest.cxx
vtkRANSACPlane.cxx)
TARGET_LINK_LIBRARIES(vtkRANSACPlaneTest vtkHybrid )
ADD_TEST(vtkRANSACPlaneTest vtkRANSACPlaneTest)

