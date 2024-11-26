ADD_SUBDIRECTORY(External)

ADD_COMPILE_DEFINITIONS(DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

LINK_LIBRARIES(
    tinynurbs::tinynurbs
    doctest::doctest
)

ADD_EXECUTABLE(TestFindSpan TestFindSpan.cpp)
ADD_EXECUTABLE(TestBSplineBasis TestBSplineBasis.cpp)
ADD_EXECUTABLE(TestBSplineDerBasis TestBSplineDerBasis.cpp)
ADD_EXECUTABLE(TestCurvePoint TestCurvePoint.cpp)
ADD_EXECUTABLE(TestCurveDerivatives TestCurveDerivatives.cpp)
ADD_EXECUTABLE(TestSurfacePoint TestSurfacePoint.cpp)
ADD_EXECUTABLE(TestSurfaceDerivatives TestSurfaceDerivatives.cpp)
ADD_EXECUTABLE(TestSurfaceNormal TestSurfaceNormal.cpp)
