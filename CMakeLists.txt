cmake_minimum_required (VERSION 3.2)
project (CovarianceTracing CXX)

# For debug build
set(CMAKE_BUILD_TYPE DEBUG)

if ( CMAKE_COMPILER_IS_GNUCC )
   add_definitions ("-Wall -pedantic -fopenmp")
endif ( CMAKE_COMPILER_IS_GNUCC )
add_definitions ("-Og -ffast-math -Wno-deprecated -Wno-deprecated-declarations")

include_directories ("include" "modules" ".")

# Add dependencies
find_package(Threads REQUIRED)
find_package(OpenGL)
find_package(GLUT)

# Add main test suite
add_executable (TestCovariance4D tests/Covariance4D.cpp)
target_compile_features(TestCovariance4D PRIVATE cxx_range_for)

enable_testing()
add_test(TestCovariance4D TestCovariance4D)

add_executable (SmallCov tutorials/tutorial1.cpp)
target_compile_features(SmallCov PRIVATE cxx_range_for)
add_executable (CovFiltering tutorials/tutorial2.cpp)
target_compile_features(CovFiltering PRIVATE cxx_range_for)

if(OPENGL_FOUND AND GLUT_FOUND)
   message("${OPENGL_INCLUDE_DIR}/Headers")
   include_directories("${OPENGL_INCLUDE_DIR}/Headers" ${GLUT_INCLUDE_DIR})
   add_executable(CovFilteringGL tutorials/tutorial2gl.cpp)
   target_compile_features(CovFilteringGL PRIVATE cxx_range_for)
   target_link_libraries(CovFilteringGL ${OPENGL_gl_LIBRARY} ${GLUT_glut_LIBRARY})
endif(OPENGL_FOUND AND GLUT_FOUND)
