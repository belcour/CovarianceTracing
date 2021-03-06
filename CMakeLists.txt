cmake_minimum_required (VERSION 3.2)
project (CovarianceTracing CXX)

# For debug build
set(CMAKE_BUILD_TYPE RELEASE)

if ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
   add_definitions ("-Wall -pedantic  -Wno-deprecated -Wno-deprecated-declarations")
   add_definitions ("-Ofast -ffast-math -fopenmp")
endif ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

include_directories ("include" "modules" ".")

# Add dependencies
find_package(OpenGL)
find_package(GLUT)
find_package(OpenMP)

# Add main test suite
add_executable (TestCovariance4D    tests/Covariance4D.cpp)
add_executable (TestInvCovariance4D tests/InvCovariance4D.cpp)
target_compile_features(TestCovariance4D    PRIVATE cxx_range_for)
target_compile_features(TestInvCovariance4D PRIVATE cxx_range_for)

enable_testing()
add_test(TestCovariance4D    TestCovariance4D)
add_test(TestInvCovariance4D TestInvCovariance4D)

add_executable (Tutorial1 tutorials/tutorial1.cpp)
target_compile_features(Tutorial1 PRIVATE cxx_range_for)
add_executable (Tutorial2 tutorials/tutorial2.cpp)
target_compile_features(Tutorial2 PRIVATE cxx_range_for)

if(OPENGL_FOUND AND GLUT_FOUND)
   add_definitions("-DGL_GLEXT_PROTOTYPES")
   include_directories("${OPENGL_INCLUDE_DIR}/Headers" ${GLUT_INCLUDE_DIR})
   add_executable(Tutorial2gl tutorials/tutorial2gl.cpp)
   target_compile_features(Tutorial2gl PRIVATE cxx_range_for)
   target_link_libraries(Tutorial2gl ${OPENGL_gl_LIBRARY} ${GLUT_glut_LIBRARY})
endif(OPENGL_FOUND AND GLUT_FOUND)

# Add OpenMP support
if ( OPENMP_FOUND )
	target_compile_options(Tutorial1   PUBLIC ${OpenMP_CXX_FLAGS})
	target_compile_options(Tutorial2   PUBLIC ${OpenMP_CXX_FLAGS})
	target_compile_options(Tutorial2gl PUBLIC ${OpenMP_CXX_FLAGS})
	if (NOT WIN32 )
		set_target_properties(Tutorial1   PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
		set_target_properties(Tutorial2   PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
		set_target_properties(Tutorial2gl PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif (NOT WIN32 )
endif( OPENMP_FOUND )