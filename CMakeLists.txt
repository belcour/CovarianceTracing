cmake_minimum_required (VERSION 3.2)
project (CovarianceTracing CXX)

# For debug build
set(CMAKE_BUILD_TYPE DEBUG)

if ( CMAKE_COMPILER_IS_GNUCC )
   add_definitions ("-Wall -pedantic -fopenmp")
endif ( CMAKE_COMPILER_IS_GNUCC )
add_definitions ("-O3 -ffast-math")

include_directories ("include" "modules" ".")

# Add dependencies
find_package(Threads REQUIRED)

# Add main test suite
add_executable (TestCovariance4D tests/Covariance4D.cpp)
target_compile_features(TestCovariance4D PRIVATE cxx_range_for)

enable_testing()
add_test(TestCovariance4D TestCovariance4D)

add_executable (SmallCov tutorials/tutorial1.cpp)
target_compile_features(SmallCov PRIVATE cxx_range_for)
add_executable (CovFiltering tutorials/tutorial2.cpp)
target_compile_features(CovFiltering PRIVATE cxx_range_for)
