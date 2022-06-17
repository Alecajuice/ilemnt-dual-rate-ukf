cmake_minimum_required(VERSION 3.0.0)
project(ilemt-dual-rate-ukf VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(ilemt-dual-rate-ukf dual_rate_ukf.cpp low_rate_ukf.cpp high_rate_ukf.cpp)

include_directories(include)

# target_link_libraries (ilemt-dual-rate-ukf "${CMAKE_CURRENT_SOURCE_DIR}/include/libz.a")
# target_link_libraries (ilemt-dual-rate-ukf "${CMAKE_CURRENT_SOURCE_DIR}/include/libmatio.a")
target_link_libraries (ilemt-dual-rate-ukf z)
target_link_libraries (ilemt-dual-rate-ukf matio) # (built with zlib)
# target_link_libraries (ilemt-dual-rate-ukf "/usr/local/lib/libz.a")
# target_link_libraries (ilemt-dual-rate-ukf "/usr/local/lib/libmatio.a")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
