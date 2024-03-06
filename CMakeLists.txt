cmake_minimum_required(VERSION 3.17)
project(publicMatrix)

set(CMAKE_CXX_STANDARD 14)

add_executable(test md2.c lhe.c lhep.c bn_ext.c pi2.c test.c)

target_link_libraries(test libgmp.so libflint.so librelic.so librelic_s.so)