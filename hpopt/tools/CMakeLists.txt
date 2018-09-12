
#-----------------------------------------------------------#

include_directories(${PROJECT_SOURCE_DIR}/include)

#-----------------------------------------------------------#

set(LIBRARIES hpopt m)

#-----------------------------------------------------------#

add_executable(workdist workdist.cpp)
add_executable(func_generator func_generator.cpp)
add_executable(compDist compDist.cpp)
#-----------------------------------------------------------#

target_link_libraries(workdist ${LIBRARIES})
target_link_libraries(compDist ${LIBRARIES})

#-----------------------------------------------------------#

install(TARGETS workdist
                func_generator
                compDist
                DESTINATION tools)

#-----------------------------------------------------------#
