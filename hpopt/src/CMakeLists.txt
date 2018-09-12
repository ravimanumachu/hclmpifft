#-----------------------------------------------------------#

set (HPOPTHEADERS 
             ${PROJECT_SOURCE_DIR}/include/dirhelper.h
             ${PROJECT_SOURCE_DIR}/include/filehelper.h
             ${PROJECT_SOURCE_DIR}/include/memorization.h
             ${PROJECT_SOURCE_DIR}/include/wdistmethods.h
             ${PROJECT_SOURCE_DIR}/include/functionclass.h
             ${PROJECT_SOURCE_DIR}/include/hcl_hpopta.h)
             
#-----------------------------------------------------------#

set (HPOPTSOURCES 
            dirhelper.cpp
            filehelper.cpp
            memorization.cpp
            wdistmethods.cpp
            functionclass.cpp
            hcl_hpopta.cpp)
             
#-----------------------------------------------------------#

add_library(hpopt
            ${HPOPTHEADERS} 
            ${HPOPTSOURCES})

target_link_libraries(hpopt)

#-----------------------------------------------------------#

install(TARGETS hpopt DESTINATION lib)
install(FILES ${HPOPTHEADERS} DESTINATION include)

#-----------------------------------------------------------#
