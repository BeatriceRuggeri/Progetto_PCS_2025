###########src################



# message(${CMAKE_SOURCE_DIR})
# message(${CMAKE_CURRENT_SOURCE_DIR})
# message(${CMAKE_BINARY_DIR})
# message(${CMAKE_CURRENT_BINARY_DIR})
# message(${PROJECT_NAME})
# message(${CMAKE_PROJECT_NAME})


##########################################################################

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralMesh.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)

##########################################################################

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/supporto_triangolazione.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/triangolazione_I.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/triangolazione_II.cpp)

###########################################################################

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/genesis.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MeshClean.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/funzTopo.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/ExportMesh.cpp)

##########################################################################
##per dopo ##
##########################################################################

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/duale.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/grafo.cpp)

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_headers ${src_headers} PARENT_SCOPE)
set(src_sources ${src_sources} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)