cmake_minimum_required(VERSION 3.5)

# Set the project name
project (SmartCFS)


# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    SmartSFC.cpp
    stdafx.cpp
	stdafx.h
	targetver.h
)

# Create a sources variable with a link to all cpp files to compile
set(CGAL_SOURCES
	cgal/cgalpackage.cpp
	cgal/cgalpackage.h
	cgal/clipper.cpp
	cgal/clipper.hpp
	cgal/dllmain.cpp
	cgal/IntegralCurvature.cpp
	cgal/IO.cpp
	cgal/kdtree.cpp
	cgal/kdtree.h
	cgal/NewtonApple_hull3D.cpp
	cgal/NewtonApple_hull3D.h
)

set(PYMESH_SOURCES
	pymesh/Exception.h
	pymesh/MshLoader.cpp
	pymesh/MshLoader.h
	pymesh/MshSaver.cpp
	pymesh/MshSaver.h
)

SOURCE_GROUP("cgal" FILES ${CGAL_SOURCES})

SOURCE_GROUP("pymesh" FILES ${PYMESH_SOURCES})


set(CFS_SOURCES
	cfs/Base.cpp
	cfs/CFSCNC.cpp
	cfs/CFSCNC.h
	cfs/Circuit.h
	cfs/Constant.h
	cfs/MathHelper.h
	cfs/Output.cpp
	cfs/RoughMachining.cpp
	cfs/RoughMachining.h
	cfs/scallop_path.cpp
	cfs/Source.cpp
	cfs/Strip.h
	cfs/ToolPathTimeEstimator.cpp
	cfs/ToolPathTimeEstimator.hpp
	cfs/Tree.h
	cfs/Vector2D.cpp
	cfs/Vector2D.hpp
	cfs/AM.cpp
	cfs/AM.h
)

SOURCE_GROUP("cfs" FILES ${CFS_SOURCES})


# Add an executable with the above sources
add_executable(SmartCFS ${SOURCES} ${CGAL_SOURCES} ${PYMESH_SOURCES} ${CFS_SOURCES})

set(CFS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/cfs)
set(CGAL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/cgal)
set(PYMESH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/pymesh)
set(ThirdPartyLibs_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/include)



include(ExternalProject)
ExternalProject_Add(
    pgl
    PREFIX ${CMAKE_BINARY_DIR}/pgl
    GIT_REPOSITORY https://github.com/haisenzhao/personal-geom-lib.git
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    )
ExternalProject_Get_Property(pgl source_dir)

set(PglIncludeDir ${source_dir})
set(PglGlmIncludeDir ${source_dir}/glm)


#MESSAGE(STATUS ${CGAL_INCLUDE_DIR} " found")
include_directories(${PROJECT_SOURCE_DIR} ${CGAL_INCLUDE_DIR} ${ThirdPartyLibs_INCLUDE_DIR} ${PYMESH_INCLUDE_DIR} ${CFS_INCLUDE_DIR} ${PglIncludeDir} ${PglGlmIncludeDir})

TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/CGAL-vc120-mt-4.7.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/CGAL-vc120-mt-gd-4.7.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/CGAL_Core-vc120-mt-4.7.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/CGAL_Core-vc120-mt-gd-4.7.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/CGAL_ImageIO-vc120-mt-4.7.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/CGAL_ImageIO-vc120-mt-gd-4.7.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/libboost_chrono-vc120-mt-1_59.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/libboost_date_time-vc120-mt-1_59.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/libboost_system-vc120-mt-1_59.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/libboost_thread-vc120-mt-1_59.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/libgmp-10.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/libmpfr-4.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/Wm5Core.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/Wm5Imagics.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/Wm5Mathematics.lib)
TARGET_LINK_LIBRARIES(SmartCFS ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib/Wm5WglGraphics.lib)
	
	
add_dependencies(${PROJECT_NAME} pgl)


#find_library(LIBHELLO_PATH CGAL_Core-vc120-mt-4.7 ${PROJECT_SOURCE_DIR}/ThirdPartyLibs/lib)
#IF(NOT LIBHELLO_PATH)
#MESSAGE(FATAL_ERROR "libhello not found")
#ENDIF(NOT LIBHELLO_PATH) 
#
#MESSAGE(STATUS ${LIBHELLO_PATH} " found")
#TARGET_LINK_LIBRARIES(SmartCFS ${LIBHELLO_PATH})


