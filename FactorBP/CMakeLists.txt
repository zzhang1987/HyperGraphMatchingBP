Project(HyperGMBP)
cmake_minimum_required(VERSION 3.4)


macro(ADD_USER_FRAMEWORK fwname appname)
    find_library(FRAMEWORK_${fwname}
    NAMES ${fwname}
    PATHS /Library/Frameworks
    NO_DEFAULT_PATH
    REQUIRED)
    MARK_AS_ADVANCED(${fwname})
    if( ${FRAMEWORK_${fwname}} STREQUAL FRAMEWORK_${fwname}-NOTFOUND)
        MESSAGE(ERROR ": Framework ${fwname} not found")
    else()
        #INCLUDE_DIRECTORIES(${fwname}_INCLUDE_DIR)
	#MESSAGE(${fwname}_INCLUDE_DIR)
        TARGET_LINK_LIBRARIES(${appname} "${FRAMEWORK_${fwname}}/${fwname}")
        MESSAGE(STATUS "Framework ${fwname} found at ${FRAMEWORK_${fwname}}")
    endif()
endmacro(ADD_USER_FRAMEWORK)




set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -std=c++1y") 

SET(LIBSRCS Auction.cpp
	    FactorCreator.cpp
	    Factor.cpp
	    FactorGraph.cpp
	    SubTourFactor.cpp
	    Auction.h
	    FactorCreator.h
	    FactorStore.h
	    BaBTypes.h
	    FactorGraph.h
	    PRTypes.h
	    Factor.h
	    FactorGraphStore.h
	    SubTourFactor.h
	    )

set (HEADERS  Auction.h
	    FactorCreator.h
	    FactorStore.h
	    BaBTypes.h
	    FactorGraph.h
	    PRTypes.h
	    Factor.h
	    FactorGraphStore.h
	    SubTourFactor.h)


add_library(GMBP STATIC ${LIBSRCS})
ADD_USER_FRAMEWORK(boost GMBP)
add_executable(testMain testMain.cpp)

target_link_libraries(testMain ${EXTRA_LIBS} GMBP)
ADD_USER_FRAMEWORK(boost testMain)
