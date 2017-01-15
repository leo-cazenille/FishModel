project(FishModel)
cmake_minimum_required(VERSION 3.2.0)
enable_testing()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(OpenCV 3 REQUIRED)
find_package(OpenMP REQUIRED)

# SYSTEM is used to suppress GCC warnings in libraries headers
include_directories(SYSTEM ${OpenCV_INCLUDE_DIR})

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.4
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
        message(FATAL_ERROR "GCC version must be at least 4.4")
    else()
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
            message("GCC version is later than 5.0")
            set(MORE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 \
            -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual \
            -Wredundant-decls -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default \
            -Wno-unused-parameter -Wold-style-cast -Wfloat-equal -Wundef -Wshadow -Wsuggest-override -Wsuggest-final-methods")
        else()
            message("GCC version is earlier than 5.0")
            set(MORE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 \
            -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual \
            -Wredundant-decls -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default \
            -Wno-unused-parameter -Wold-style-cast -Wfloat-equal -Wundef -Wshadow")
        endif()
    endif()
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#    # require at least clang 3.2
#    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
#        message(FATAL_ERROR "Clang version must be at least 3.2!")
#    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MORE_CXX_FLAGS}")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MORE_CXX_FLAGS} -Wnull-dereference")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MORE_CXX_FLAGS}")
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
endif()


set(srcs
    bm.cpp
    factory.cpp
    model.cpp
    random.cpp
    testModel.cpp
)

add_executable(model ${srcs})
target_link_libraries(model ${OpenCV_LIBS}
                        ${QTGSTREAMER_LIBRARY} ${QTGSTREAMER_LIBRARIES}
                        ${QTGSTREAMER_UTILS_LIBRARY} ${QTGSTREAMER_UTILS_LIBRARIES}
                        ${QTGSTREAMER_UTILS_LIBRARY} ${QTGSTREAMER_UTILS_LIBRARIES}
                        ${QTGLIB_LIBRARIES} ${QTGSTREAMER_UI_LIBRARY} ${QTGSTREAMER_UI_LIBRARIES}
                        ${GSTREAMER_LIBRARIES} ${GSTREAMER_VIDEO_LIBRARIES})
