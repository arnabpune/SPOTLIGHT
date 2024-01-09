cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(spotlight)

find_package(Torch REQUIRED)
find_package(PkgConfig)

pkg_check_modules(PC_ZeroMQ QUIET zmq)

## use the hint from above to find where 'zmq.hpp' is located
find_path(ZeroMQ_INCLUDE_DIR
        NAMES zmq.hpp
        PATHS ${PC_ZeroMQ_INCLUDE_DIRS}
        )

## use the hint from above to find the location of libzmq
find_library(ZeroMQ_LIBRARY
        NAMES zmq
        PATHS ${PC_ZeroMQ_LIBRARY_DIRS}
        )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -I . -g -lzmq -fopenmp -O3")

add_executable(production_noRL production_final.cpp)
target_include_directories(production_noRL  PUBLIC ${ZeroMQ_INCLUDE_DIR})
target_link_libraries(production_noRL  PUBLIC ${ZeroMQ_LIBRARY})
set_property(TARGET production_noRL  PROPERTY CXX_STANDARD 14)

add_executable(prodvacPT prodvac_pytorch.cpp)
target_include_directories(prodvacPT PUBLIC ${ZeroMQ_INCLUDE_DIR})
target_link_libraries(prodvacPT PUBLIC ${ZeroMQ_LIBRARY} ${TORCH_LIBRARIES})
set_property(TARGET prodvacPT PROPERTY CXX_STANDARD 14)

add_executable(parprodvacPT prodvac_pytorch_parallelgen.cpp)
target_include_directories(parprodvacPT PUBLIC ${ZeroMQ_INCLUDE_DIR})
target_link_libraries(parprodvacPT PUBLIC ${ZeroMQ_LIBRARY} ${TORCH_LIBRARIES})
set_property(TARGET parprodvacPT PROPERTY CXX_STANDARD 14)

add_executable(productionPT production_pytorch_prot.cpp)
target_include_directories(productionPT PUBLIC ${ZeroMQ_INCLUDE_DIR})
target_link_libraries(productionPT PUBLIC ${ZeroMQ_LIBRARY} ${TORCH_LIBRARIES})
set_property(TARGET productionPT PROPERTY CXX_STANDARD 14)

add_executable(prodprunePT prodvac_pytorch_parprunegen.cpp)
target_include_directories(prodprunePT PUBLIC ${ZeroMQ_INCLUDE_DIR})
target_link_libraries(prodprunePT PUBLIC ${ZeroMQ_LIBRARY} ${TORCH_LIBRARIES})
set_property(TARGET prodprunePT PROPERTY CXX_STANDARD 14)

add_executable(productionprunePT prodprune_pytorch_prot.cpp)
target_include_directories(productionprunePT PUBLIC ${ZeroMQ_INCLUDE_DIR})
target_link_libraries(productionprunePT PUBLIC ${ZeroMQ_LIBRARY} ${TORCH_LIBRARIES})
set_property(TARGET productionprunePT PROPERTY CXX_STANDARD 14)

add_executable(prodvacPT_trainable prodvac_pytorch_trainable.cpp)
target_include_directories(prodvacPT_trainable  PUBLIC ${ZeroMQ_INCLUDE_DIR})
target_link_libraries(prodvacPT_trainable  PUBLIC ${ZeroMQ_LIBRARY} ${TORCH_LIBRARIES})
set_property(TARGET prodvacPT_trainable  PROPERTY CXX_STANDARD 14)

include_directories(BEFORE ../dnv)

# The following code block is suggested to be used on Windows.
# According to https://github.com/pytorch/pytorch/issues/25457,
# the DLLs need to be copied to avoid memory errors.
if (MSVC)
  file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
  add_custom_command(TARGET example-app
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${TORCH_DLLS}
                     $<TARGET_FILE_DIR:app1>)
endif (MSVC)

