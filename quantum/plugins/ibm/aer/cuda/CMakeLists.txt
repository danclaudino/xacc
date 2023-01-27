set(LIBRARY_NAME aer_custatevec)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)

enable_language(CUDA)
file(GLOB SRC aer_custatevec.cu)

set(CUDA_ARCH_BIN "70")
set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -shared -forward-unknown-to-host-compiler -std=c++14 -gencode arch=compute_${CUDA_ARCH_BIN},code=sm_${CUDA_ARCH_BIN} -rdc=true --compiler-options -fPIC ${CMAKE_CXX_FLAGS}")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -forward-unknown-to-host-compiler -std=c++14 -gencode arch=compute_${CUDA_ARCH_BIN},code=sm_${CUDA_ARCH_BIN} -rdc=true --compiler-options -fPIC ${CMAKE_CXX_FLAGS}")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -shared -std=c++14 -gencode arch=compute_${CUDA_ARCH_BIN},code=sm_${CUDA_ARCH_BIN} -rdc=true --compiler-options -fPIC")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DAER_THRUST_CUDA -isystem -use_fast_math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAER_THRUST_CUDA -DAER_THRUST_SUPPORTED")

#cuda_add_library(${LIBRARY_NAME} SHARED ${SRC})
add_library(${LIBRARY_NAME} STATIC ${SRC})
#add_library(${LIBRARY_NAME} SHARED ${SRC})

set_target_properties(${LIBRARY_NAME}
                      PROPERTIES
                      CUDA_SEPARABLE_COMPILATION ON
                      CUDA_RESOLVE_DEVICE_SYMBOLS ON)

set_target_properties(${LIBRARY_NAME} PROPERTIES DEFINE_SYMBOL "")

target_include_directories(${LIBRARY_NAME} PRIVATE 
                                            $ENV{CUQUANTUM_INSTALL_PATH}/include
                                            ../src
                                            #${AER_SOURCE_PATH}/src
                                            #${CMAKE_SOURCE_DIR}/tpls/spdlog
                                            ${CMAKE_SOURCE_DIR}/tpls
                                            ${CMAKE_SOURCE_DIR}/tpls/spd/include
  )
set(AER_COMPILER_DEFINITIONS ${AER_COMPILER_DEFINITIONS} THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA)
set(AER_COMPILER_DEFINITIONS ${AER_COMPILER_DEFINITIONS} AER_CUSTATEVEC)
set(AER_COMPILER_DEFINITIONS ${AER_COMPILER_DEFINITIONS} AER_THRUST_SUPPORTED=TRUE)

target_link_libraries(${LIBRARY_NAME} PRIVATE $ENV{CUQUANTUM_INSTALL_PATH}/lib/libcustatevec.so)
target_compile_definitions(${LIBRARY_NAME} PRIVATE ${AER_COMPILER_DEFINITIONS})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${LIBRARY_NAME} PUBLIC OpenMP::OpenMP_CXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

install(TARGETS ${LIBRARY_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
