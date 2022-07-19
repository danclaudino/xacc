configure_file(pds_vqs.in.cpp ${CMAKE_BINARY_DIR}/quantum/examples/pds_vqs/pds_vqs.cpp)
add_executable(pds_vqs ${CMAKE_BINARY_DIR}/quantum/examples/pds_vqs/pds_vqs.cpp)
target_link_libraries(pds_vqs PRIVATE xacc xacc-quantum-gate)

configure_file(adapt_pds_vqs.in.cpp ${CMAKE_BINARY_DIR}/quantum/examples/pds_vqs/adapt_pds_vqs.cpp)
add_executable(adapt_pds_vqs ${CMAKE_BINARY_DIR}/quantum/examples/pds_vqs/adapt_pds_vqs.cpp)
target_link_libraries(adapt_pds_vqs PRIVATE xacc xacc-quantum-gate)

configure_file(adapt_vqe.in.cpp ${CMAKE_BINARY_DIR}/quantum/examples/pds_vqs/adapt_vqe.cpp)
add_executable(adapt_vqe ${CMAKE_BINARY_DIR}/quantum/examples/pds_vqs/adapt_vqe.cpp)
target_link_libraries(adapt_vqe PRIVATE xacc xacc-quantum-gate)