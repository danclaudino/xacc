include_directories(${CMAKE_BINARY_DIR})
add_xacc_test(QWC)
target_link_libraries(QWCTester xacc xacc-quantum-gate)