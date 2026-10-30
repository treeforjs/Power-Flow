# CMake generated Testfile for 
# Source directory: /home/treefor/projects/raikou-mfem
# Build directory: /home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test("raikou.mfem_electrostatics" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_mfem_electrostatics")
set_tests_properties("raikou.mfem_electrostatics" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;186;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.mfem_gap_pic" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_mfem_gap_pic")
set_tests_properties("raikou.mfem_gap_pic" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;189;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.mfem_q1_scalar_gather" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_mfem_q1_scalar_gather")
set_tests_properties("raikou.mfem_q1_scalar_gather" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;192;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.core" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_core")
set_tests_properties("raikou.core" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;250;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.c_api" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_c_api")
set_tests_properties("raikou.c_api" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;250;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.backend" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_backend")
set_tests_properties("raikou.backend" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;250;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.distributed" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_distributed")
set_tests_properties("raikou.distributed" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;250;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.adaptive_time" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_adaptive_time")
set_tests_properties("raikou.adaptive_time" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;250;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.distributed_mpi_1" "/usr/sbin/mpiexec" "-n" "1" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_distributed" "--mpi")
set_tests_properties("raikou.distributed_mpi_1" PROPERTIES  TIMEOUT "120" _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;298;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.distributed_mpi_2" "/usr/sbin/mpiexec" "-n" "2" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_distributed" "--mpi")
set_tests_properties("raikou.distributed_mpi_2" PROPERTIES  TIMEOUT "120" _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;298;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.distributed_mpi_4" "/usr/sbin/mpiexec" "-n" "4" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_distributed" "--mpi")
set_tests_properties("raikou.distributed_mpi_4" PROPERTIES  TIMEOUT "120" _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;298;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.mpi" "/usr/sbin/mpiexec" "-n" "2" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou_test_backend" "--mpi")
set_tests_properties("raikou.mpi" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;305;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.cli" "/usr/sbin/python3.14" "/home/treefor/projects/raikou-mfem/tests/test_cli.py" "/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/raikou")
set_tests_properties("raikou.cli" PROPERTIES  _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;368;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
add_test("raikou.python_legacy" "/usr/sbin/python3.14" "/home/treefor/projects/raikou-mfem/tests/test_python.py")
set_tests_properties("raikou.python_legacy" PROPERTIES  ENVIRONMENT "PYTHONPATH=/home/treefor/projects/raikou-mfem/python;RAIKOU_LIBRARY=/home/treefor/projects/runtimes/pic-repo-integration-20261029/mpi-build/libraikou_c.so.1.0.0" _BACKTRACE_TRIPLES "/home/treefor/projects/raikou-mfem/CMakeLists.txt;369;add_test;/home/treefor/projects/raikou-mfem/CMakeLists.txt;0;")
