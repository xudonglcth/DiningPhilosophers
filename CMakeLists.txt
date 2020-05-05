cmake_minimum_required(VERSION 3.14)
project(DiningPhilosophers)

set(CMAKE_CXX_STANDARD 14)
find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(BoostTest main.cpp GroupDining.h IncSyncBuffer.h twoPhils.h Coreach.h diningTest2.h ReachBoolTest.h)
target_link_libraries(BoostTest ${Boost_LIBRARIES})
add_executable(DiningPhilosophers main.cpp System.h
        DiningPhilosophers.h readVLTS.h DiningPhilosophersNew.h SyncInT.h DiningPhilosophersEff.h test.h GVAlgorithmClass.h)
include_directories(/Applications/MATLAB_R2019b.app/extern/include)