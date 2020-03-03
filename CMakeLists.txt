cmake_minimum_required(VERSION 3.14)
project(DiningPhilosophers)

set(CMAKE_CXX_STANDARD 17)

add_executable(DiningPhilosophers main.cpp System.h DiningPhilosophers.h readVLTS.h DiningPhilosophersNew.h SyncInT.h DiningPhilosophersEff.h test.h OutputData.h)