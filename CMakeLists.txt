cmake_minimum_required(VERSION 3.8)
project(NRPsMatcher)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp
        Aminoacid/Aminoacids.cpp
        Aminoacid/Aminoacids.h
        NRPsPrediction/AminoacidPrediction.cpp
        NRPsPrediction/AminoacidPrediction.h
        NRPsPrediction/NRPsPart.cpp
        NRPsPrediction/NRPsPart.h
        NRPsPrediction/NRPsPrediction.cpp
        NRPsPrediction/NRPsPrediction.h)
add_executable(NRPsMatcher ${SOURCE_FILES})