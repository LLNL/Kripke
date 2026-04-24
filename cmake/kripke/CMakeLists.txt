cmake_minimum_required(VERSION 3.23)

cmake_policy(SET CMP0057 NEW)
project(KRIPKE LANGUAGES CXX)

if(NOT DEFINED KRIPKE_SOURCE_ROOT)
  message(FATAL_ERROR "KRIPKE_SOURCE_ROOT is not set")
endif()

include("${KRIPKE_SOURCE_ROOT}/Kripke.cmake")
