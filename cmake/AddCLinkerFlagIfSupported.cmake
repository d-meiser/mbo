include(CheckCCompilerFlag)
macro(add_c_linker_flag_if_supported FLAG)
  CHECK_C_COMPILER_FLAG(${FLAG} FLAG_SUPPORTED)
  if(FLAG_SUPPORTED)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${FLAG}")
  endif()
endmacro(add_c_linker_flag_if_supported)


