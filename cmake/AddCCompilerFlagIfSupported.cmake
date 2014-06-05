include(CheckCCompilerFlag)
macro(add_c_compiler_flag_if_supported FLAG)
  CHECK_C_COMPILER_FLAG(${FLAG} FLAG_SUPPORTED)
  if(FLAG_SUPPORTED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FLAG}")
  endif()
endmacro(add_c_compiler_flag_if_supported)

