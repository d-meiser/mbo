include(CheckCCompilerFlag)
function(add_c_compiler_flag_if_supported FLAG FLAGS)
  CHECK_C_COMPILER_FLAG(${FLAG} FLAG_SUPPORTED)
  if(FLAG_SUPPORTED)
    set(FLAGS "${FLAGS} ${FLAG}")
  endif()
endfunction(add_c_compiler_flag_if_supported)

