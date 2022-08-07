# If we don't need RTTI or EH, there's no reason to export anything
# from the hello plugin.


add_llvm_component_library(LLVMHello
  Hello.cpp

  ADDITIONAL_HEADER_DIRS
  ${LLVM_MAIN_INCLUDE_DIR}/llvm/Transforms

  DEPENDS
  intrinsics_gen
  )
