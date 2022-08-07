//===- CalledValuePropagation.h - Propagate called values -------*- C++ -*-===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file implements a transformation that attaches !callees metadata to
// indirect call sites. For a given call site, the metadata, if present,
// indicates the set of functions the call site could possibly target at
// run-time. This metadata is added to indirect call sites when the set of
// possible targets can be determined by analysis and is known to be small. The
// analysis driving the transformation is similar to constant propagation and
// makes uses of the generic sparse propagation solver.
//
//===----------------------------------------------------------------------===//

#ifndef LLVM_TRANSFORMS_IPO_CALLEDVALUEPROPAGATION_H1
#define LLVM_TRANSFORMS_IPO_CALLEDVALUEPROPAGATION_H1

#include "llvm/IR/PassManager.h"
#include "llvm/Support/CommandLine.h"

namespace llvm {
//extern cl::opt<std::string> CallerName;
//extern cl::opt<std::string> CalleeName;
class CalledValuePropagationPass1
    : public PassInfoMixin<CalledValuePropagationPass1> {
public:
  PreservedAnalyses run(Module &M, ModuleAnalysisManager &);
};
} // namespace llvm

#endif // LLVM_TRANSFORMS_IPO_CALLEDVALUEPROPAGATION_H1
