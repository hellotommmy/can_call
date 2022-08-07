//===- Hello.cpp - Example code from "Writing an LLVM Pass" ---------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file implements two versions of the LLVM "Hello World" pass described
// in docs/WritingAnLLVMPass.html
//
//===----------------------------------------------------------------------===//

#include "llvm/Transforms/Hello/CalledValuePropagation1.h"
#include "llvm/Analysis/SparsePropagation.h"
#include "llvm/Analysis/ValueLatticeUtils.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/MDBuilder.h"
#include "llvm/InitializePasses.h"
#include "llvm/Pass.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Transforms/IPO.h"
#include"./GraphDataStruct.h"

using namespace llvm;
#include<iostream>
using std::cout;
#define DEBUG_TYPE "called-value-propagation1"




vector<Function *> trackCandidateCallees(CallBase* callI){
  vector<Function *> candidates;
  errs() << "in trackCandidate!\n";
  Instruction * parentI= dyn_cast<Instruction >(callI->getCalledOperand());
  int trackPathLen = 0;
  while(trackPathLen < 100 && parentI){
    switch (parentI->getOpcode()) {
    //case Instruction::Call:
    //case Instruction::Invoke:
    //  return visitCallBase(cast<CallBase>(I), ChangedValues, SS);
    case Instruction::Load:
      {
      LoadInst * parentIL = dyn_cast<LoadInst>(parentI);
      parentI = dyn_cast<Instruction>(parentIL->getPointerOperand());
      break;
      }
    case Instruction::Store:
      {

        errs() << "Did it ever get here?\n";
      StoreInst * parentIS = dyn_cast<StoreInst >(parentI);
      Value * v = parentIS -> getValueOperand();
      Type * t = v->getType();
      if(t->isFunctionTy()){
        candidates.push_back(dyn_cast<Instruction>(v)->getFunction());
        return candidates;
      }
      else{
        return candidates;
      }
      break;
      }
    default:
      return candidates; 
    }
    trackPathLen++;
  }
  return candidates;
}



static bool runCVP1(Module &M) {

  DiGraph callGraph; 
  for(auto & F: M){
      llvm::StringRef callerName = F.getName();
      errs() << "functions in file:\n";
      errs().write_escaped(F.getName()) << '\n';
      callGraph.addVertex(callerName.str());
      for (auto& B : F) {
        for (auto& I : B) {
          if (auto* callI = dyn_cast<CallBase>(&I)) {
            // Insert at the point where the instruction `op` appears.


            Function * callee = callI->getCalledFunction();
            llvm::StringRef calleeName;
            if(callee == nullptr){
//              errs() << "call with func ptr\n";
//              errs() << callI->getType();
              auto callees = trackCandidateCallees(callI);
              for(auto callee : callees){
                calleeName = callee->getName();
                callGraph.addVertex(calleeName.str());
              }
            }
            else{
              calleeName = callee->getName();
              callGraph.addVertex(calleeName.str());
              callGraph.addEdgeTo(callerName.str(), calleeName.str());

            }
          }
        }
      }
  }

      //errs() << *callGraph[0] << "\n";

      //errs() << "Hello!!!!!\n";

  //string callerInput(CallerName.c_str());
  //string calleeInput(CalleeName.c_str());
  errs() << "Outputing all possibilities of caller--->callee:\n";
  for(auto & v: callGraph.vertices) {
    for(auto &  w: callGraph.vertices) {
      if(callGraph.isPathDFS(v->label, w->label)){
          errs() << *v << "--->" ;
          errs() << *w;
          errs() << "\n";
          errs() << "true\n";
      }
      else{
          errs() << *v << "--->" ;
          errs() << *w;
          errs() << "\n";
          errs() << "false\n";
      }
    }
  }

     // errs() << callGraph;
     // if(callGraph.isPathDFS("main", "bizz"))
     //   errs() << "There is a path from main to bizz\n";
     // else
     //   errs() << "no path from main to bizz\n";

  return false;
}

PreservedAnalyses CalledValuePropagationPass1::run(Module &M,
                                                  ModuleAnalysisManager &) {
  runCVP1(M);
  return PreservedAnalyses::all();
}

namespace {
class CalledValuePropagationLegacyPass1 : public ModulePass {
public:
  static char ID;

  void getAnalysisUsage(AnalysisUsage &AU) const override {
    AU.setPreservesAll();
  }

  CalledValuePropagationLegacyPass1() : ModulePass(ID) {
    initializeCalledValuePropagationLegacyPass1Pass(
        *PassRegistry::getPassRegistry());
  }

  bool runOnModule(Module &M) override {
    if (skipModule(M))
    {
      return false;
    }
    return runCVP1(M);
  }
};
} // namespace

char CalledValuePropagationLegacyPass1::ID = 0;
INITIALIZE_PASS(CalledValuePropagationLegacyPass1, "called-value-propagation1",
                "Called Value Propagation1", false, false)

ModulePass *llvm::createCalledValuePropagationPass1() {
  return new CalledValuePropagationLegacyPass1();
}

//#include "llvm/ADT/Statistic.h"
//#include "llvm/IR/Function.h"
//#include "llvm/Pass.h"
//#include "llvm/Support/raw_ostream.h"
//
//#include "llvm/IR/LegacyPassManager.h"
//#include "llvm/IR/InstrTypes.h"
//#include "llvm/IR/IRBuilder.h"
//#include "llvm/Transforms/Utils/BasicBlockUtils.h"
//#include "llvm/Transforms/IPO/PassManagerBuilder.h"
//
//
//
//#include "llvm/Transforms/IPO/Attributor.h"
//#include "llvm/Analysis/AliasAnalysis.h"
//#include "llvm/Analysis/CallGraph.h"
//#include "llvm/Analysis/CallGraphSCCPass.h"
//
//using namespace llvm;
//
//#include"./GraphDataStruct.h"
//
//#define DEBUG_TYPE "hello"
//
//STATISTIC(HelloCounter, "Counts number of functions greeted");
//
//namespace {
//  // Hello - The first implementation, without getAnalysisUsage.
//  struct Hello : public FunctionPass {
//    static char ID; // Pass identification, replacement for typeid
//    Hello() : FunctionPass(ID) {}
//
//    bool runOnFunction(Function &F) override {
//      ++HelloCounter;
//      errs() << "Hello: ";
//      errs().write_escaped(F.getName()) << '\n';
//      return false;
//    }
//  };
//}
//
//char Hello::ID = 0;
//static RegisterPass<Hello> X("hello", "Hello World Pass");
//
//namespace {
//  // Hello2 - The second implementation with getAnalysisUsage implemented.
//  struct Hello2 : public FunctionPass {
//    static char ID; // Pass identification, replacement for typeid
//    Hello2() : FunctionPass(ID) {}
//
//    bool runOnFunction(Function &F) override {
//
//#include"./GraphDataStruct.h"
//      llvm::StringRef callerName = F.getName();
//      if(callerName.str() == "baz")
//        errs() << F;
//      errs().write_escaped(F.getName()) << '\n';
//      DiGraph callGraph; 
//      callGraph.addVertex(callerName.str());
//      for (auto& B : F) {
//        for (auto& I : B) {
//          if (auto* callI = dyn_cast<CallBase>(&I)) {
//            // Insert at the point where the instruction `op` appears.
//            IRBuilder<> builder(callI);
//
//
//            Function * callee = callI->getCalledFunction();
//            if(callee == nullptr){
//              errs() << "call with func ptr\n";
//              errs() << callI->getType();
//              errs() << "\n intrinsicID is " << callI->getIntrinsicID();
//            }
//            else{
//              llvm::StringRef calleeName = callee->getName();
//              callGraph.addVertex(calleeName.str());
//              callGraph.addEdgeTo(callerName.str(), calleeName.str());
//
//              if(callee)
//                errs() << "callee is " <<  callee->getName() << "\n";
//              //callGraph.insert(callee->getName)
//              else{
//                errs() << "cannot extract callee\n";
//              }
//            }
//          }
//        }
//      }
//
//      //errs() << *callGraph[0] << "\n";
//
//      //errs() << "Hello!!!!!\n";
//
//      errs() << callGraph;
//      if(callGraph.isPathDFS("main", "bizz"))
//        errs() << "There is a path from main to bizz\n";
//      else
//        errs() << "no path from main to bizz\n";
//
//      return false;
//    }
//  };
//
//
//
//}
//
//char Hello2::ID = 0;
//static RegisterPass<Hello2>
//Y("hello2", "Hello World Pass (with getAnalysisUsage implemented)");
//
//namespace {
//
//struct Hello3: public CallGraphSCCPass {
//  static char ID;
//
//  Hello3() : CallGraphSCCPass(ID) {
//    //initializeAttributorCGSCCLegacyPassPass(*PassRegistry::getPassRegistry());
//  }
//
//  bool runOnSCC(CallGraphSCC &SCC) override {
//    if (skipSCC(SCC))
//    {
//      errs() << "SCC pass skipped hello3";
//      return false;
//    }
//
//    errs() << "*****************************dumping graph started\n";
//    SCC.getCallGraph().dump();
//    errs() << "*****************************dumping graph ended\n";
//    errs() << "**************************************************************************************************\n";
//    SCC.getCallGraph().getExternalCallingNode()->dump(); 
//    errs() << "**************************************************************************************************\n";
//    SetVector<Function *> Functions;
//    errs() << "*************************************************\n";
//    errs() << "***********************dumping nodes started\n";
//    for(auto const node:  SCC){
//      errs() << "\tnode:\n";
//      node->dump();
//    }
//    errs() << "***********************dumping nodes ended\n";
//    for (CallGraphNode *CGN : SCC)
//    {
//      //errs() << "node: " << CGN << "\n";
//      if (Function *Fn = CGN->getFunction())
//      {
//        //errs() << "Function corresponding to node in call graph "  << " :\n " << *Fn << "\n";
//        for(const Use &U : Fn->uses()){
//          if(const auto * CB = dyn_cast<CallBase>(U.getUser()))
//          {
//            //errs() << "callbase: " << *CB << "\n";
//          }
//          else 
//          {
//            //errs() << U << "\n";
//            //errs() << "indirect calls?\n";
//          }
//        }
//        if (!Fn->isDeclaration()){
//          //errs() << "this is not a function declaration\n";
//          Functions.insert(Fn);
//        }
//      }
//
//    }
//
//    if (Functions.empty())
//      return false;
//
//    AnalysisGetter AG;
//    CallGraph &CG = const_cast<CallGraph &>(SCC.getCallGraph());
//    CallGraphUpdater CGUpdater;
//    CGUpdater.initialize(CG, SCC);
//    Module &M = *Functions.back()->getParent();
//    BumpPtrAllocator Allocator;
//    InformationCache InfoCache(M, AG, Allocator, /* CGSCC */ &Functions);
//    //return runAttributorOnFunctions(InfoCache, Functions, AG, CGUpdater,
//    //                                /* DeleteFns */ false,
//    //                                /* IsModulePass */ false);
//  }
//
//  void getAnalysisUsage(AnalysisUsage &AU) const override {
//    // FIXME: Think about passes we will preserve and add them here.
//    AU.addRequired<TargetLibraryInfoWrapperPass>();
//    CallGraphSCCPass::getAnalysisUsage(AU);
//  }
//};
//}
//char Hello3::ID = 0;    
//static RegisterPass<Hello3> Z("hello3", "CallGraph Pass");
//
//namespace {
//
//struct Hello4: public CallGraphSCCPass {
//  static char ID;
//
//  Hello4() : CallGraphSCCPass(ID) {
//    //initializeAttributorCGSCCLegacyPassPass(*PassRegistry::getPassRegistry());
//  }
//
//  bool runOnSCC(CallGraphSCC &SCC) override {
//    if (skipSCC(SCC))
//    {
//      errs() << "SCC pass skipped hello3";
//      return false;
//    }
//
//    SetVector<Function *> Functions;
//    for (CallGraphNode *CGN : SCC)
//    {
//      errs() << "node: " << CGN << "\n";
//      if (Function *Fn = CGN->getFunction())
//      {
//        errs() << "uses: " << Fn->uses();
//        errs() << "\t" << "getFunction:\n " << *Fn << "\n";
//        if (!Fn->isDeclaration()){
//          errs() << "this is not a function declaration\n";
//          Functions.insert(Fn);
//        }
//      }
//
//    }
//
//    if (Functions.empty())
//      return false;
//
//    AnalysisGetter AG;
//    CallGraph &CG = const_cast<CallGraph &>(SCC.getCallGraph());
//    CallGraphUpdater CGUpdater;
//    CGUpdater.initialize(CG, SCC);
//    Module &M = *Functions.back()->getParent();
//    BumpPtrAllocator Allocator;
//    InformationCache InfoCache(M, AG, Allocator, /* CGSCC */ &Functions);
//    //return runAttributorOnFunctions(InfoCache, Functions, AG, CGUpdater,
//    //                                /* DeleteFns */ false,
//    //                                /* IsModulePass */ false);
//  }
//
//  void getAnalysisUsage(AnalysisUsage &AU) const override {
//    // FIXME: Think about passes we will preserve and add them here.
//    AU.addRequired<TargetLibraryInfoWrapperPass>();
//    CallGraphSCCPass::getAnalysisUsage(AU);
//  }
//};
//}
//char Hello4::ID = 0;    
//static RegisterPass<Hello4> W("hello4", "CallGraph Pass");
//
