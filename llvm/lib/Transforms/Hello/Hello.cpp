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

/// The maximum number of functions to track per lattice value. Once the number
/// of functions a call site can possibly target exceeds this threshold, it's
/// lattice value becomes overdefined. The number of possible lattice values is
/// bounded by Ch(F, M), where F is the number of functions in the module and M
/// is MaxFunctionsPerValue1. As such, this value should be kept very small. We
/// likely can't do anything useful for call sites with a large number of
/// possible targets, anyway.
static cl::opt<unsigned> MaxFunctionsPerValue1(
    "cvp-max-functions-per-value1", cl::Hidden, cl::init(4),
    cl::desc("The maximum number of functions to track per lattice value"));

namespace {
/// To enable interprocedural analysis, we assign LLVM values to the following
/// groups. The register group represents SSA registers, the return group
/// represents the return values of functions, and the memory group represents
/// in-memory values. An LLVM Value can technically be in more than one group.
/// It's necessary to distinguish these groups so we can, for example, track a
/// global variable separately from the value stored at its location.
enum class IPOGrouping1 { Register, Return, Memory };

/// Our LatticeKeys are PointerIntPairs composed of LLVM values and groupings.
using CVPLatticeKey1 = PointerIntPair<Value *, 2, IPOGrouping1>;

/// The lattice value type used by our custom lattice function. It holds the
/// lattice state, and a set of functions.
class CVPLatticeVal1 {
public:
  /// The states of the lattice values. Only the FunctionSet state is
  /// interesting. It indicates the set of functions to which an LLVM value may
  /// refer.
  enum CVPLatticeStateTy { Undefined, FunctionSet, Overdefined, Untracked };

  /// Comparator for sorting the functions set. We want to keep the order
  /// deterministic for testing, etc.
  struct Compare {
    bool operator()(const Function *LHS, const Function *RHS) const {
      return LHS->getName() < RHS->getName();
    }
  };

  CVPLatticeVal1() = default;
  CVPLatticeVal1(CVPLatticeStateTy LatticeState) : LatticeState(LatticeState) {}
  CVPLatticeVal1(std::vector<Function *> &&Functions)
      : LatticeState(FunctionSet), Functions(std::move(Functions)) {
    assert(llvm::is_sorted(this->Functions, Compare()));
  }

  /// Get a reference to the functions held by this lattice value. The number
  /// of functions will be zero for states other than FunctionSet.
  const std::vector<Function *> &getFunctions() const {
    return Functions;
  }

  /// Returns true if the lattice value is in the FunctionSet state.
  bool isFunctionSet() const { return LatticeState == FunctionSet; }

  bool operator==(const CVPLatticeVal1 &RHS) const {
    return LatticeState == RHS.LatticeState && Functions == RHS.Functions;
  }

  bool operator!=(const CVPLatticeVal1 &RHS) const {
    return LatticeState != RHS.LatticeState || Functions != RHS.Functions;
  }

private:
  /// Holds the state this lattice value is in.
  CVPLatticeStateTy LatticeState = Undefined;

  /// Holds functions indicating the possible targets of call sites. This set
  /// is empty for lattice values in the undefined, overdefined, and untracked
  /// states. The maximum size of the set is controlled by
  /// MaxFunctionsPerValue1. Since most LLVM values are expected to be in
  /// uninteresting states (i.e., overdefined), CVPLatticeVal1 objects should be
  /// small and efficiently copyable.
  // FIXME: This could be a TinyPtrVector and/or merge with LatticeState.
  std::vector<Function *> Functions;
};

/// The custom lattice function used by the generic sparse propagation solver.
/// It handles merging lattice values and computing new lattice values for
/// constants, arguments, values returned from trackable functions, and values
/// located in trackable global variables. It also computes the lattice values
/// that change as a result of executing instructions.
class CVPLatticeFunc1
    : public AbstractLatticeFunction<CVPLatticeKey1, CVPLatticeVal1> {
public:
  CVPLatticeFunc1()
      : AbstractLatticeFunction(CVPLatticeVal1(CVPLatticeVal1::Undefined),
                                CVPLatticeVal1(CVPLatticeVal1::Overdefined),
                                CVPLatticeVal1(CVPLatticeVal1::Untracked)) {}

  /// Compute and return a CVPLatticeVal1 for the given CVPLatticeKey1.
  CVPLatticeVal1 ComputeLatticeVal(CVPLatticeKey1 Key) override {
    switch (Key.getInt()) {
    case IPOGrouping1::Register:
      if (isa<Instruction>(Key.getPointer())) {
        return getUndefVal();
      } else if (auto *A = dyn_cast<Argument>(Key.getPointer())) {
        if (canTrackArgumentsInterprocedurally(A->getParent()))
          return getUndefVal();
      } else if (auto *C = dyn_cast<Constant>(Key.getPointer())) {
        return computeConstant(C);
      }
      return getOverdefinedVal();
    case IPOGrouping1::Memory:
    case IPOGrouping1::Return:
      if (auto *GV = dyn_cast<GlobalVariable>(Key.getPointer())) {
        if (canTrackGlobalVariableInterprocedurally(GV))
          return computeConstant(GV->getInitializer());
      } else if (auto *F = cast<Function>(Key.getPointer()))
        if (canTrackReturnsInterprocedurally(F))
          return getUndefVal();
    }
    return getOverdefinedVal();
  }

  /// Merge the two given lattice values. The interesting cases are merging two
  /// FunctionSet values and a FunctionSet value with an Undefined value. For
  /// these cases, we simply union the function sets. If the size of the union
  /// is greater than the maximum functions we track, the merged value is
  /// overdefined.
  CVPLatticeVal1 MergeValues(CVPLatticeVal1 X, CVPLatticeVal1 Y) override {
    if (X == getOverdefinedVal() || Y == getOverdefinedVal())
      return getOverdefinedVal();
    if (X == getUndefVal() && Y == getUndefVal())
      return getUndefVal();
    std::vector<Function *> Union;
    errs() << "Merging\n";
    for(auto x: X.getFunctions()){
      

      x->print(errs());

    }

    for(auto y: Y.getFunctions()){
      y->print(errs()); 
    }
    std::set_union(X.getFunctions().begin(), X.getFunctions().end(),
                   Y.getFunctions().begin(), Y.getFunctions().end(),
                   std::back_inserter(Union), CVPLatticeVal1::Compare{});
    if (Union.size() > MaxFunctionsPerValue1)
      return getOverdefinedVal();
    return CVPLatticeVal1(std::move(Union));
  }

  /// Compute the lattice values that change as a result of executing the given
  /// instruction. The changed values are stored in \p ChangedValues. We handle
  /// just a few kinds of instructions since we're only propagating values that
  /// can be called.
  void ComputeInstructionState(
      Instruction &I, DenseMap<CVPLatticeKey1, CVPLatticeVal1> &ChangedValues,
      SparseSolver<CVPLatticeKey1, CVPLatticeVal1> &SS) override {
    switch (I.getOpcode()) {
    case Instruction::Call:
    case Instruction::Invoke:
      return visitCallBase(cast<CallBase>(I), ChangedValues, SS);
    case Instruction::Load:
      return visitLoad(*cast<LoadInst>(&I), ChangedValues, SS);
    case Instruction::Ret:
      return visitReturn(*cast<ReturnInst>(&I), ChangedValues, SS);
    case Instruction::Select:
      return visitSelect(*cast<SelectInst>(&I), ChangedValues, SS);
    case Instruction::Store:
      return visitStore(*cast<StoreInst>(&I), ChangedValues, SS);
    default:
      return visitInst(I, ChangedValues, SS);
    }
  }

  /// Print the given CVPLatticeVal1 to the specified stream.
  void PrintLatticeVal(CVPLatticeVal1 LV, raw_ostream &OS) override {
    if (LV == getUndefVal())
      OS << "Undefined  ";
    else if (LV == getOverdefinedVal())
      OS << "Overdefined";
    else if (LV == getUntrackedVal())
      OS << "Untracked  ";
    else
      OS << "FunctionSet";
  }

  /// Print the given CVPLatticeKey1 to the specified stream.
  void PrintLatticeKey(CVPLatticeKey1 Key, raw_ostream &OS) override {
    if (Key.getInt() == IPOGrouping1::Register)
      OS << "<reg> ";
    else if (Key.getInt() == IPOGrouping1::Memory)
      OS << "<mem> ";
    else if (Key.getInt() == IPOGrouping1::Return)
      OS << "<ret> ";
    if (isa<Function>(Key.getPointer()))
      OS << Key.getPointer()->getName();
    else
      OS << *Key.getPointer();
  }

  /// We collect a set of indirect calls when visiting call sites. This method
  /// returns a reference to that set.
  SmallPtrSetImpl<CallBase *> &getIndirectCalls() { return IndirectCalls; }

private:
  /// Holds the indirect calls we encounter during the analysis. We will attach
  /// metadata to these calls after the analysis indicating the functions the
  /// calls can possibly target.
  SmallPtrSet<CallBase *, 32> IndirectCalls;

  /// Compute a new lattice value for the given constant. The constant, after
  /// stripping any pointer casts, should be a Function. We ignore null
  /// pointers as an optimization, since calling these values is undefined
  /// behavior.
  CVPLatticeVal1 computeConstant(Constant *C) {
    if (isa<ConstantPointerNull>(C))
      return CVPLatticeVal1(CVPLatticeVal1::FunctionSet);
    if (auto *F = dyn_cast<Function>(C->stripPointerCasts()))
      return CVPLatticeVal1({F});
    return getOverdefinedVal();
  }

  /// Handle return instructions. The function's return state is the merge of
  /// the returned value state and the function's return state.
  void visitReturn(ReturnInst &I,
                   DenseMap<CVPLatticeKey1, CVPLatticeVal1> &ChangedValues,
                   SparseSolver<CVPLatticeKey1, CVPLatticeVal1> &SS) {
    Function *F = I.getParent()->getParent();
    if (F->getReturnType()->isVoidTy())
      return;
    auto RegI = CVPLatticeKey1(I.getReturnValue(), IPOGrouping1::Register);
    auto RetF = CVPLatticeKey1(F, IPOGrouping1::Return);
    ChangedValues[RetF] =
        MergeValues(SS.getValueState(RegI), SS.getValueState(RetF));
  }

  /// Handle call sites. The state of a called function's formal arguments is
  /// the merge of the argument state with the call sites corresponding actual
  /// argument state. The call site state is the merge of the call site state
  /// with the returned value state of the called function.
  void visitCallBase(CallBase &CB,
                     DenseMap<CVPLatticeKey1, CVPLatticeVal1> &ChangedValues,
                     SparseSolver<CVPLatticeKey1, CVPLatticeVal1> &SS) {
    Function *F = CB.getCalledFunction();
    auto RegI = CVPLatticeKey1(&CB, IPOGrouping1::Register);

    // If this is an indirect call, save it so we can quickly revisit it when
    // attaching metadata.
    if (!F){
      errs() << "This should happen for the func ptr call\n";
      IndirectCalls.insert(&CB);
    }

    // If we can't track the function's return values, there's nothing to do.
    if (!F || !canTrackReturnsInterprocedurally(F)) {
      errs() << "can't track return values\n";
      // Void return, No need to create and update CVPLattice state as no one
      // can use it.
      if (CB.getType()->isVoidTy())
        return;
      ChangedValues[RegI] = getOverdefinedVal();
      return;
    }

    // Inform the solver that the called function is executable, and perform
    // the merges for the arguments and return value.
    SS.MarkBlockExecutable(&F->front());
    auto RetF = CVPLatticeKey1(F, IPOGrouping1::Return);
    for (Argument &A : F->args()) {
      auto RegFormal = CVPLatticeKey1(&A, IPOGrouping1::Register);
      auto RegActual =
          CVPLatticeKey1(CB.getArgOperand(A.getArgNo()), IPOGrouping1::Register);
      ChangedValues[RegFormal] =
          MergeValues(SS.getValueState(RegFormal), SS.getValueState(RegActual));
    }

    // Void return, No need to create and update CVPLattice state as no one can
    // use it.
    if (CB.getType()->isVoidTy())
      return;

    ChangedValues[RegI] =
        MergeValues(SS.getValueState(RegI), SS.getValueState(RetF));
  }

  /// Handle select instructions. The select instruction state is the merge the
  /// true and false value states.
  void visitSelect(SelectInst &I,
                   DenseMap<CVPLatticeKey1, CVPLatticeVal1> &ChangedValues,
                   SparseSolver<CVPLatticeKey1, CVPLatticeVal1> &SS) {
    auto RegI = CVPLatticeKey1(&I, IPOGrouping1::Register);
    auto RegT = CVPLatticeKey1(I.getTrueValue(), IPOGrouping1::Register);
    auto RegF = CVPLatticeKey1(I.getFalseValue(), IPOGrouping1::Register);
    ChangedValues[RegI] =
        MergeValues(SS.getValueState(RegT), SS.getValueState(RegF));
  }

  /// Handle load instructions. If the pointer operand of the load is a global
  /// variable, we attempt to track the value. The loaded value state is the
  /// merge of the loaded value state with the global variable state.
  void visitLoad(LoadInst &I,
                 DenseMap<CVPLatticeKey1, CVPLatticeVal1> &ChangedValues,
                 SparseSolver<CVPLatticeKey1, CVPLatticeVal1> &SS) {
    auto RegI = CVPLatticeKey1(&I, IPOGrouping1::Register);
    if (auto *GV = dyn_cast<GlobalVariable>(I.getPointerOperand())) {
      auto MemGV = CVPLatticeKey1(GV, IPOGrouping1::Memory);
      ChangedValues[RegI] =
          MergeValues(SS.getValueState(RegI), SS.getValueState(MemGV));
    } else {
      ChangedValues[RegI] = getOverdefinedVal();
    }
  }

  /// Handle store instructions. If the pointer operand of the store is a
  /// global variable, we attempt to track the value. The global variable state
  /// is the merge of the stored value state with the global variable state.
  void visitStore(StoreInst &I,
                  DenseMap<CVPLatticeKey1, CVPLatticeVal1> &ChangedValues,
                  SparseSolver<CVPLatticeKey1, CVPLatticeVal1> &SS) {
    auto *GV = dyn_cast<GlobalVariable>(I.getPointerOperand());
    if (!GV){
      errs() << "do not track non-global values\n";
      I.print(errs());

      errs() << "operand is:" << "\n";
      I.getValueOperand()->dump();
      errs() << "pt operand is : " << "\n";
      I.getPointerOperand()->dump();
      return;
    }
    auto RegI = CVPLatticeKey1(I.getValueOperand(), IPOGrouping1::Register);
    auto MemGV = CVPLatticeKey1(GV, IPOGrouping1::Memory);
    ChangedValues[MemGV] =
        MergeValues(SS.getValueState(RegI), SS.getValueState(MemGV));
  }

  /// Handle all other instructions. All other instructions are marked
  /// overdefined.
  void visitInst(Instruction &I,
                 DenseMap<CVPLatticeKey1, CVPLatticeVal1> &ChangedValues,
                 SparseSolver<CVPLatticeKey1, CVPLatticeVal1> &SS) {
    // Simply bail if this instruction has no user.
    if (I.use_empty())
      return;
    auto RegI = CVPLatticeKey1(&I, IPOGrouping1::Register);
    ChangedValues[RegI] = getOverdefinedVal();
  }
};
} // namespace

namespace llvm {
/// A specialization of LatticeKeyInfo for CVPLatticeKeys. The generic solver
/// must translate between LatticeKeys and LLVM Values when adding Values to
/// its work list and inspecting the state of control-flow related values.
template <> struct LatticeKeyInfo<CVPLatticeKey1> {
  static inline Value *getValueFromLatticeKey(CVPLatticeKey1 Key) {
    return Key.getPointer();
  }
  static inline CVPLatticeKey1 getLatticeKeyFromValue(Value *V) {
    return CVPLatticeKey1(V, IPOGrouping1::Register);
  }
};
} // namespace llvm

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
  // Our custom lattice function and generic sparse propagation solver.
  //CVPLatticeFunc1 Lattice;
  //SparseSolver<CVPLatticeKey1, CVPLatticeVal1> Solver(&Lattice);

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
// For each function in the module, if we can't track its arguments, let the
// generic solver assume it is executable.
  //for (Function &F : M)
  //  if (!F.isDeclaration() && !canTrackArgumentsInterprocedurally(&F))
  //    Solver.MarkBlockExecutable(&F.front());

  // Solver our custom lattice. In doing so, we will also build a set of
  // indirect call sites.
  //Solver.Solve();

  // Attach metadata to the indirect call sites that were collected indicating
  // the set of functions they can possibly target.
  //bool Changed = false;
  //MDBuilder MDB(M.getContext());
  //for (CallBase *C : Lattice.getIndirectCalls()) {
  //  errs() << C->getFunctionType();
  //  errs() << C->getCalledOperand();
  //  auto RegI = CVPLatticeKey1(C->getCalledOperand(), IPOGrouping1::Register);
  //  CVPLatticeVal1 LV = Solver.getExistingValueState(RegI);
    
  //  for (auto f: LV.getFunctions()) {

  //    errs() << "HI!!!!!\n";
  //    errs() << f;
  //  }
    
    
  //  if (!LV.isFunctionSet() || LV.getFunctions().empty())
  //    continue;
  //  MDNode *Callees = MDB.createCallees(LV.getFunctions());
  //  C->setMetadata(LLVMContext::MD_callees, Callees);
  //  Changed = true;
  //}

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
