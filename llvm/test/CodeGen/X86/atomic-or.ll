; NOTE: Assertions have been autogenerated by utils/update_llc_test_checks.py
; RUN: llc < %s -mtriple=x86_64-- -verify-machineinstrs | FileCheck %s

; rdar://9692967

define void @t1(ptr %p, i32 %b) nounwind {
; CHECK-LABEL: t1:
; CHECK:       # %bb.0: # %entry
; CHECK-NEXT:    movq %rdi, -{{[0-9]+}}(%rsp)
; CHECK-NEXT:    movl $2147483648, %eax # imm = 0x80000000
; CHECK-NEXT:    lock orq %rax, (%rdi)
; CHECK-NEXT:    retq
entry:
  %p.addr = alloca ptr, align 8
  store ptr %p, ptr %p.addr, align 8
  %tmp = load ptr, ptr %p.addr, align 8
  %0 = atomicrmw or ptr %tmp, i64 2147483648 seq_cst
  ret void
}

define void @t2(ptr %p, i32 %b) nounwind {
; CHECK-LABEL: t2:
; CHECK:       # %bb.0: # %entry
; CHECK-NEXT:    movq %rdi, -{{[0-9]+}}(%rsp)
; CHECK-NEXT:    lock orq $2147483644, (%rdi) # imm = 0x7FFFFFFC
; CHECK-NEXT:    retq
entry:
  %p.addr = alloca ptr, align 8
  store ptr %p, ptr %p.addr, align 8
  %tmp = load ptr, ptr %p.addr, align 8
  %0 = atomicrmw or ptr %tmp, i64 2147483644 seq_cst
  ret void
}
