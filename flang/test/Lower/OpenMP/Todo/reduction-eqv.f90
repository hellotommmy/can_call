! RUN: %not_todo_cmd bbc -emit-fir -fopenmp -o - %s 2>&1 | FileCheck %s
! RUN: %not_todo_cmd %flang_fc1 -emit-fir -fopenmp -o - %s 2>&1 | FileCheck %s

! CHECK: not yet implemented: Reduction of some intrinsic operators is not supported
subroutine reduction_eqv(y)
  logical :: x, y(100)
  !$omp parallel
  !$omp do reduction(.eqv.:x)
  do i=1, 100
    x = x .eqv. y(i)
  end do
  !$omp end do
  !$omp end parallel
  print *, x
end subroutine
