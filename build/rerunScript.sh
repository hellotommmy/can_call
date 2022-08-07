export PATH="~/llvm-project/build/bin:${PATH}"
#clang -O0 -emit-llvm ~/llvm-project/llvm/lib/Transforms/Hello/hi.c -c -o ~/llvm-project/llvm/lib/Transforms/Hello/hi.bc
cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DLLVM_USE_LINKER=gold ../llvm
cmake --build .
./bin/clang -O0 -disable-llvm-passes -disable-llvm-optzns hi.c
./bin/clang -O1 -disable-llvm-passes -disable-llvm-optzns hi.c
./bin/clang -O2 -disable-llvm-passes -disable-llvm-optzns hi.c
./bin/clang -O3 -disable-llvm-passes -disable-llvm-optzns hi.c
#bin/opt -enable-new-pm=0 -load lib/LLVMHello.so -passes=-called-value-propagation1 < ~/llvm-project/llvm/lib/Transforms/Hello/hi.bc  > /dev/null 
#bin/opt -enable-new-pm=0 -load lib/LLVMHello.so -passes=-called-value-propagation1 < ~/llvm-project/llvm/lib/Transforms/Hello/hi.bc  > /dev/null 
#bin/opt -enable-new-pm=0 -load lib/LLVMHello.so -passes=-called-value-propagation1 < ~/llvm-project/llvm/lib/Transforms/Hello/hi.bc  > /dev/null 
