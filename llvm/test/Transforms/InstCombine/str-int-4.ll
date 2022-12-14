; NOTE: Assertions have been autogenerated by utils/update_test_checks.py
; Verify that calls to strtol and strtoll are interpreted correctly even
; in corner cases (or not folded).
;
; RUN: opt < %s -passes=instcombine -S | FileCheck %s

declare i32 @strtol(i8*, i8**, i32)
declare i64 @strtoll(i8*, i8**, i32)


; All POSIX whitespace characters.
@ws = constant [7 x i8] c"\09\0d\0a\0b\0c \00"

; A negative and positive number preceded by all POSIX whitespace.
@ws_im123 = constant [11 x i8] c"\09\0d\0a\0b\0c -123\00"
@ws_ip234 = constant [11 x i8] c"\09\0d\0a\0b\0c +234\00"

@wsplus = constant [3 x i8] c" +\00"
@i0 = constant [3 x i8] c" 0\00"
@i8 = constant [3 x i8] c" 8\00"
@i9 = constant [3 x i8] c" 9\00"
@ia = constant [3 x i8] c" a\00"
@i08 = constant [3 x i8] c"08\00"
@x0x = constant [3 x i8] c"0x\00"
@wsplusws0 = constant [5 x i8] c" + 0\00"
@i19azAZ = constant [7 x i8] c"19azAZ\00"
@i32min = constant [13 x i8] c" -2147483648\00"
@i32min_m1 = constant [13 x i8] c" -2147483649\00"
@o32min = constant [15 x i8] c" +020000000000\00"
@mo32min = constant [15 x i8] c" -020000000000\00"
@x32min = constant [13 x i8] c" +0x80000000\00"
@mx32min = constant [13 x i8] c" -0x80000000\00"

@i32max = constant [12 x i8] c" 2147483647\00"
@x32max = constant [12 x i8] c" 0x7fffffff\00"
@i32max_p1 = constant [12 x i8] c" 2147483648\00"

@ui32max = constant [12 x i8] c" 4294967295\00"
@ui32max_p1 = constant [12 x i8] c" 4294967296\00"


; Exercise folding calls to 32-bit strtol.

define void @fold_strtol(i32* %ps) {
; CHECK-LABEL: @fold_strtol(
; CHECK-NEXT:    store i32 -123, i32* [[PS:%.*]], align 4
; CHECK-NEXT:    [[PS1:%.*]] = getelementptr i32, i32* [[PS]], i64 1
; CHECK-NEXT:    store i32 234, i32* [[PS1]], align 4
; CHECK-NEXT:    [[PS2:%.*]] = getelementptr i32, i32* [[PS]], i64 2
; CHECK-NEXT:    store i32 0, i32* [[PS2]], align 4
; CHECK-NEXT:    [[PS3:%.*]] = getelementptr i32, i32* [[PS]], i64 3
; CHECK-NEXT:    store i32 9, i32* [[PS3]], align 4
; CHECK-NEXT:    [[PS4:%.*]] = getelementptr i32, i32* [[PS]], i64 4
; CHECK-NEXT:    store i32 10, i32* [[PS4]], align 4
; CHECK-NEXT:    [[PS5:%.*]] = getelementptr i32, i32* [[PS]], i64 5
; CHECK-NEXT:    store i32 76095035, i32* [[PS5]], align 4
; CHECK-NEXT:    [[PS6:%.*]] = getelementptr i32, i32* [[PS]], i64 6
; CHECK-NEXT:    store i32 -2147483648, i32* [[PS6]], align 4
; CHECK-NEXT:    [[PS7:%.*]] = getelementptr i32, i32* [[PS]], i64 7
; CHECK-NEXT:    store i32 -2147483648, i32* [[PS7]], align 4
; CHECK-NEXT:    [[PS8:%.*]] = getelementptr i32, i32* [[PS]], i64 8
; CHECK-NEXT:    store i32 -2147483648, i32* [[PS8]], align 4
; CHECK-NEXT:    [[PS9:%.*]] = getelementptr i32, i32* [[PS]], i64 9
; CHECK-NEXT:    store i32 -2147483648, i32* [[PS9]], align 4
; CHECK-NEXT:    [[PS10:%.*]] = getelementptr i32, i32* [[PS]], i64 10
; CHECK-NEXT:    store i32 2147483647, i32* [[PS10]], align 4
; CHECK-NEXT:    [[PS11:%.*]] = getelementptr i32, i32* [[PS]], i64 11
; CHECK-NEXT:    store i32 2147483647, i32* [[PS11]], align 4
; CHECK-NEXT:    ret void
;
; Fold a valid sequence with leading POSIX whitespace and a minus to -123.
  %pwsm123 = getelementptr [11 x i8], [11 x i8]* @ws_im123, i32 0, i32 0
  %im123 = call i32 @strtol(i8* %pwsm123, i8** null, i32 10)
  %ps0 = getelementptr i32, i32* %ps, i32 0
  store i32 %im123, i32* %ps0

; Fold a valid sequence with leading POSIX whitespace and a plus to +234.
  %pwsp234 = getelementptr [11 x i8], [11 x i8]* @ws_ip234, i32 0, i32 0
  %ip234 = call i32 @strtol(i8* %pwsp234, i8** null, i32 10)
  %ps1 = getelementptr i32, i32* %ps, i32 1
  store i32 %ip234, i32* %ps1

; Fold " 0" in base 0 to verify correct base autodetection.
  %psi0 = getelementptr [3 x i8], [3 x i8]* @i0, i32 0, i32 0
  %i0 = call i32 @strtol(i8* %psi0, i8** null, i32 0)
  %ps2 = getelementptr i32, i32* %ps, i32 2
  store i32 %i0, i32* %ps2

; Fold " 9" in base 0 to verify correct base autodetection.
  %psi9 = getelementptr [3 x i8], [3 x i8]* @i9, i32 0, i32 0
  %i9 = call i32 @strtol(i8* %psi9, i8** null, i32 0)
  %ps3 = getelementptr i32, i32* %ps, i32 3
  store i32 %i9, i32* %ps3

; Fold " a" in base 16 to 10.
  %psia = getelementptr [3 x i8], [3 x i8]* @ia, i32 0, i32 0
  %ia = call i32 @strtol(i8* %psia, i8** null, i32 16)
  %ps4 = getelementptr i32, i32* %ps, i32 4
  store i32 %ia, i32* %ps4

; Fold "19azAZ" in base 36 to 76095035.
  %psi19azAZ = getelementptr [7 x i8], [7 x i8]* @i19azAZ, i32 0, i32 0
  %i19azAZ = call i32 @strtol(i8* %psi19azAZ, i8** null, i32 36)
  %ps5 = getelementptr i32, i32* %ps, i32 5
  store i32 %i19azAZ, i32* %ps5

; Fold INT32_MIN.
  %psmin = getelementptr [13 x i8], [13 x i8]* @i32min, i32 0, i32 0
  %min = call i32 @strtol(i8* %psmin, i8** null, i32 10)
  %ps6 = getelementptr i32, i32* %ps, i32 6
  store i32 %min, i32* %ps6

; Fold -INT32_MIN in octal.
  %psmo32min = getelementptr [15 x i8], [15 x i8]* @mo32min, i32 0, i32 0
  %mo32min = call i32 @strtol(i8* %psmo32min, i8** null, i32 0)
  %ps7 = getelementptr i32, i32* %ps, i32 7
  store i32 %mo32min, i32* %ps7

; Fold -INT32_MIN in hex and base 0.
  %psmx32min = getelementptr [13 x i8], [13 x i8]* @mx32min, i32 0, i32 0
  %mx32min_0 = call i32 @strtol(i8* %psmx32min, i8** null, i32 0)
  %ps8 = getelementptr i32, i32* %ps, i32 8
  store i32 %mx32min_0, i32* %ps8

; Fold -INT32_MIN in hex and base 16.
  %mx32min_16 = call i32 @strtol(i8* %psmx32min, i8** null, i32 16)
  %ps9 = getelementptr i32, i32* %ps, i32 9
  store i32 %mx32min_16, i32* %ps9

; Fold INT32_MAX.
  %psmax = getelementptr [12 x i8], [12 x i8]* @i32max, i32 0, i32 0
  %max = call i32 @strtol(i8* %psmax, i8** null, i32 10)
  %ps10 = getelementptr i32, i32* %ps, i32 10
  store i32 %max, i32* %ps10

; Fold INT32_MAX in hex.
  %psxmax = getelementptr [12 x i8], [12 x i8]* @x32max, i32 0, i32 0
  %xmax = call i32 @strtol(i8* %psxmax, i8** null, i32 0)
  %ps11 = getelementptr i32, i32* %ps, i32 11
  store i32 %xmax, i32* %ps11

  ret void
}


; Exercise not folding calls to 32-bit strtol.

define void @call_strtol(i32* %ps) {
; CHECK-LABEL: @call_strtol(
; CHECK-NEXT:    [[MINM1:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([13 x i8], [13 x i8]* @i32min_m1, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    store i32 [[MINM1]], i32* [[PS:%.*]], align 4
; CHECK-NEXT:    [[MAXP1:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([12 x i8], [12 x i8]* @i32max_p1, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    [[PS1:%.*]] = getelementptr i32, i32* [[PS]], i64 1
; CHECK-NEXT:    store i32 [[MAXP1]], i32* [[PS1]], align 4
; CHECK-NEXT:    [[IPLUS:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([3 x i8], [3 x i8]* @wsplus, i64 0, i64 0), i8** null, i32 0)
; CHECK-NEXT:    [[PS2:%.*]] = getelementptr i32, i32* [[PS]], i64 2
; CHECK-NEXT:    store i32 [[IPLUS]], i32* [[PS2]], align 4
; CHECK-NEXT:    [[IA:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([3 x i8], [3 x i8]* @ia, i64 0, i64 0), i8** null, i32 0)
; CHECK-NEXT:    [[PS3:%.*]] = getelementptr i32, i32* [[PS]], i64 3
; CHECK-NEXT:    store i32 [[IA]], i32* [[PS3]], align 4
; CHECK-NEXT:    [[I8:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([3 x i8], [3 x i8]* @i8, i64 0, i64 0), i8** null, i32 8)
; CHECK-NEXT:    [[PS4:%.*]] = getelementptr i32, i32* [[PS]], i64 4
; CHECK-NEXT:    store i32 [[I8]], i32* [[PS4]], align 4
; CHECK-NEXT:    [[I0X:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([3 x i8], [3 x i8]* @x0x, i64 0, i64 0), i8** null, i32 0)
; CHECK-NEXT:    [[PS5:%.*]] = getelementptr i32, i32* [[PS]], i64 5
; CHECK-NEXT:    store i32 [[I0X]], i32* [[PS5]], align 4
; CHECK-NEXT:    [[IWSPWS0:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([5 x i8], [5 x i8]* @wsplusws0, i64 0, i64 0), i8** null, i32 0)
; CHECK-NEXT:    [[PS6:%.*]] = getelementptr i32, i32* [[PS]], i64 6
; CHECK-NEXT:    store i32 [[IWSPWS0]], i32* [[PS6]], align 4
; CHECK-NEXT:    [[I19AZAZ:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([7 x i8], [7 x i8]* @i19azAZ, i64 0, i64 0), i8** null, i32 35)
; CHECK-NEXT:    [[PS7:%.*]] = getelementptr i32, i32* [[PS]], i64 7
; CHECK-NEXT:    store i32 [[I19AZAZ]], i32* [[PS7]], align 4
; CHECK-NEXT:    [[O32MIN:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([15 x i8], [15 x i8]* @o32min, i64 0, i64 0), i8** null, i32 0)
; CHECK-NEXT:    [[PS8:%.*]] = getelementptr i32, i32* [[PS]], i64 8
; CHECK-NEXT:    store i32 [[O32MIN]], i32* [[PS8]], align 4
; CHECK-NEXT:    [[X32MIN:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([13 x i8], [13 x i8]* @x32min, i64 0, i64 0), i8** null, i32 0)
; CHECK-NEXT:    [[PS9:%.*]] = getelementptr i32, i32* [[PS]], i64 9
; CHECK-NEXT:    store i32 [[X32MIN]], i32* [[PS9]], align 4
; CHECK-NEXT:    [[X32MIN_10:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([13 x i8], [13 x i8]* @x32min, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    [[PS10:%.*]] = getelementptr i32, i32* [[PS]], i64 10
; CHECK-NEXT:    store i32 [[X32MIN_10]], i32* [[PS10]], align 4
; CHECK-NEXT:    [[NWS:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([7 x i8], [7 x i8]* @ws, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    [[PS11:%.*]] = getelementptr i32, i32* [[PS]], i64 11
; CHECK-NEXT:    store i32 [[NWS]], i32* [[PS11]], align 4
; CHECK-NEXT:    [[NWSP6:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([7 x i8], [7 x i8]* @ws, i64 0, i64 6), i8** null, i32 10)
; CHECK-NEXT:    [[PS12:%.*]] = getelementptr i32, i32* [[PS]], i64 12
; CHECK-NEXT:    store i32 [[NWSP6]], i32* [[PS12]], align 4
; CHECK-NEXT:    [[I0B1:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([3 x i8], [3 x i8]* @i0, i64 0, i64 0), i8** null, i32 1)
; CHECK-NEXT:    [[PS13:%.*]] = getelementptr i32, i32* [[PS]], i64 13
; CHECK-NEXT:    store i32 [[I0B1]], i32* [[PS13]], align 4
; CHECK-NEXT:    [[I0B256:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([3 x i8], [3 x i8]* @i0, i64 0, i64 0), i8** null, i32 256)
; CHECK-NEXT:    [[PS14:%.*]] = getelementptr i32, i32* [[PS]], i64 14
; CHECK-NEXT:    store i32 [[I0B256]], i32* [[PS14]], align 4
; CHECK-NEXT:    ret void
;

; Do not fold the result of conversion that's less than INT32_MIN.
  %psminm1 = getelementptr [13 x i8], [13 x i8]* @i32min_m1, i32 0, i32 0
  %minm1 = call i32 @strtol(i8* %psminm1, i8** null, i32 10)
  %ps0 = getelementptr i32, i32* %ps, i32 0
  store i32 %minm1, i32* %ps0

; Do not fold the result of conversion that's greater than INT32_MAX.
  %psmaxp1 = getelementptr [12 x i8], [12 x i8]* @i32max_p1, i32 0, i32 0
  %maxp1 = call i32 @strtol(i8* %psmaxp1, i8** null, i32 10)
  %ps1 = getelementptr i32, i32* %ps, i32 1
  store i32 %maxp1, i32* %ps1

; Do not fold " +".
  %psplus = getelementptr [3 x i8], [3 x i8]* @wsplus, i32 0, i32 0
  %iplus = call i32 @strtol(i8* %psplus, i8** null, i32 0)
  %ps2 = getelementptr i32, i32* %ps, i32 2
  store i32 %iplus, i32* %ps2

; Do not fold " a" in base 0.
  %psia = getelementptr [3 x i8], [3 x i8]* @ia, i32 0, i32 0
  %ia = call i32 @strtol(i8* %psia, i8** null, i32 0)
  %ps3 = getelementptr i32, i32* %ps, i32 3
  store i32 %ia, i32* %ps3

; Do not fold " 8" in base 8.
  %psi8 = getelementptr [3 x i8], [3 x i8]* @i8, i32 0, i32 0
  %i8 = call i32 @strtol(i8* %psi8, i8** null, i32 8)
  %ps4 = getelementptr i32, i32* %ps, i32 4
  store i32 %i8, i32* %ps4

; Do not fold the "0x" alone in base 0 that some implementations (e.g.,
; BSD and Darwin) set errno to EINVAL for.
  %psx0x = getelementptr [3 x i8], [3 x i8]* @x0x, i32 0, i32 0
  %i0x = call i32 @strtol(i8* %psx0x, i8** null, i32 0)
  %ps5 = getelementptr i32, i32* %ps, i32 5
  store i32 %i0x, i32* %ps5

; Do not fold " + 0".
  %pwspws0 = getelementptr [5 x i8], [5 x i8]* @wsplusws0, i32 0, i32 0
  %iwspws0 = call i32 @strtol(i8* %pwspws0, i8** null, i32 0)
  %ps6 = getelementptr i32, i32* %ps, i32 6
  store i32 %iwspws0, i32* %ps6

; Do not fold "19azAZ" in base 35.
  %psi19azAZ = getelementptr [7 x i8], [7 x i8]* @i19azAZ, i32 0, i32 0
  %i19azAZ = call i32 @strtol(i8* %psi19azAZ, i8** null, i32 35)
  %ps7 = getelementptr i32, i32* %ps, i32 7
  store i32 %i19azAZ, i32* %ps7

; Do not fold INT32_MIN in octal.
  %pso32min = getelementptr [15 x i8], [15 x i8]* @o32min, i32 0, i32 0
  %o32min = call i32 @strtol(i8* %pso32min, i8** null, i32 0)
  %ps8 = getelementptr i32, i32* %ps, i32 8
  store i32 %o32min, i32* %ps8

; Do not fold INT32_MIN in hex.
  %psx32min = getelementptr [13 x i8], [13 x i8]* @x32min, i32 0, i32 0
  %x32min = call i32 @strtol(i8* %psx32min, i8** null, i32 0)
  %ps9 = getelementptr i32, i32* %ps, i32 9
  store i32 %x32min, i32* %ps9

; Do not fold INT32_MIN in hex in base 10.
  %x32min_10 = call i32 @strtol(i8* %psx32min, i8** null, i32 10)
  %ps10 = getelementptr i32, i32* %ps, i32 10
  store i32 %x32min_10, i32* %ps10

; Do not fold a sequence consisting of just whitespace characters.
  %psws = getelementptr [7 x i8], [7 x i8]* @ws, i32 0, i32 0
  %nws = call i32 @strtol(i8* %psws, i8** null, i32 10)
  %ps11 = getelementptr i32, i32* %ps, i32 11
  store i32 %nws, i32* %ps11

; Do not fold an empty sequence.
  %pswsp6 = getelementptr [7 x i8], [7 x i8]* @ws, i32 0, i32 6
  %nwsp6 = call i32 @strtol(i8* %pswsp6, i8** null, i32 10)
  %ps12 = getelementptr i32, i32* %ps, i32 12
  store i32 %nwsp6, i32* %ps12

; Do not fold the invalid base 1.
  %psi0 = getelementptr [3 x i8], [3 x i8]* @i0, i32 0, i32 0
  %i0b1 = call i32 @strtol(i8* %psi0, i8** null, i32 1)
  %ps13 = getelementptr i32, i32* %ps, i32 13
  store i32 %i0b1, i32* %ps13

; Do not fold the invalid base 256.
  %i0b256 = call i32 @strtol(i8* %psi0, i8** null, i32 256)
  %ps14 = getelementptr i32, i32* %ps, i32 14
  store i32 %i0b256, i32* %ps14

  ret void
}


@i64min = constant [22 x i8] c" -9223372036854775808\00"
@i64min_m1 = constant [22 x i8] c" -9223372036854775809\00"

@i64max = constant [21 x i8] c" 9223372036854775807\00"
@i64max_p1 = constant [21 x i8] c" 9223372036854775808\00"

@ui64max = constant [22 x i8] c" 18446744073709551615\00"
@ui64max_p1 = constant [22 x i8] c" 18446744073709551616\00"


; Exercise folding calls to the 64-bit strtoll.

define void @fold_strtoll(i64* %ps) {
; CHECK-LABEL: @fold_strtoll(
; CHECK-NEXT:    store i64 -123, i64* [[PS:%.*]], align 4
; CHECK-NEXT:    [[PS1:%.*]] = getelementptr i64, i64* [[PS]], i64 1
; CHECK-NEXT:    store i64 234, i64* [[PS1]], align 4
; CHECK-NEXT:    [[PS2:%.*]] = getelementptr i64, i64* [[PS]], i64 2
; CHECK-NEXT:    store i64 -9223372036854775808, i64* [[PS2]], align 4
; CHECK-NEXT:    [[PS3:%.*]] = getelementptr i64, i64* [[PS]], i64 3
; CHECK-NEXT:    store i64 9223372036854775807, i64* [[PS3]], align 4
; CHECK-NEXT:    ret void
;
; Fold a valid sequence with leading POSIX whitespace and a minus to -123.
  %pwsm123 = getelementptr [11 x i8], [11 x i8]* @ws_im123, i32 0, i32 0
  %im123 = call i64 @strtoll(i8* %pwsm123, i8** null, i32 10)
  %ps0 = getelementptr i64, i64* %ps, i32 0
  store i64 %im123, i64* %ps0

; Fold a valid sequence with leading POSIX whitespace and a plus to +234.
  %pwsp234 = getelementptr [11 x i8], [11 x i8]* @ws_ip234, i32 0, i32 0
  %ip234 = call i64 @strtoll(i8* %pwsp234, i8** null, i32 10)
  %ps1 = getelementptr i64, i64* %ps, i32 1
  store i64 %ip234, i64* %ps1

; Fold INT64_MIN.
  %psmin = getelementptr [22 x i8], [22 x i8]* @i64min, i32 0, i32 0
  %min = call i64 @strtoll(i8* %psmin, i8** null, i32 10)
  %ps2 = getelementptr i64, i64* %ps, i32 2
  store i64 %min, i64* %ps2

; Fold INT64_MAX.
  %psmax = getelementptr [21 x i8], [21 x i8]* @i64max, i32 0, i32 0
  %max = call i64 @strtoll(i8* %psmax, i8** null, i32 10)
  %ps3 = getelementptr i64, i64* %ps, i32 3
  store i64 %max, i64* %ps3

  ret void
}


; Exercise not folding calls to the 64-bit strtoll.

define void @call_strtoll(i64* %ps) {
; CHECK-LABEL: @call_strtoll(
; CHECK-NEXT:    [[MINM1:%.*]] = call i64 @strtoll(i8* nocapture getelementptr inbounds ([22 x i8], [22 x i8]* @i64min_m1, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    store i64 [[MINM1]], i64* [[PS:%.*]], align 4
; CHECK-NEXT:    [[MAXP1:%.*]] = call i64 @strtoll(i8* nocapture getelementptr inbounds ([21 x i8], [21 x i8]* @i64max_p1, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    [[PS1:%.*]] = getelementptr i64, i64* [[PS]], i64 1
; CHECK-NEXT:    store i64 [[MAXP1]], i64* [[PS1]], align 4
; CHECK-NEXT:    [[NWS:%.*]] = call i64 @strtoll(i8* nocapture getelementptr inbounds ([7 x i8], [7 x i8]* @ws, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    [[PS2:%.*]] = getelementptr i64, i64* [[PS]], i64 2
; CHECK-NEXT:    store i64 [[NWS]], i64* [[PS2]], align 4
; CHECK-NEXT:    [[NWSP6:%.*]] = call i64 @strtoll(i8* nocapture getelementptr inbounds ([7 x i8], [7 x i8]* @ws, i64 0, i64 6), i8** null, i32 10)
; CHECK-NEXT:    [[PS3:%.*]] = getelementptr i64, i64* [[PS]], i64 3
; CHECK-NEXT:    store i64 [[NWSP6]], i64* [[PS3]], align 4
; CHECK-NEXT:    ret void
;
; Do not fold the result of conversion that's less than INT64_MIN.
  %psminm1 = getelementptr [22 x i8], [22 x i8]* @i64min_m1, i32 0, i32 0
  %minm1 = call i64 @strtoll(i8* %psminm1, i8** null, i32 10)
  %ps0 = getelementptr i64, i64* %ps, i32 0
  store i64 %minm1, i64* %ps0

; Do not fold the result of conversion that's greater than INT64_MAX.
  %psmaxp1 = getelementptr [21 x i8], [21 x i8]* @i64max_p1, i32 0, i32 0
  %maxp1 = call i64 @strtoll(i8* %psmaxp1, i8** null, i32 10)
  %ps1 = getelementptr i64, i64* %ps, i32 1
  store i64 %maxp1, i64* %ps1

; Do not fold a sequence consisting of just whitespace characters.
  %psws = getelementptr [7 x i8], [7 x i8]* @ws, i32 0, i32 0
  %nws = call i64 @strtoll(i8* %psws, i8** null, i32 10)
  %ps2 = getelementptr i64, i64* %ps, i32 2
  store i64 %nws, i64* %ps2

; Do not fold an empty sequence.
  %pswsp6 = getelementptr [7 x i8], [7 x i8]* @ws, i32 0, i32 6
  %nwsp6 = call i64 @strtoll(i8* %pswsp6, i8** null, i32 10)
  %ps3 = getelementptr i64, i64* %ps, i32 3
  store i64 %nwsp6, i64* %ps3

  ret void
}

@i_1_2_3_ = constant [9 x i8] c" 1 2\09\3\0a\00";

; Verify that strings of digits that are followed by whitespace are not
; folded (the whitespace could be interpreted in locales other than C
; as part of the leading digits, such as "123 456" is interepreted as
; 123456 in the French locale).

define void @call_strtol_trailing_space(i32* %ps) {
; CHECK-LABEL: @call_strtol_trailing_space(
; CHECK-NEXT:    [[N1:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([9 x i8], [9 x i8]* @i_1_2_3_, i64 0, i64 0), i8** null, i32 10)
; CHECK-NEXT:    [[PS1:%.*]] = getelementptr i32, i32* [[PS:%.*]], i64 1
; CHECK-NEXT:    store i32 [[N1]], i32* [[PS1]], align 4
; CHECK-NEXT:    [[N2:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([9 x i8], [9 x i8]* @i_1_2_3_, i64 0, i64 2), i8** null, i32 10)
; CHECK-NEXT:    [[PS2:%.*]] = getelementptr i32, i32* [[PS]], i64 2
; CHECK-NEXT:    store i32 [[N2]], i32* [[PS2]], align 4
; CHECK-NEXT:    [[N3:%.*]] = call i32 @strtol(i8* nocapture getelementptr inbounds ([9 x i8], [9 x i8]* @i_1_2_3_, i64 0, i64 4), i8** null, i32 10)
; CHECK-NEXT:    [[PS3:%.*]] = getelementptr i32, i32* [[PS]], i64 3
; CHECK-NEXT:    store i32 [[N3]], i32* [[PS3]], align 4
; CHECK-NEXT:    ret void
;
  %p1 = getelementptr [9 x i8], [9 x i8]* @i_1_2_3_, i32 0, i32 0
  %n1 = call i32 @strtol(i8* %p1, i8** null, i32 10)
  %ps1 = getelementptr i32, i32* %ps, i32 1
  store i32 %n1, i32* %ps1

  %p2 = getelementptr [9 x i8], [9 x i8]* @i_1_2_3_, i32 0, i32 2
  %n2 = call i32 @strtol(i8* %p2, i8** null, i32 10)
  %ps2 = getelementptr i32, i32* %ps, i32 2
  store i32 %n2, i32* %ps2

  %p3 = getelementptr [9 x i8], [9 x i8]* @i_1_2_3_, i32 0, i32 4
  %n3 = call i32 @strtol(i8* %p3, i8** null, i32 10)
  %ps3 = getelementptr i32, i32* %ps, i32 3
  store i32 %n3, i32* %ps3

  ret void
}
