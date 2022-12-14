import lldb
from intelpt_testcase import *
from lldbsuite.test.lldbtest import *
from lldbsuite.test import lldbutil
from lldbsuite.test.decorators import *

class TestTraceEvents(TraceIntelPTTestCaseBase):

    @testSBAPIAndCommands
    def testCPUEvents(self):
      trace_description_file_path = os.path.join(self.getSourceDir(), "intelpt-multi-core-trace", "trace_missing_threads.json")
      self.traceLoad(traceDescriptionFilePath=trace_description_file_path, substrs=["intel-pt"])

      self.expect("thread trace dump instructions 3 -e --forward -c 5",
        substrs=['''thread #3: tid = 3497496
    0: (event) HW clock tick [40450075477621505]
    1: (event) CPU core changed [new CPU=51]
    2: (event) HW clock tick [40450075477657246]
  m.out`foo() + 65 at multi_thread.cpp:12:21'''])

      self.expect("thread trace dump instructions 3 -e --forward -c 5 -J",
        substrs=['''{
    "id": 0,
    "event": "HW clock tick",
    "hwClock": 40450075477621505
  },
  {
    "id": 1,
    "event": "CPU core changed",
    "cpuId": 51
  }'''])

    @testSBAPIAndCommands
    def testPauseEvents(self):
      '''
        Everytime the target stops running on the CPU, a 'disabled' event will
        be emitted, which is represented by the TraceCursor API as a 'paused'
        event.
      '''
      self.expect("target create " +
            os.path.join(self.getSourceDir(), "intelpt-trace-multi-file", "a.out"))
      self.expect("b 12")
      self.expect("r")
      self.traceStartThread()
      self.expect("n")
      self.expect("n")
      self.expect("si")
      self.expect("si")
      self.expect("si")
      # We ensure that the paused events are printed correctly forward
      self.expect("thread trace dump instructions -e -f",
        patterns=[f'''thread #1: tid = .*
  a.out`main \+ 23 at main.cpp:12
    0: {ADDRESS_REGEX}    movl .*
    1: \(event\) software disabled tracing
    2: {ADDRESS_REGEX}    addl .*
    3: {ADDRESS_REGEX}    movl .*
    4: \(event\) software disabled tracing
  a.out`main \+ 34 \[inlined\] inline_function\(\) at main.cpp:4
    5: {ADDRESS_REGEX}    movl .*
  a.out`main \+ 41 \[inlined\] inline_function\(\) \+ 7 at main.cpp:5
    6: {ADDRESS_REGEX}    movl .*
    7: {ADDRESS_REGEX}    addl .*
    8: {ADDRESS_REGEX}    movl .*
  a.out`main \+ 52 \[inlined\] inline_function\(\) \+ 18 at main.cpp:6
    9: {ADDRESS_REGEX}    movl .*
  a.out`main \+ 55 at main.cpp:14
    10: {ADDRESS_REGEX}    movl .*
    11: {ADDRESS_REGEX}    addl .*
    12: {ADDRESS_REGEX}    movl .*
    13: \(event\) software disabled tracing
  a.out`main \+ 63 at main.cpp:16
    14: {ADDRESS_REGEX}    callq  .* ; symbol stub for: foo\(\)
    15: \(event\) software disabled tracing
  a.out`symbol stub for: foo\(\)
    16: {ADDRESS_REGEX}    jmpq'''])

      # We ensure that the paused events are printed correctly backward
      self.expect("thread trace dump instructions -e --id 16",
        patterns=[f'''thread #1: tid = .*
  a.out`symbol stub for: foo\(\)
    16: {ADDRESS_REGEX}    jmpq .*
    15: \(event\) software disabled tracing
  a.out`main \+ 63 at main.cpp:16
    14: {ADDRESS_REGEX}    callq  .* ; symbol stub for: foo\(\)
    13: \(event\) software disabled tracing
  a.out`main \+ 60 at main.cpp:14
    12: {ADDRESS_REGEX}    movl .*
    11: {ADDRESS_REGEX}    addl .*
    10: {ADDRESS_REGEX}    movl .*
  a.out`main \+ 52 \[inlined\] inline_function\(\) \+ 18 at main.cpp:6
    9: {ADDRESS_REGEX}    movl .*
  a.out`main \+ 49 \[inlined\] inline_function\(\) \+ 15 at main.cpp:5
    8: {ADDRESS_REGEX}    movl .*
    7: {ADDRESS_REGEX}    addl .*
    6: {ADDRESS_REGEX}    movl .*
  a.out`main \+ 34 \[inlined\] inline_function\(\) at main.cpp:4
    5: {ADDRESS_REGEX}    movl .*
    4: \(event\) software disabled tracing
  a.out`main \+ 31 at main.cpp:12
    3: {ADDRESS_REGEX}    movl .*
    2: {ADDRESS_REGEX}    addl .*
    1: \(event\) software disabled tracing
    0: {ADDRESS_REGEX}    movl .*'''])
