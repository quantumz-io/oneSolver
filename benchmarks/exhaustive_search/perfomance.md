# Performance tests of exhaustive search

The performance test were running on the `Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz` processor.
The results were reported by using `time` linux program for measuring the execution time of each model.

Obtained results are:
|Test model    | Real [s]| User [s]|Sys [s]|
|------------- |---------|---------|-------|
|simple.qubo   | 0.506   | 0.446   | 0.079 |
|test1.qubo    | 0.514   | 0.813   | 0.044 |
|test2.qubo    | 0.565   | 0.949   | 0.048 |
|csp5.qubo     | 0.555   | 0.674   | 0.088 |
|csp7.qubo     | 0.558   | 0.697   | 0.061 |
|csp13.qubo    | 0.569   | 0.763   | 0.081 |
|dwave_doc.qubo| 0.010   | 0.000   | 0.010 |


Where the measured times are following meaning:

**Real** is wall clock time - time from start to finish of the call.
This is all elapsed time including time slices used by other processes
and time the process spends blocked (for example if it is waiting for
I/O to complete).

**User** is the amount of CPU time spent in user-mode code (outside the
kernel) within the process. This is only actual CPU time used in
executing the process. Other processes and time the process spends
blocked do not count towards this figure.

**Sys** is the amount of CPU time spent in the kernel within the
process. This means executing CPU time spent in system calls within the
kernel, as opposed to library code, which is still running in
user-space. Like 'user', this is only CPU time used by the process. See
below for a brief description of kernel mode (also known as 'supervisor'
mode) and the system call mechanism.