# How to read the example file 

In the same directory as this readme file, there is an `example.txt` file that records the running log of the Mode-aware search strategy for a `qft_10` circuit on the IBM Tokyo device when $\delta$ is set to 2. A part of this log is used to illustrate Figure 4 in Section IV-B of the paper. By reading the `example.txt` file, you can further understand the running process of the Mode-aware search strategy.

## Start part  
The beginning of the file is four independent lines, respectively representing the parameter settings of this run, the random seed, the number of logic qubits of the `qft_10` circuit, and the number of CNOT gates. Ignoring them will not affect the reader's understanding of the running process of the Mode-aware search strategy.
```
[INFO] parameter: IN_LIMIT: 3 MAX_LOOK_SIZE: 4 C_MAX: 500
Seed: 0
The number of logic qubits:10
The number of CNOT gates:90
```

## The main part
The main part of the file is divided into two parts: Exploration Mode and Exploitation Mode. `EffectiveQM` switches between Exploration Mode and Exploitation Mode. Since `EffectiveQM` always starts with Exploration Mode, the beginning part of this log is Exploration Mode.

### Exploration Mode
Exploration Mode starts with 
```
Local optima found,  budget -1
```
which means `EffectiveQM` consumes a budget and switches to Exploration Mode to generate a physical circuit.This may be because `nonUpdCnt` exceeds the limit, and the `MD_Search_strategy` function determines that a local optimum has occurred, so it switches to Exploration Mode. It may also be because this is the first time a physical circuit is generated, so it switches to Exploration Mode.

After that, the log of Exploration Mode is similar to
```
Exploration Mode begin
PC* size is 162
nonUpdCnt is 0
```
This represents the start of Exploration Mode, the size of `PC*` is 162, and `nonUpdCnt` is 0. The symbols in these logs have specific meanings, please refer to Section IV-B of the manuscript.

### Exploitation Mode
The log of Exploitation Mode is similar to
```
Exploitation Mode begin
PC* size is 141
nonUpdCnt is 0
```
This represents the start of Exploitation Mode, the size of `PC*` is 141, and `nonUpdCnt` is 0. Exploitation Mode may appear after Exploration Mode, or after Exploitation Mode, depending on whether `nonUpdCnt` exceeds the limit.

### local optima found
When `nonUpdCnt` exceeds the limit, the `MD_Search_strategy` function determines that a local optimum has occurred, and the log is recorded as above.
```
nonUpdCnt out of limit, local optima found
```

## End part
In the end part, the following log will appear
```
[INFO] check over. Circult has no fault.
Cnot gate count: 27
Size of the best physical circuit:117
Direction of best mapping: 1
[INFO] restart count: 5494
[INFO] total bridge decision count: 3059
[INFO] best decision bridge count: 0
Time cost :19.784 seconds.
```
This marks some statistical information. For readers, only the `Cnot gate count` line and the `Time cost` line need to be understood. For example, in this example, `EffectiveQM` inserted 27 auxiliary CNOT gates, and the total running time was 19.784 seconds.