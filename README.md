# owards Effective Local Search for Qubit Mapping Problem

*EffectiveQM* is an effective local search algorithm for solving qubit mapping problems. This repository contains the entire code implementation of *EffectiveQM*, the device coupling graph used in the experiment, the quantum circuits, raw data of experiments, and data tables after preliminary statistics.

## Instructions for Building *EffectiveQM*

```bash
sh build.sh	
```

By executing the script file, the user can build an executable of *EffectiveQM* in the repository root directory. Please note that this script should be run on 64-bit GNU/Linux OS. *EffectiveQM* is cross-platform, for other OS, please follow the *Makefile* in the *src* directory to build it by yourself.

## Instructions for Running *EffectiveQM*

The command to run EffectiveQM once is as follows:

```bash
./EffectiveQM [RANDOM_SEED] [RESULT_FILE_PATH] [DEVICE_NAME] [DELTA] [THETA] [LAMBDA] [QASM_FILE_PATH]
```

DEVICE_NAME is the parameter used to specify the quantum device, and its enumerated values are listed in the following table:
| value      | device in paper |
| ----------- | ----------- |
| `TOKYO`      | IBM Tokyo       |
| `Q16`   | Q16        |
| `GUADALUPE`      | IBM Guadalupe       |
| `ROCHESTER`   | IBM Rochester        |
| `torino`      | IBM Torino       |
| `SYCAMORE`   | Google Sycamore        |

## Example Command for Running *EffectiveQM*

```bash
./EffectiveQM 1 log.txt TOKYO 5 5 500 ./benchmarks/4gt13_92.qasm
```

Running this command will call the `EffectiveQM` program, using `0` as random seed and set the hyperparameters $\delta=5$, $\theta=5$, $\lambda=500$, map `./benchmarks/4gt13_92.qasm` quantum circuit on the `IBM Tokyo` device, and store the resulting in the `./log.txt` file.
## Implementation of *EffectiveQM*

The directory named `src/` includes the implementation of *EffectiveQM*. 

## Testing Benchmarks for Evaluating *EffectiveQM*

The directory named `benchmarks/` contains all testing benchmarks' qasm files. 

## Implementation of Main Competitors of *EffectiveQM*

As mentioned in the paper, the main competitors of EffectiveQM are ILS, SAHS, FiDLS, Qiskit and Tket. Their codes can be obtained from the following link.

- *ILS*: https://github.com/joyofly/ILS-QuantumCircuitMapper

- *SAHS*: https://github.com/BensonZhou1991/circuittransform/

- *FiDLS*: https://github.com/BensonZhou1991/Qubit-Mapping-Subgraph-Isomorphismr

- *Qiskit*: https://github.com/Qiskit/qiskit-terra

- *Tket*: https://github.com/CQCL/tket

Among them, for the three stochastic algorithms ILS, SAHS and Qiskit, we have slightly modified the implementation of the random number generation part of them for generating reproducible experimental results.

## Experimental Results

The directory `result/` contains all the raw data tables after preliminary statistics. 

The directory is organized as follows:

- [original_data](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data): Provides the original experimental results of *EffectiveQM*, *ILS*, *SAHS*, *FiDLS*, *Qiskit* and *Tket* on all devices.
  - [RevLib](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/RevLib): the original experimental results on benchmark RevLib.
  - [QV](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/QV): the original experimental results of benchmark QV.
  - [QUEKNO](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/QUEKNO): the original experimental results of benchmark QUEKNO.
  - [RW](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/RW): the original experimental results of benchmark RW.
- [table](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table): Provides collated experimental data against tables in the paper.
  - [comp](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/comp): Comparison data between *EffectiveQM* and all competitors.
  - [abalation](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/abalation): Provides data for ablation analysis experiments.
    - [abalation_for_Potential_Guided](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/abalation/abalation_for_Potential_Guided): The experimental data of the ablation analysis of Potential-guided scoring Function, that is, the comparative data of the *Alt-2* version in the paper.
    - [abalation_for_Mode_Aware](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/abalation/abalation_for_Mode_Aware): The experimental data of the ablation analysis of Mode Aware, that is, the comparative data of the *Alt-1* version in the paper.
  - [hyper-para](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/hyper-para): Experimental data for hyperparameter analysis experiments.
    - [lambda](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/hyper-para/lambda): Analyzing the experimental data on the hyperparameters of $\lambda$.
    - [delta](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/hyper-para/delta): Analyzing the experimental data on the hyperparameters of $\delta$.
    - [theta](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/hyper-para/theta)：Analyzing the experimental data on the hyperparameters of $\theta$.
  - [discussion](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/discussion): Experimental data for discussion experiments in paper.
    - [discussion_on_bridge](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/discussion/discussion_on_bridge): Experimental data for discussion of Bridge effectiveness.
    - [discussion_on_more_benchmark](https://github.com/chuanluocs/EffectiveQM/tree/main/result/table/discussion/discussion_on_more_benchmark): Experimental data discussed for more benchmark sets (i.e., RW, QUEKNO, and QV).



## The running examples mentioned in Section IV-B
The directory `example` contains examples of mode-aware search strategies from Section IV-B of the manuscript.

The directory is organized as follows:
- [example.txt](https://github.com/chuanluocs/EffectiveQM/tree/main/example/example.txt): Provide the complete examples mentioned in Section IV-B of the manuscript.
- [readme.md](https://github.com/chuanluocs/EffectiveQM/tree/main/example/readme.md): Provide a detailed readme document to explain how to interpret the examples.
