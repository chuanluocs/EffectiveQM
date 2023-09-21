# EffectiveQM: Towards Effective Local Search for Qubit Mapping Problem

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
| `19X19`      | 19X19       |
| `SYCAMORE`   | Google Sycamore        |

## Example Command for Running *EffectiveQM*

```bash
./EffectiveQM 1 log.txt TOKYO 5 5 500 ./benchmarks/4gt13_92.qasm
```

Running this command will call the `EffectiveQM` program, using `0` as random seed and set the hyperparameters $\delta=5$, $\theta=5$, $\lambda=500$, map `./benchmarks/4gt13_92.qasm` quantum circuit on the `IBM Tokyo` device, and store the resulting in the `./log.txt` file.
## Implementation of *EffectiveQM*

The directory named `src/` includes the implementation of *EffectiveQM*. 

## Testing Benchmarks for Evaluating *EffectiveQM*

The directory named `benchmarks/` contains all testing benchmarks‘ qasm files. 

## Implementation of Main Competitors of *EffectiveQM*

As mentioned in the paper, the main competitors of EffectiveQM are ILS,SAD2 and SID3.Their codes can be obtained from the following link.

- *ILS*: https://github.com/joyofly/ILS-QuantumCircuitMapper

- *SAD2*: https://github.com/BensonZhou1991/circuittransform/

- *SID3*: https://github.com/BensonZhou1991/Qubit-Mapping-Subgraph-Isomorphismr

Among them, for the two stochastic algorithms ILS and SAD2, we have slightly modified the implementation of the random number generation part of them for generating reproducible experimental results.

## Experimental Results

The directory `result/` contains all the raw data tables after preliminary statistics. 

The directory is organized as follows:

- [original_data](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data): Provides the original experimental results of *EffectiveQM*, *ILS*, *SAD2*, *SID3* on all devices.
  - [EffectiveQM](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/EffectiveQM): the original experimental results of *EffectiveQM*.
  - [ILS](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/ILS): the original experimental results of *ILS*.
  - [SAD2](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/SAD2): the original experimental results of *SAD2*.
  - [SID3](https://github.com/chuanluocs/EffectiveQM/tree/main/result/original_data/SID3): the original experimental results of *SID3*.
- [excel](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel): Provides collated experimental data against tables in the paper.
  - [comp](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/comp): Comparison data between *EffectiveQM* and all competitors.
  - [cutoff](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/cutoff): The comparative data after taking strict time constraints (refer to the paper), that is the comparative data of the *EffectiveQM-short* version in the paper.
  - [abalation](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/abalation): Provides data for ablation analysis experiments.
    - [abalation_for_Potential_Guided](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/abalation/abalation_for_Potential_Guided): The experimental data of the ablation analysis of Potential-guided scoring Function, that is, the comparative data of the *Alt-2* version in the paper.
    - [abalation_for_State_Aware](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/abalation/abalation_for_State_Aware): The experimental data of the ablation analysis of State Aware, that is, the comparative data of the *Alt-1* version in the paper.
  - [hyper-para](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/hyper-para): Experimental data for hyperparameter analysis experiments.
    - [lambda](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/hyper-para/lambda): Analyzing the experimental data on the hyperparameters of $\lambda$.
    - [delta](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/hyper-para/delta): Analyzing the experimental data on the hyperparameters of $\delta$.
    - [theta](https://github.com/chuanluocs/EffectiveQM/tree/main/result/excel/hyper-para/theta)：Analyzing the experimental data on the hyperparameters of $\theta$.
