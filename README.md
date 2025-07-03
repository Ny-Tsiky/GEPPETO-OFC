# GEPPETO-OFC
GEPPETO-OFC is a speech motor control model that simulates the biomechanics of tongue motion during speech production. It uses an Optimal Feedback Control strategy to control a reduced biomechanical model of the tongue, implemented as a compiled MEX function in MATLAB.

This repository includes a simulation script for producing the vowel-consonant-vowel sequence: /ə/ – /t/ – /e/.

## Requirements
At this stage, GEPPETO-OFC has specific platform and software dependencies:

### Platform
Supported: macOS on Apple Silicon (M1, M2, M3 chips)

### MATLAB Version
MATLAB R2021a or later is required.

### MATLAB Toolboxes
**The following MATLAB toolboxes must be installed:**

Curve Fitting Toolbox

Deep Learning Toolbox

MATLAB Coder

Optimization Toolbox

Parallel Computing Toolbox

Signal Processing Toolbox

Statistics and Machine Learning Toolbox

### External Dependencies
A working C# compiler is required.

## Installation
1. Clone this repository:
```

bash

git clone https://github.com/yourusername/geppeto-ofc.git
cd geppeto-ofc

```

2. Open the project in MATLAB.

3. Build Mex functions:
```
matlab

mex_lstm.m
```

4. Run the main simulation script:
```
matlab

run_main_simulation.m
```

5. Visualize the results:
```
matlab

generate_plot.m
```

## Notes
This is an early version of the code.

Compatibility with other platforms (macOS Intel, Windows, Linux) is planned and will require different compilation instructions.

## License
This project is licensed under the AGPL-3.0 license.
