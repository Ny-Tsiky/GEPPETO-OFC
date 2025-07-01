# GEPPETO-OFC
GEPPETO-OFC is a speech motor control model that simulates the biomechanics of tongue motion during speech production. It uses an Optimal Feedback Control strategy to control a reduced biomechanical model of the tongue, implemented in MATLAB.

## Requirements
At this stage, GEPPETO-OFC has specific platform and software dependencies:

### Platform
✅ Supported: macOS on Apple Silicon (M1, M2, M3 chips)

❌ Not yet supported: macOS on Intel, Windows, Linux
(Instructions for these platforms will be provided in future updates.)

### MATLAB Version
MATLAB R2021a or later is required.

### MATLAB Toolboxes
The following MATLAB toolboxes must be installed:

Curve Fitting Toolbox

Deep Learning Toolbox

MATLAB Coder

Optimization Toolbox

Parallel Computing Toolbox

Signal Processing Toolbox

Statistics and Machine Learning Toolbox

Image Processing Toolbox

### External Dependencies
A working C# compiler (for macOS) is required.

On macOS, you can use the Mono C# compiler or install .NET via Homebrew (brew install dotnet-sdk).

## Installation
### Clone this repository:

bash
Copy
Edit
git clone https://github.com/yourusername/geppeto-ofc.git
cd geppeto-ofc

### Open the project in MATLAB R2021a or newer on an Apple Silicon Mac.

### Make sure all required toolboxes are installed (see list above).

### Ensure your system has a working C# compiler. For example:

bash
Copy
Edit
dotnet --version

### Run the main simulation script:

matlab
Copy
Edit
run_main_simulation.m

## Notes
This is an early version of the code intended for researchers or developers familiar with MATLAB.

Compatibility with other platforms (macOS Intel, Windows, Linux) is planned and will require different compilation instructions.

## License
This project is licensed under the MIT License.
