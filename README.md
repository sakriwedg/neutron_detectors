# neutron_detectors

`neutron_detectors` is a Python toolkit for reading neutron detector data files and generating quick reports on detector performance.

Its main objective is to provide a unified analysis framework for multiple neutron detector architectures, making comparisons between detectors and across different measurement campaigns straightforward and reproducible.

## Features

* Unified interface for analysing different neutron detectors.
* Automatic generation of detector performance reports.
* Extensible detector definitions.
* Reusable analysis functions compatible with all supported detectors.

## Project Structure

```text
neutron_detectors/
├── data/            # Input data files
├── detectors/       # Detector definitions
├── docs/            # Documentation
├── reports/         # Generated reports
├── scripts/         # Analysis scripts
├── src/             # Generic analysis functions
├── tests/           # Unit tests
└── requirements.txt # Python dependencies
```

## Installation

It is recommended to work inside a Python virtual environment.

Install the required libraries with:

```bash
python -m pip install -r requirements.txt
```

## Adding Support for a New Detector

To add a new detector architecture:

1. Copy `detectors/ILL.py`.
2. Rename the file according to your detector.
3. Adapt the implementation to match your detector's data format and characteristics.

All functions in `src/` are designed to work with any detector that follows this interface.

## Running the Analysis

Analysis scripts are located in the `scripts/` directory.

Execute the appropriate script to process your data and generate reports. Generated reports will be written to:

```text
reports/
```

## Goal

The long-term goal of this project is to establish a common methodology for analysing neutron detector data, enabling consistent performance evaluation and simplifying comparisons between different detector technologies and over time.
