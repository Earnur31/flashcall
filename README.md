# BigWig Analysis Script

This repository contains a Python script for analyzing BigWig files, calculating statistics, and identifying enriched intervals.

## Requirements

- Python 3.x
- pyBigWig
- pybedtools
- scipy
- numpy

## Installation

Install the required packages using pip:

```sh
pip install pyBigWig pybedtools scipy numpy


USAGE
python script.py -b /path/to/bigwigfile.bw -o /path/to/outputfile.bed -g human -p 0.05 -t 1.0 -d 50 -a 5.0
