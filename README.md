# RE_EXTREME

This code accompanies the paper "Trade-offs between hydrogen production and hydropower flexibility in a Swedish power system with high shares of renewables" currently available as preprint at ![engrxiv](https://engrxiv.org/kgxr7/).

To run a set of simulations automatically run the helper tool rungams as follows:
```
python3 rungams.py -s min_cost_elh2.gms -y csets.yml -p 28 --run --clean
```

rungams.py supports following options:

```
usage: rungams.py [-h] [-s SIM] [-y YML] [-p MAX_PROCS] [--clean] [--resume]
                  [--run]

runs gams processes in parallel

optional arguments:
  -h, --help            show this help message and exit
  -s SIM, --sim SIM     defines the yaml file
  -y YML, --yaml YML    defines the yaml file
  -p MAX_PROCS, --procs MAX_PROCS
                        how many processors should be used, default and
                        maximum is number of processors - 1
  --clean
  --resume
  --run
```
