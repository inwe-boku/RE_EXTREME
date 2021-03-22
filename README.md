# RE_EXTREME

These scripts provide:

- the simulation of power generation (Nuclear, Fossil, Hydro, Wind, SolarPV) in Sweden using GAMS
- the possibility of simulating a multitude of simulations automatically according to a configuration file
- the semi-automatic analysis of results from the simulations

## Scientific Articles

Holtinger, S; Mikovits, C; Schmidt, J; Baumgartner, J; Arheimer, B; Lindstrom, G; Wetterlund, E.
(2019): The impact of climatic extreme events on the feasibility of fully renewable power systems: A case study for Sweden
ENERGY. 2019; 178: 695-713. http://dx.doi.org/10.1016/j.energy.2019.04.128

Mikovits, C; Wetterlund, E; Wehrle, S; Baumgartner, J; Schmidt, J.
(2021): Stronger together: Multi-annual variability of hydrogen production supported by wind power in Sweden
APPL ENERG. 2021; 282, 116082 http://dx.doi.org/10.1016/j.apenergy.2020.116082

## Result Downloads

Results can be found here: https://zenodo.org/record/3712940

## Accknowledgements

The Swedish Research Council Formas (dnr. 2016–20118) and Bio4Energy, Sweden financially supported this work. We also gratefully acknowledge support from the European Research Council (“reFUEL” ERC2017-STG 758149) and by CLIM2POWER, Sweden. Project CLIM2POWER is part of ERA4CS, an ERA-NET initiated by JPI Climate and funded by FORMAS (SE), BMBF (DE), BMWFW (AT), FCT (PT), EPA (IE), ANR (FR) with co-funding by the European Union (Grant 690462). We would also like to thank SMHI, the Swedish Meteorological and Hydrological Institute, for providing the S-Hype data.

## Headstart

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
