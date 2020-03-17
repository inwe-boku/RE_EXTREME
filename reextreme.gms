$title ReExtreme hydro-thermal distribution model

*******************
** Base settings **
*******************

$offlisting
$offsymxref offsymlist

options
  limrow  = 0,
  limcol  = 0,
  threads = 1
;

****************************
** Commandline parameters **
****************************

$if not set WF_HYDRO $set WF_HYDRO    0.1
$if not set WF_LOSSL $set WF_LOSSL 1500
$if not set WF_EXPRT $set WF_EXPRT 1000
$if not set WF_ELH   $set WF_ELH    500
$if not set P_MINHY  $set P_MINHY    62
$if not set P_HY0    $set P_HY0      62
$if not set CP_ELH   $set CP_ELH      0
$if not set PR_ELH   $set PR_ELH      0
$if not set HY_RMP   $set HY_RMP      5

***************************************
** Definition of sets and PARAMETERs **
***************************************

SETS
  scenario    /scen1*scen4/
  year        /y1986*y2014/
  month       /m1*m12/
  hour        /h0001*h8784/
  month_hour(month,hour)
              /m1  .(h0001*h0744)
               m2  .(h0745*h1416)
               m3  .(h1417*h2160)
               m4  .(h2161*h2880)
               m5  .(h2881*h3624)
               m6  .(h3625*h4344)
               m7  .(h4345*h5088)
               m8  .(h5089*h5832)
               m9  .(h5833*h6552)
               m10 .(h6553*h7296)
               m11 .(h7297*h8016)
               m12 .(h8017*h8784) /
   type       
   fuel
;

* INFLOW hydro in m3 MW
PARAMETER TS_INFLOW(year, hour)/
$include ../Data/gams/shype_unlimited_full.csv
/;

* bias correction factor reduce potential production by 10%
TS_INFLOW(year, hour) = TS_INFLOW(year, hour) / 1.1;

* NUCLEAR production per szenario 
PARAMETER PROD_NUCLEAR(scenario, year, hour)/
$include ../Data/gams/ts_nuclear_annually.csv
/;

* THERMAL production
PARAMETER PROD_THERMAL(type<, fuel<)/
$include ../Data/gams/ts_thermal_annually.csv
/;

PARAMETER PROD_THERMAL_TYPE;
PROD_THERMAL_TYPE(type) = sum(fuel, PROD_THERMAL(type, fuel));

* THERMAL costs, from -128 / 24 to 100
PARAMETER COST_THERMAL(type)/
$include ../Data/gams/cost_thermal_co2.csv
/;

* THERMAL minimum production
PARAMETER THERMAL_MIN_MONTH(month, type)/
$include ../Data/gams/thermal_min_monthtype.csv
/;

* THERMAL maximum production
PARAMETER THERMAL_MAX_MONTH(month, type)/
$include ../Data/gams/thermal_max_monthtype.csv
/;


PARAMETERS
  NET_LOAD(hour)    hourly load
  INFLOW(hour)      hourly inflow to reservoir
  IRE(hour)         hourly renewable generation
  NUCLEAR(hour)     hourly nuclear generation
  TOTALSTORAGE      total reservoir capacity available
  MAXSTORAGE        maximum reservoir capacity
  MINSTORAGE        minimum reservoir level
  STARTSTORAGE      reservoir level at the start of the year
  MINENDSTORAGE     minimum reservoir level at the end of the year
  MAXRAMP           maximum hourly ramps of hydro power generation
  MAXRAMPTHERMAL    maximum hourly ramps of thermal power generation
  MAXRAMPELH2       maximum hourly ramps of hydrolysis
  MINFLOW           minimum hourly flow
  HYDROCAPACITY     installed hydro power capacity
  CF_HYDRO          maximum capacity factor for hydro power
  CF_THERMAL        maximum capacity factor for thermal power
  WF_HYDRO          cost hydro power generation
  WF_LOSSL          cost loss of load
  WF_EXPRT          cost for exports
  WF_ELH            cost for electrolysis
  P_MINHY           percentage minimal reservoir level
  P_HY0             percentage start value first year reservoir level
  CP_ELH            electrolysis capacity
  PR_ELH            electrolysis ramping percentage (of total capacity per hour)
  HY_RMP             scenario number for hydro ramping and minimum flow
;

* set values
  WF_HYDRO = %WF_HYDRO%;
  WF_LOSSL = %WF_LOSSL%;
  WF_EXPRT = %WF_EXPRT%;
  WF_ELH   = %WF_ELH%;
  P_MINHY  = %P_MINHY%;
  P_HY0    = %P_HY0%;
  CP_ELH   = %CP_ELH%;
  PR_ELH   = %PR_ELH%;
  HY_RMP    = %HY_RMP%;
  
  TOTALSTORAGE    = 33.7E6;
  MAXSTORAGE      = TOTALSTORAGE * 0.98;
  MINSTORAGE      = TOTALSTORAGE * 0.05;
  MINENDSTORAGE   = MAXSTORAGE * P_MINHY / 100;
  MAXRAMP         = 4000;
  MAXRAMPTHERMAL  = 1500;
* ramping of hydrolysis: xx % of total capacity
  MAXRAMPELH2     = PR_ELH * CP_ELH / 100;
* http://www.svenskenergi.se/Global/Statistik/El%C3%A5ret/QUICK%20FACTS%20ABOUT%20SWEDEN%20AND%20ENERGY_2014.pdf
  HYDROCAPACITY   = 16155;
  CF_HYDRO        = 0.85;
  CF_THERMAL      = 0.98;
  HYDROCAPACITY   = HYDROCAPACITY * CF_HYDRO;
  PROD_THERMAL_TYPE(type) =   PROD_THERMAL_TYPE(type) * CF_THERMAL
  
* net load per szenario in MW
PARAMETER TS_NET_LOAD(scenario, year, hour)/
$include ../Data/Biofuel/loadelh2/net_load_full-%CP_ELH%.csv
/;

* IRE production per szenario MW
PARAMETER PROD_IRE(scenario, year, hour)/
$include ../Data/Biofuel/windelh2/ts_ire_annually-%CP_ELH%.csv
/;

* HYDRO minimum flow, monthly values
PARAMETER HYDRO_MIF_MONTH(month)/
$include ../Data/gams/Hydro/MIF%HY_RMP%.csv
/;

* HYDRO minimum flow remapped to hours
PARAMETER HYDRO_MIF_HOUR(hour);
HYDRO_MIF_HOUR(hour) = sum(month $month_hour(month, hour), HYDRO_MIF_MONTH(month));

* HYDRO maximum ramping per month
PARAMETER HYDRO_MRR_MONTH(month)/
$include ../Data/gams/Hydro/MRR%HY_RMP%.csv
/;

* HYDRO maximum ramping remapped to hours
PARAMETER HYDRO_MRR_HOUR(hour);
HYDRO_MRR_HOUR(hour) = sum(month $month_hour(month, hour), HYDRO_MRR_MONTH(month));

* THERMAL minimum generation per typed remapped to hourly values
PARAMETER THERMAL_MIN_HOUR(hour, type);
THERMAL_MIN_HOUR(hour, type) = sum(month $month_hour(month, hour), (THERMAL_MIN_MONTH(month, type) * PROD_THERMAL_TYPE(type) / 100));

* THERMAL maximum generation per typed remapped to hourly values
PARAMETER THERMAL_MAX_HOUR(hour, type);
THERMAL_MAX_HOUR(hour, type) = sum(month $month_hour(month, hour), ((100-THERMAL_MAX_MONTH(month, type)) * PROD_THERMAL_TYPE(type) / 100));

display HYDRO_MIF_HOUR;
display THERMAL_MAX_HOUR;
*******************************************************************************
***                    Variables and equations
*******************************************************************************

positive VARIABLES
  x_hydro_gen(hour)
  x_thermal_gen(hour, type)
  x_curtailment(hour)
  x_loss_of_load(hour)
  x_export(hour)
  x_reslevel(hour)
  x_spillhydro(hour)
  x_elh2(hour)
;

free VARIABLES
  total_cost;
  x_reslevel.lo(hour)    = MINSTORAGE;
  x_reslevel.up(hour)    = MAXSTORAGE;
  x_hydro_gen.up(hour)   = HYDROCAPACITY;
  x_spillhydro.up(hour)  = 50000;
  
EQUATIONS
  objective
  equ_thermal_gen
  equ_load_balance
  reslevel
  reslevelend
  equ_minflow
  equ_maxFlow
  equ_hydrochangeup
  equ_hydrochangedown
  equ_minthermalprod
  equ_maxthermalprod
  maxcurtail
  rampup
  rampdown
  equ_elh2
  rampupthermal
  rampdownthermal
  rampupelh2
  rampdownelh2
;

* minimize balancing costs
objective..
              total_cost
              =E=
              sum(hour, (x_spillhydro(hour) * WF_HYDRO/2)
                + (x_hydro_gen(hour) * WF_HYDRO)
                + (sum(type, x_thermal_gen(hour, type) * COST_THERMAL(type))) 
                + (x_loss_of_load(hour) * WF_LOSSL)
                + (x_export(hour) * WF_EXPRT)
                - (x_elh2(hour) * WF_ELH) )
              ;
              
equ_thermal_gen(hour, type)..
             x_thermal_gen(hour, type)
             =L=
             PROD_THERMAL_TYPE(type)
             ;
             
equ_elh2(hour)..
            x_elh2(hour)
            =L=
            CP_ELH
            ;

equ_load_balance(hour)..
              NET_LOAD(hour)
              =E=
              x_hydro_gen(hour)
              + sum(type, x_thermal_gen(hour, type))
              + x_loss_of_load(hour)
              - x_curtailment(hour)
              - x_export(hour)
              - x_elh2(hour)
              ;

* control that reservoir level equals level of past period + INFLOW - hydroprod
reslevel(hour)$(ord(hour)>1)..
              x_reslevel(hour)
              =E=
              x_reslevel(hour-1)
              + INFLOW(hour)
              - x_hydro_gen(hour)
              - x_spillhydro(hour)
              ;

* Reservoir Endlevel is equal or higher than the defined MINENDSTORAGE
reslevelend..
              x_reslevel("H8760")
              =G=
              MINENDSTORAGE
              ;

*  minimum flow constraint (historic minimum production)
equ_minflow(hour)..
              x_hydro_gen(hour) + x_spillhydro(hour)
              =G=
              HYDRO_MIF_HOUR(hour)
              ;

*  maximum flow constraint (historic maximum production)
equ_maxFlow(hour)..
              x_hydro_gen(hour) + x_spillhydro(hour)
              =L=
              HYDROCAPACITY
              ;

*  maximum flow rate change by month
equ_hydrochangeup(hour)$(ord(hour) > 1)..
              x_hydro_gen(hour) + x_spillhydro(hour) -
              x_hydro_gen(hour-1) - x_spillhydro(hour-1)
              =L=
              HYDRO_MRR_HOUR(hour)
              ;

equ_hydrochangedown(hour)$(ord(hour) > 1)..
              x_hydro_gen(hour-1) + x_spillhydro(hour-1) -
              x_hydro_gen(hour) - x_spillhydro(hour)
              =L=
              HYDRO_MRR_HOUR(hour)
              ;

*  hourly minimum thermal generation for each month
equ_minthermalprod(type, hour)..
              x_thermal_gen(hour, type)
              =G=
              THERMAL_MIN_HOUR(hour, type)
              ;

*  hourly maximum thermal generation for each month
equ_maxthermalprod(type, hour)..
              x_thermal_gen(hour, type)
              =L=
              THERMAL_MAX_HOUR(hour, type)
              ;

* maximum RE curtailment
maxcurtail(hour)..
              x_curtailment(hour)
              =L=
              IRE(hour) + NUCLEAR(hour)
              ;

* ramping restrictions hydro
rampup(hour)$(ord(hour) > 1)..
              x_hydro_gen(hour)  - x_hydro_gen(hour-1)
              =L=
              MAXRAMP
              ;

rampdown(hour)$(ord(hour) > 1)..
              x_hydro_gen(hour-1)  - x_hydro_gen(hour)
              =L=
              MAXRAMP
              ;

* ramping restrictions thermal
rampupthermal(hour)$(ord(hour) > 1)..
              sum(type, x_thermal_gen(hour, type)) - sum(type, (x_thermal_gen(hour - 1, type)))
              =L=
              MAXRAMPTHERMAL
              ;

rampdownthermal(hour)$(ord(hour) > 1)..
              sum(type, (x_thermal_gen(hour - 1, type))) - sum(type, (x_thermal_gen(hour, type)))
              =L=
              MAXRAMPTHERMAL
              ;

rampupelh2(hour)$(ord(hour) > 1)..
              x_elh2(hour)  - x_elh2(hour-1)
              =L=
              MAXRAMPELH2
              ;

rampdownelh2(hour)$(ord(hour) > 1)..
              x_elh2(hour-1)  - x_elh2(hour)
              =L=
              MAXRAMPELH2
              ;


MODEL OptNetLoad /all/;

FILE outfile / result.csv /;
PUT outfile;

outfile.lw=0;
outfile.nw=0;
outfile.pw=1000;
outfile.pc=5;
outfile.nd=3;

* CSV Header
PUT "model_stat" "solver_stat" "res_used"
PUT  "scenario" "year" "hour"
PUT "total_cost" "NET_LOAD" "VREgen" "NUCLEAR_gen" "hydro_generation" "thermal_generation"  "curtailment"  "loss_of_load" "res_level" "hydro_spill"
LOOP(type,
  PUT type.tl 
);
PUT "elh2_generation"
PUT "export"
PUT /

LOOP(scenario$(ord(scenario) = 3),
  STARTSTORAGE = MAXSTORAGE * P_HY0 / 100;
  
  LOOP(year$(ord(year) < 99),

    STARTSTORAGE = MAXSTORAGE * P_HY0 / 100;
    x_reslevel.fx("H0001") = STARTSTORAGE;

    NET_LOAD(hour)              = 0;
    INFLOW(hour)                = 0;
    x_hydro_gen.L(hour)         = 0;
    x_reslevel.L(hour)          = 0;
    LOOP(type,
      x_thermal_gen.L(hour, type) = 0
    );
    x_spillhydro.L(hour)        = 0;
    x_curtailment.L(hour)       = 0;
    NET_LOAD(hour)              = TS_NET_LOAD(scenario, year, hour);
    INFLOW(hour)                = TS_INFLOW(year, hour);
    IRE(hour)                   = PROD_IRE(scenario, year, hour);
    NUCLEAR(hour)               = PROD_NUCLEAR(scenario, year, hour);

* solve Model

    SOLVE OptNetLoad using lp minimizing total_cost;


    LOOP(hour$NET_LOAD(hour),
      PUT OptNetLoad.modelstat
      PUT OptNetLoad.solvestat
      PUT OptNetLoad.resusd
      PUT scenario.tl
      PUT year.tl
      PUT hour.tl
      PUT total_cost.L
      PUT NET_LOAD(hour)
      PUT IRE(hour)
      PUT NUCLEAR(hour)
      PUT x_hydro_gen.L(hour)
      PUT sum(type, x_thermal_gen.L(hour, type))
      PUT x_curtailment.L(hour)
      PUT x_loss_of_load.L(hour)
      PUT x_reslevel.L(hour)
      PUT x_spillhydro.L(hour)
      LOOP(type,
        PUT x_thermal_gen.L(hour, type)
      );
      PUT x_elh2.L(hour)
      PUT x_export.L(hour)
      PUT /
    );

*    STARTSTORAGE = x_reslevel.L("H8760");
  );
);

PUTCLOSE;
