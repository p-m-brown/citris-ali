/******************************************************************************
Author, date: pmbrown, feb2020
Program name: figure1.sas
Description:  create figure showing power estimates against assumed sample size
              (developed in SAS 9.4)

Reference:    jama 2019: Effect of Vitamin C Infusion on Organ Failure and Biomarkers
                    of Inflammation and Vascular Injury in Patients With Sepsis
                    and Severe Acute Respiratory Failure
Changes:      

******************************************************************************/


********************************************************************;
***                   preamble                            ***;
********************************************************************;

*options;
options nodate nonumber nofmterr nomprint ; 

*data libraries;
libname der "T:\alle\sepsisforskning\pmb\vitC\der";

*clear temp datasets;
proc datasets lib=work kill memtype=data;
run;
quit;

*define formats;
proc format;
  value $compos_f
    'CO'='Decision rule'
    'COmiss'='Decision rule (assuming missing data if patient dies)'
    'GR'='Global rank';
run;

*define macros;
%include 'T:\alle\sepsisforskning\pmb\vitC\NearestCorr.sas';
%include 'T:\alle\sepsisforskning\pmb\vitC\simul_data.sas';
%include 'T:\alle\sepsisforskning\pmb\vitC\derive_gr.sas';

*save log for records;
proc printto log='T:\alle\sepsisforskning\pmb\vitC\output\figure1.log' new ;
run;
proc printto print='T:\alle\sepsisforskning\pmb\vitC\output\figure1.lst' new;
run;

********************************************************************;
***           generate data                   ***;
********************************************************************;

%macro varyn();

%let total_t0 = %sysfunc(datetime()); 

data der.power;
  set _null_;
run;


%do varyn=50 %to 300 %by 10;   /*even numbers ie total n, n/2 per group*/  
                               /*CITRIS-ALI used n=170 thus vary around this*/

%iterat_simul(n_=&varyn,
              out=der.randsamp,
              out_iter=der.iterations); 


********************************************************************;
***       derive composite and decision rule               ***;
********************************************************************;

data der.randsamp2;
  set der.randsamp;
  if var1 lt 0 then var1=0;
  if var2 lt 0 then var2=0;
  if var3 lt 0 then var3=0;
run;

%derive_GR(indata=der.randsamp2,outdata=der.globrnk,outpval=der.gr_pval);
%include 'T:\alle\sepsisforskning\pmb\vitC\derive_co.sas';

********************************************************************;
***      estimate power              ***;
********************************************************************;

data der.pvals;
  format compos $compos_f.;
  set der.gr_pval (in=gr) der.co_pval (in=co);
  if gr then do;
    compos='GR';
    if prob lt 0.05 then sig=1;
    else sig=0;
  end;
/*decision rule: only if the smallest p-value was less than .02, the
                  second smallest less than .03, and the largest less than .05
                  was simulation considered a success*/
  else if co then do;
    if newvar=1 then compos='CO';
    else if newvar=2 then compos='COmiss';
    if pval1 lt 0.02 and pval2 lt 0.03 and pval3 lt 0.05 then sig=1;
    else sig=0;
  end;
  keep samp compos sig;
run;

proc sort data=der.pvals;
  by compos;
run;

ods output onewayfreqs=pctsig (where=(sig=1) keep=sig compos percent);
proc freq data=der.pvals;
  ods select onewayfreqs;
  tables sig;
  by compos;
run;

********************************************************************;
***      collect ests              ***;
********************************************************************;

%let total_speed=%sysfunc(datetime())-&total_t0;

data der.power;
  set der.power pctsig (in=a);
  if a then do;
    power=percent;
    totaln=&varyn;
    precision=&gprec; /*remove those which dont converge sufficiently*/
    total_speed=put((&total_speed)/60,5.2)||' mins';
  end;
  drop percent;
run;

%end; 
%mend varyn;

%varyn();

*reset log and output;
proc printto;
run;

********************************************************************;
***     create figure             ***;
********************************************************************;

ods listing gpath="T:\alle\sepsisforskning\pmb\vitC\output";
ods graphics on / imagefmt=tiff reset imagename="figure1" width=2500px 
height=2000px outputfmt=tiff maxobs=2074738 antialiasmax=9800 ;
ods html close;

title1 j=c h=3 "Figure 1. Estimated power by assumed total sample size";
title2;
title3;

footnote1 j=l h=1 "Decision rule: Holm-Bonferroni adjustment used in CITRIS-ALI.";
footnote2;
footnote3;
footnote4;

proc sgplot data=der.power;
  format compos $compos_f.;
  styleattrs datasymbols=(circlefilled trianglefilled squarefilled );
  series x=totaln y=power / group=compos lineattrs=(thickness=4);
  scatter x=totaln y=power / group=compos 
    filledoutlinedmarkers  markerattrs=(size=15) name='s';
  yaxis values=(0 to 80 by 20)
    label='Estimated statistical power (%)' grid gridattrs=(pattern=longdash thickness=2)
    offsetmax=0.05 offsetmin=0.05
    labelattrs=(color=black family=Arial Size=18 style=normal )
    valueattrs=(color=black family=Arial Size=15 style=normal ) ;
  xaxis values=(50 to 300 by 50)
    label='Total sample size (n)' offsetmax=0.05 offsetmin=0.05
    labelattrs=(color=black family=Arial Size=18 style=normal )
    valueattrs=(color=black family=Arial Size=15 style=normal ) ;
  keylegend 's' / title="" linelength=20 noborder autoitemsize
    autooutline 
    valueattrs=(color=black family=Arial Size=15 style=normal )
    titleattrs=(color=black family=Arial Size=18 style=normal )
    location=outside position=bottom;
run;

quit;
ods graphics off;



********************************************************************;
