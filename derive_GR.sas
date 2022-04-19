/******************************************************************************
Author, date: pmbrown, feb2020
Macro name:   derive_GR
Description:  derive the global rank endpoint on simulated data of up to 5 outcomes. 
              (developed in SAS 9.4)
Reference:    A Global Rank End Point for Clinical Trials in Acute Heart Failure, 
              Felker GM, Maisel AS, Circulation: Heart Failure, 2010
Validation    Susanne Bruun on 19th April at 15:49. 
******************************************************************************/

%macro derive_GR(indata=finalsamp,outdata=out_gr,outpval=pval_gr,
                 cut1=10,cut2=100,cut3=10,cut4=28,cut5=, 
         /*cut: cut-offs defining failure (for dichotomous outcome cut-off = 1)*/
                 desc1=d,desc2=d,desc3=d,desc4=,desc5=,
        /*desc: indicate whether high values are good or bad, ie how to sort data*/
                 ord1=2,ord2=3,ord3=4,ord4=1,ord5=); /*mortality is first*/
           /*ord: indicates position of outcome in hierarchy*/     

********************************************************************;
***              Identify and rank failures                      ***;
********************************************************************;

%do h=1 %to &gnumvars; 

proc sort data=&indata out=out&h (keep=samp subjno trt var&h);
  %if &&desc&h=d %then %do;
    by samp descending var&h;
  %end;
  %else %do;
    by samp var&h; 
  %end;
run;

data out2&h;
  format param $20. valuen best.;
  retain rank 0; 
  set out&h; 
  %if &&desc&h=d %then %do; 
    by samp descending var&h;
    if var&h ge &&cut&h then fail=1; 
  %end;
  %else %do;
    by samp var&h; 
    if var&h le &&cut&h then fail=1; 
  %end;
  ord=&&ord&h; 
  param="&&gtype&h"||"&h";
  valuen=var&h;
  if first.samp then rank=0; 
  if first.var&h then rank=rank+1; 
  keep samp subjno trt ord rank param valuen fail; 
run; 

proc sort data=out2&h ;
  by samp subjno;
run;

%end;

********************************************************************;
***                   derive global rank                         ***;
********************************************************************;

*calculate temporary score to be ranked below;
data globrnk;
  set out21-out2&gnumvars;
  if fail=1 and ord ne &gnumvars then 
    globrnk_tmp=ord*1000+rank; /*1000=arbitrarily large number, must be >> n*/
  /*rank non-fails on last outcome*/
  else if ord=&gnumvars then globrnk_tmp=&gnumvars*1000+rank; 
  /*remove data not failing and not belonging to last endpt*/
  else delete; 
  label ord='Position in hierarchy' rank='Rank within outcome'
        param='Outcome ranked on' valuen='Value';
run;

proc sort data=globrnk;
  by samp subjno globrnk_tmp;
run;

proc sort data=globrnk out=globrnk2 nodupkey;
  by samp subjno ;
run;

proc sort data=globrnk2;
  by samp globrnk_tmp;
run;

*assign global ranks;
data globrnk3;
  retain globrnk 0 ;
  set globrnk2;
  by samp globrnk_tmp;
  if first.samp then globrnk=0;
  if first.globrnk_tmp then globrnk=globrnk+1;
  label globrnk='Global Rank';
run;

proc sort data=globrnk3;
  by samp subjno;
run;

*put variables in a suitable order ;
data &outdata;
  retain samp subjno trt param valuen rank globrnk ;
  set globrnk3;
  keep samp subjno trt param valuen rank globrnk ;
run;

********************************************************************;
***                   p-value                         ***;
********************************************************************;

ods output KruskalWallisTest=&outpval ;
proc npar1way data=&outdata wilcoxon plots=none ;
  ods select KruskalWallisTest ; 
  class trt;
  var globrnk ;
  by samp;
quit;

%mend derive_GR;


***end****************************************************************;


