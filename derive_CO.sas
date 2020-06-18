/******************************************************************************
Author, date: pmbrown, feb2020
Macro name:   derive_CO
Description:  implement decision rule as per CITRIS-ALI trial
Notes:        developed in SAS 9.4
Reference:    jama: Only if the smallest P value was less than .02, the
                  second smallest less than .03, and the largest less than .05
                  was simulation considered a success

******************************************************************************/

*set to missing if died;
data randsampmiss;
  set der.randsamp2;
  if var42 lt 28 then do; 
    var1_=.;
    var2_=.;
    var3_=.;
  end;
  else do;
    var1_=var1;
    var2_=var2;
    var3_=var3;
  end;
run;

ods output KruskalWallisTest=co_pval ;
proc npar1way data=randsampmiss wilcoxon plots=none ;
  ods select KruskalWallisTest ; 
  class trt;
  var var1 var2 var3
      var1_ var2_ var3_;
  by samp;
quit;

data co_pval2;
  set co_pval;
  if index(variable,'_')=0 then newvar=1;
  else newvar=2;
run;

*sort p-values from smallest to largest;
proc sort data=co_pval2;
  by newvar samp prob;
run;

*create var indicating sort order and use in subsequent proc transpose ;
data co_pval3;
  retain count;
  set co_pval2;
  by newvar samp prob;
  if first.samp then count=1;
  else count=count+1;
run;

*thus: pval1 is smallest etc;
proc transpose data=co_pval3 out=der.co_pval (drop=_name_ _label_) prefix=pval;
  var prob;
  id count;
  by newvar samp;
run;


***end***************************************************************************;
