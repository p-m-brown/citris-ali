/******************************************************************************
Author, date: Paul Brown, feb2020
Macro name:   simul_data, 
              iterat_simul (run simul_data iteratively in order to
              obtain the desired correlations between the outcomes), 
              nearestcorr
Description:  simulate random samples from a multivariate normal distribution,
              transform variables to required types. Create a dataset containing 
              num random samples indicated by the variable samp=1,2,...num.  
Reference:    
Notes:        _co=control, _ac=active
              co=control grp, the base response
              df=treatment difference, defines response for 
                              ac=active grp
              _sd: if var1_tp=norm or lognorm (sd for log(post/base))
              tte_fu: follow-up time for censoring of TTE endpts, corresponds
                               to the response estimates assumed for TTE outcomes
              assumed effects sizes: see email sent on 22feb2020
                rough sd obtained as: (q3-q1)/1.35, from fig2 jama
Validation:   
******************************************************************************/


********************************************************************;
***               define macro parameters                       ***;
********************************************************************;

/*!if update _df here also update in macro call below!*/

%macro simul_data(n=100,     /*total sample size ie n/2 in each grp*/
             num=1000,       /*no. of random samples*/
             numvars=4,      /*number of variables*/
             tte_fu=28,      /*days, corresponds to assumed response for var4*/
             var1_tp=NORM, var1_co=7, var1_df=-1, var1_sd=3.7, /*sofa score 96hrs */ 
             var2_tp=NORM, var2_co=40, var2_df=-5, var2_sd=45, /*crp (ug/mL) 168hrs*/
             var3_tp=NORM, var3_co=12, var3_df=-2, var3_sd=5, /*thrombomodulin (ng/mL) 168hrs*/
             var4_tp=SURV, var4_co=0.75, var4_df=0.02, var4_sd=, /*mortality at d28*/
             var5_tp=, var5_co=, var5_df=, var5_sd=,       
             corr12=0.2, corr13=0.2, corr14=0.05, corr15=0,           /*initial working correlations */
                         corr23=0.2, corr24=0.05, corr25=0,
                                     corr34=0.05, corr35=0,
                                                  corr45=0); 

%let gtype1=%upcase(&var1_tp) ;
%let gtype2=%upcase(&var2_tp) ;
%let gtype3=%upcase(&var3_tp) ;
%let gtype4=%upcase(&var4_tp) ;
%let gtype5=%upcase(&var5_tp) ;

********************************************************************;
***               convert to normal variates           ***;
********************************************************************;

data _null_;
  %do i=1 %to &numvars;
    %if &&gtype&i=SURV %then %do;
      surv_co=&&var&i._co;
      surv_ac=&&var&i._co+&&var&i._df;
      *--hazards;
      survh_co=-log(surv_co)/&tte_fu;
      survh_ac=-log(surv_ac)/&tte_fu;
      *--log hazards;
      survz_co=log(survh_co);
      survz_ac=log(survh_ac);
      *--number of deaths;
      survd_co=(&n/2)*(1-surv_co);
      survd_ac=(&n/2)*(1-surv_ac);
      *--variance;
      survv_co=1/survd_co; 
      survv_ac=1/survd_ac;
      call symput("var&i.v_co",trim(left(put(survv_co,best.))));
      call symput("var&i.v_ac",trim(left(put(survv_ac,best.))));
      call symput("var&i.z_co",trim(left(put(survz_co,best.))));
      call symput("var&i.z_ac",trim(left(put(survz_ac,best.))));
    %end;
    %else %if &&gtype&i=BINO %then %do;
      bino_co=&&var&i._co;
      bino_ac=&&var&i._co+&&var&i._df;
      *--log odds;
      binoz_co=log(bino_co/(1-bino_co));
      binoz_ac=log(bino_ac/(1-bino_ac));
      *--variance of log odds;
      binov_co=1/((&n/2)*bino_co*(1-bino_co)); 
      binov_ac=1/((&n/2)*bino_ac*(1-bino_ac));
      call symput("var&i.v_co",trim(left(put(binov_co,best.))));
      call symput("var&i.v_ac",trim(left(put(binov_ac,best.))));
      call symput("var&i.z_co",trim(left(put(binoz_co,best.))));
      call symput("var&i.z_ac",trim(left(put(binoz_ac,best.))));
    %end;
    %else %if &&gtype&i=NORM %then %do;
      normv_co=&&var&i._sd**2;
      normv_ac=&&var&i._sd**2;
      normz_ac=&&var&i._co+&&var&i._df;
      normz_co=&&var&i._co;
      call symput("var&i.v_co",trim(left(put(normv_co,best.))));
      call symput("var&i.v_ac",trim(left(put(normv_ac,best.))));
      call symput("var&i.z_co",trim(left(put(normz_co,best.))));
      call symput("var&i.z_ac",trim(left(put(normz_ac,best.))));
    %end;
    %else %if &&gtype&i=LOGN %then %do;
      *--log ratio (post baseline/baseline);
      lognz_ac=log(&&var&i._co+&&var&i._df);
      lognz_co=log(&&var&i._co);
      *--variance of log ratio;
      lognv_co=&&var&i._sd**2;
      lognv_ac=&&var&i._sd**2;
      call symput("var&i.v_co",trim(left(put(lognv_co,best.))));
      call symput("var&i.v_ac",trim(left(put(lognv_ac,best.))));
      call symput("var&i.z_co",trim(left(put(lognz_co,best.))));
      call symput("var&i.z_ac",trim(left(put(lognz_ac,best.))));
    %end;
  %end;
run;


********************************************************************;
***          simulate data from multivariate normal         ***;
********************************************************************;

proc iml;
%nearestcorr(); 
  *specify parameters of multivariate normal dist;
  mean_co=1:&numvars; varn_co=1:&numvars;
  mean_ac=1:&numvars; varn_ac=1:&numvars;
  %do i=1 %to &numvars;
    mean_co[&i]=&&var&i.z_co; 
    varn_co[&i]=&&var&i.v_co; 
    mean_ac[&i]=&&var&i.z_ac; 
    varn_ac[&i]=&&var&i.v_ac;
  %end;
  corr_tmp_={1 &corr12 &corr13 &corr14 &corr15, 
        &corr12 1 &corr23 &corr24 &corr25, 
        &corr13 &corr23 1 &corr34 &corr35,
        &corr14 &corr24 &corr34 1 &corr45,
        &corr15 &corr25 &corr35 &corr45 1}; 
  corr_tmp=corr_tmp_[1:&numvars, 1:&numvars]; /*in case numvar < 5*/
  eigval=eigval(corr_tmp); 
  if all(eigval>0) then corr=corr_tmp; /*if positive definite*/
  else corr=NearestCorr(corr_tmp); /*otherwise find nearest pos def matrix*/
  print corr; 
  covr_co=corr#sqrt(varn_co`*varn_co); 
  covr_ac=corr#sqrt(varn_ac`*varn_ac); 
  print covr_co; print covr_ac; 
  *obtain random samples (seed=dec2014);
  call randseed(1214);           
  co=randnormal((&n/2)*&num,mean_co,covr_co); 
  ac=randnormal((&n/2)*&num,mean_ac,covr_ac);
  *create sas datasets;
  samp=colvec(repeat(T(1:&num),1,&n/2)); 
  z=samp||co;
  create rand_co from z[c={"samp" "x1" "x2" "x3" "x4" "x5"}];
    append from z;
  close rand_co;
  z=samp||ac;
  create rand_ac from z[c={"samp" "x1" "x2" "x3" "x4" "x5"}];
    append from z;
  close rand_ac;
quit;

*combine data for the 2 grps;
data randztmp;
  format trt $10.;
  set rand_co (in=c) rand_ac (in=a);
  if a then trt='Active';
  else if c then trt='Control';
run;

*sort data;
proc sort data=randztmp;
  by samp trt;
run;

*finalise data, create patient number;
data randz;
  retain subjno;
  set randztmp;
  by samp trt;
  if first.samp then subjno=1;
  else subjno=subjno+1;
run;

********************************************************************;
***               transform from normal variates                 ***;
********************************************************************;

data randsamp;
  set randz;
  %do i=1 %to &numvars;
    %if &&gtype&i=SURV %then %do;
      var&i=round(-log(ranuni(1412))/exp(x&i),1)+1;  
      var&i.2=var&i; var&i.c=0;
      if var&i gt &tte_fu then do;
        var&i.2=&tte_fu; 
        var&i.c=1;
      end;
      label var&i="Survival endpoint &i"
        var&i.2="Survival endpoint &i with censoring" 
        var&i.c="Survival endpoint &i censoring indicator";
    %end;
    %else %if &&gtype&i=BINO %then %do;
      if trt='Active' then do;
        x&i._=(x&i-&&var&i.z_ac)/sqrt(&&var&i.v_ac);
        x&i.__=tinv(1-(&&var&i._co+&&var&i._df),&n/2-1); 
        if x&i._ gt x&i.__ then var&i=1; 
        else var&i=0;
      end;
      else if trt='Control' then do;
        x&i._=(x&i-&&var&i.z_co)/sqrt(&&var&i.v_co);
        x&i.__=tinv(1-&&var&i._co,&n/2-1);
        if x&i._ gt x&i.__ then var&i=1;
        else var&i=0;
      end;
      drop x&i._ x&i.__;
      label var&i="Dichotomous endpoint &i";
    %end;
    %else %if &&gtype&i=NORM %then %do;
      var&i=x&i;
      label var&i="Continuous endpoint &i";
    %end;
    %else %if &&gtype&i=LOGN %then %do;
      var&i=100*(exp(x&i)-1); 
      label var&i="Pct change endpoint &i";
    %end;
  %end;
  drop x1-x&numvars ;
  label trt='Treatment' samp='Rand sample no.' subjno='Patient no.';
run;

proc sort data=randsamp;
  by samp trt;
run;

%mend simul_data;

********************************************************************;
***             iterations to obtain correlations                ***;
********************************************************************;

%macro iterat_simul(n_=100,num_=1000,numvars_=4,criterion=0.025,maxiter=100,
            var1_df_=-1,var2_df_=-5,var3_df_=-2,var4_df_=0.02,var5_df_=,
                    aim12=0.20, aim13=0.20, aim14=0.05, aim15=0, 
                    aim23=0.20, aim24=0.05, aim25=0, 
                    aim34=0.05, aim35=0, 
                    aim45=0,
                    out=finalsamp, out_iter=finaliter); /*output datasets*/

%global gnumvars gnum gtype1 gtype2 gtype3 gtype4 gtype5 gprec ; 
  
%let gnumvars=&numvars_;
%let gnum=&num_;
%let t0 = %sysfunc(datetime()); 

data iterations;
  set _null_;
run;

%simul_data(n=&n_,num=&num_,numvars=&numvars_,
var1_df=&var1_df_,var2_df=&var2_df_,var3_df=&var3_df_,
var4_df=&var4_df_,var5_df=&var5_df_);

%let gprec=1;
%let iterat=0;
%do %while (&gprec gt &criterion and &iterat lt &maxiter); 
                   /*prec (precision) = max difference between actual
                                   correlation and the desired value*/

%let iterat=%eval(&iterat+1);

ods output pearsonCorr=corrns1 
           (keep=variable var:
            rename=(var1-var&numvars_ = varc1-varc&numvars_));
proc corr data=randsamp pearson; /*pearson used even for binary outcomes, 
                                   produces the same corr 
                                   as the biserial point correlation*/
  var var1-var&numvars_; /*careful to exclude var[]2 (survival vars)*/
run;

proc sort data=corrns1;
  by variable;
run;

data randzc;
  set randz;
  rename x1-x&numvars_ = var1-var&numvars_;
run;

ods output pearsonCorr=corrns2 
           (keep=variable var: 
            rename=(var1-var&numvars_ = varcz1-varcz&numvars_));
proc corr data=randzc pearson;
  var var:; /*across random samples*/
run;

proc sort data=corrns2;
  by variable;
run;

data corrns;
  merge corrns1 corrns2;
  by variable;
run;

data iterations;
  format precision 7.5 ;
  retain prectmp 0;
  set iterations corrns (in=a);
  conv=0.005; /*conv=convergence limit*/
  if a then do;
    iteration=&iterat;
    %do p=1 %to 5; 
      %do b=1 %to 5;
        if variable="var&p" then do;
          %if &b gt &p and &b le &numvars_ %then %do; /*corrns above diagonal*/
            diff&p&b=abs(varc&b-&&aim&p&b);
            prectmp=max(prectmp,diff&p&b);
            if diff&p&b lt conv then c&p&b=varcz&b; /*desired corrn attained*/
            else do; /*desired value not attained, thus adjust*/
              if &&aim&p&b gt 0 then c&p&b=min((&&aim&p&b/varc&b)*varcz&b,1); 
                                                             /*+ve corrns*/
              else c&p&b=max((&&aim&p&b/varc&b)*varcz&b,-1); /*-ve corrns*/
            end;
            call symput("c&p&b",trim(left(put(c&p&b,best.)))); 
                                                     /*corrn assigned its new value*/
          %end;
        end;
          %if &b gt &p and (&p gt &numvars_ or &b gt &numvars_) %then %do;
            call symput("c&p&b",trim(left(put(0,best.)))); /*superfluous correlations*/
          %end;
      %end;
    %end;
    if variable="var&numvars_" then call symput("gprec",trim(left(put(prectmp,7.5))));
    precision=prectmp; /*this precision will be compared against the criterion 
                         and the loop terminated if satisfied*/
  end;
run;

%simul_data(n=&n_,num=&num_,numvars=&numvars_,
var1_df=&var1_df_,var2_df=&var2_df_,var3_df=&var3_df_,
var4_df=&var4_df_,var5_df=&var5_df_,
             corr12=&c12, corr13=&c13, corr14=&c14, corr15=&c15, 
                          corr23=&c23, corr24=&c24, corr25=&c25,
                                       corr34=&c34, corr35=&c35,
                                                    corr45=&c45); 
%end;

%let speed=%sysfunc(datetime())-&t0;

data &out;
  set randsamp;
run;

data &out_iter;
  retain iteration variable varc1-varc&numvars_ precision speed ; /*put vars in order*/
  format varc1-varc&numvars_ 5.2;
  set iterations;
  speed=put(&speed,6.2)||' secs';
  keep iteration variable varc1-varc&numvars_ precision speed;
run;

proc sort data=&out_iter;
  by iteration variable;
run;

%mend iterat_simul;


*** end *****************************************;
