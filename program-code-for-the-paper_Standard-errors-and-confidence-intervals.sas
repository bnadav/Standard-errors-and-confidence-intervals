%macro entire_program (	/* Macro for the entire program */
	nr				= ,	/* Number of students accepted */
	selr			= ,	/* Selection ratio = proportion accepted */
	rhoxy			= ,	/* Correlation between X and Y */
	rhoxz			= ,	/* Correlation between X and Z */
	rhoyz			= ,	/* Correlation between Y and Z */
	seed			= ,	/* Seed value for the random number generator */
	ntrials			= ,	/* Number of trials */
	nbootsamples	= ,	/* Number of bootstrap samples in each trial */
	alpha			= ,	/* Significance level */
	cond			= 	/* Special condition */
);

/* Create trials */

%let param = 0;

data sasuser.data1;
	keep x eps1 eps2 param trial rhoxy rhoxz rhoyz selr nr;
	%do a = 1 %to %sysfunc (countw (&nr., ' '));
		%do b = 1 %to %sysfunc (countw (&selr., ' '));
			%do c = 1 %to %sysfunc (countw (&rhoxy., ' '));
				%do d = 1 %to %sysfunc (countw (&rhoxz., ' '));
					%do e = 1 %to %sysfunc (countw (&rhoyz., ' '));
						%let param = %sysevalf (&param. + 1);
						%let nrr = %scan (&nr., &a., %str ( ));
						%let selrr = %scan (&selr., &b., %str ( ));
						%let rhoxyy = %scan (&rhoxy., &c., %str ( ));
						%let rhoxzz = %scan (&rhoxz., &d., %str ( ));
						%let rhoyzz = %scan (&rhoyz., &e., %str ( ));
						%let ncandidates = %sysevalf (&nrr. / &selrr.);
						selr = &selrr.;
						nr = &nrr.;
						rhoxy = &rhoxyy.;
						rhoxz = &rhoxzz.;
						rhoyz = &rhoyzz.;
						do f = 1 to &ntrials.;
							trial = f;
							do g = 1 to &ncandidates.;
								param = &param.;
								x = rannor (&seed.);
								eps1 = rannor (%sysevalf (&seed. + 111111));
								eps2 = rannor (%sysevalf (&seed. + 123456));
								output;
							end;
						end;
					%end;
				%end;
			%end;
		%end;
	%end;
run;

data sasuser.data1;
	set sasuser.data1;
	keep param trial x y z rhoxy rhoxz rhoyz selr nr;
	y = rhoxy * x + eps1 * sqrt (1 - rhoxy * rhoxy);
	z = ((rhoxz - rhoxy * rhoyz) / (1 - rhoxy * rhoxy)) * x +
		((rhoyz - rhoxy * rhoxz) / (1 - rhoxy * rhoxy)) * y +
		eps2 * sqrt (1 - ((rhoxz - rhoxy * rhoyz) / (1 - rhoxy * rhoxy)) *
		((rhoxz - rhoxy * rhoyz) / (1 - rhoxy * rhoxy)) -
		((rhoyz - rhoxy * rhoxz) / (1 - rhoxy * rhoxy)) *
		((rhoyz - rhoxy * rhoxz) / (1 - rhoxy * rhoxy)) -
		2 * ((rhoxz - rhoxy * rhoyz) / (1 - rhoxy * rhoxy)) *
		((rhoyz - rhoxy * rhoxz) / (1 - rhoxy * rhoxy)) * rhoxy);
	where &cond.;
run;

/* Select accepted */

proc sort data = sasuser.data1;
	by param trial descending z;
run;

data sasuser.data1;
	set sasuser.data1;
	by param trial descending z;
	if first.trial then obs = 1;
	else obs + 1;
run;

data sasuser.data1;
	set sasuser.data1;
	accepted = 1;
	if obs > nr then do; * Not accpeted.;
		x = .;
		y = .;
		accepted = 0;
	end;
	drop obs;
run;

/* Calculate corrected correlation in trials */

proc corr data = sasuser.data1 outp = sasuser.corr1 noprint;
	by param trial;
	var x y z;
run;

data sasuser.corr1;
	set sasuser.corr1;
	if _type_ = 'CORR' and not (_name_ = 'z');
	if _name_ = 'x' then do;
		rxy = y;
		rxz = z;
	end;
	if _name_ = 'y' then ryz = z;
run;

data sasuser.corr1;
	set sasuser.corr1;
	keep param trial rxy rxz ryz;
	lagrxy = lag (rxy);
	lagrxz = lag (rxz);
	if rxy = . then rxy = lagrxy;
	if rxz = . then rxz = lagrxz;
	if _name_ = 'y';
run;

proc means data = sasuser.data1 noprint;
	output out = sasuser.means1 var = largesz2;
	by param trial;
	var z;
run;

data sasuser.means1;
	set sasuser.means1;
	keep param trial largesz2;
run;

proc means data = sasuser.data1 noprint;
	output out = sasuser.means2 var = smallsz2;
	by param trial;
	where accepted = 1;
	var z;
run;

data sasuser.means2;
	set sasuser.means2;
	keep param trial smallsz2 nstud;
	nstud = _freq_;
run;

data sasuser.corr2;
	merge sasuser.corr1 sasuser.means1 sasuser.means2;
	by param trial;
	rc3 = (rxy + rxz * ryz * (largesz2 / smallsz2 - 1)) /
		  (sqrt ((1 + rxz ** 2 * (largesz2 / smallsz2 - 1)) * (1 + ryz ** 2 * (largesz2 / smallsz2 - 1))));
run;

proc means data = sasuser.corr2 noprint;
	output out = sasuser.means3 mean = rc3mean std = sde;
	by param;
	var rc3;
run;

data sasuser.means3;
	set sasuser.means3;
	keep param rc3mean sde;
run;

/* Calculate standard deviation according to formula */

data sasuser.var1;
	merge sasuser.corr1 sasuser.means1 sasuser.means2;
	by param trial;
	varrxy = ((1 - rxy ** 2) ** 2) / (nstud - 1);
	varrxz = ((1 - rxz ** 2) ** 2) / (nstud - 1);
	varryz = ((1 - ryz ** 2) ** 2) / (nstud - 1);
	covrxyrxz = (ryz * (1 - rxy ** 2 - rxz ** 2) - 0.5 * rxy * rxz * (1 - rxy ** 2 - rxz ** 2 - ryz ** 2))
				/ (nstud - 1);
	covrxyryz = (rxz * (1 - rxy ** 2 - ryz ** 2) - 0.5 * rxy * ryz * (1 - rxy ** 2 - rxz ** 2 - ryz ** 2))
				/ (nstud - 1);
	covrxzryz = (rxy * (1 - rxz ** 2 - ryz ** 2) - 0.5 * rxz * ryz * (1 - rxy ** 2 - rxz ** 2 - ryz ** 2))
				/ (nstud - 1);
	drxyrxy = ((1 + (largesz2 / smallsz2 - 1) * rxz ** 2) ** (-.5)) * ((1 + (largesz2 / smallsz2 - 1) *
			  ryz ** 2) ** (-.5));
	drxyrxz = (largesz2 / smallsz2 - 1) * (ryz - rxy * rxz) * ((1 + (largesz2 / smallsz2 - 1) * ryz ** 2)
			  ** (-.5)) * ((1 + (largesz2 / smallsz2 - 1) * rxz ** 2) ** (-1.5));
	drxyryz = (largesz2 / smallsz2 - 1) * (rxz - rxy * ryz) * ((1 + (largesz2 / smallsz2 - 1) * rxz ** 2)
			  ** (-.5)) * ((1 + (largesz2 / smallsz2 - 1) * ryz ** 2) ** (-1.5));
	varfor = (drxyrxy ** 2) * varrxy + (drxyrxz ** 2) * varrxz + (drxyryz ** 2) * varryz +
			 2 * drxyrxy * drxyrxz * covrxyrxz + 2 * drxyrxy * drxyryz * covrxyryz +
			 2 * drxyrxz * drxyryz * covrxzryz; * Variance from formula. ;
	sef = sqrt (varfor);
run;

proc means data = sasuser.var1 noprint;
	output out = sasuser.var2 mean = varfor;
	by param;
	var varfor;
run;

data sasuser.var2;
	set sasuser.var2;
	keep param sefmean;
	sefmean = sqrt (varfor);
run;

/* Create bootstrap samples */

proc sort data = sasuser.data1;
	by param trial accepted;
run;

proc surveyselect data = sasuser.data1 outhits seed = %sysevalf (&seed. + 654321) reps = &nbootsamples.
	method = urs out = sasuser.data2 samprate = 1 noprint;
	strata param trial accepted;
run;

proc sort data = sasuser.data2 (drop = numberhits expectedhits samplingweight);
	by param trial replicate;
run;

/* Calculate corrected correlation in bootstrap samples */

proc corr data = sasuser.data2 outp = sasuser.corr3 noprint;
	by param trial replicate;
	var x y z;
run;

data sasuser.corr3;
	set sasuser.corr3;
	if _type_ = 'CORR' and not (_name_ = 'z');
	if _name_ = 'x' then do;
		rxy = y;
		rxz = z;
	end;
	if _name_ = 'y' then ryz = z;
run;

data sasuser.corr3;
	set sasuser.corr3;
	keep param trial replicate rxy rxz ryz;
	lagrxy = lag (rxy);
	lagrxz = lag (rxz);
	if rxy = . then rxy = lagrxy;
	if rxz = . then rxz = lagrxz;
	if _name_ = 'y';
run;

proc means data = sasuser.data2 noprint;
	output out = sasuser.means4 var = largesz2;
	by param trial replicate;
	var z;
run;

data sasuser.means4;
	set sasuser.means4;
	keep param trial replicate largesz2;
run;

proc means data = sasuser.data2 noprint;
	output out = sasuser.means5 var = smallsz2;
	by param trial replicate;
	where accepted = 1;
	var z;
run;

data sasuser.means5;
	set sasuser.means5;
	keep param trial replicate smallsz2;
run;

data sasuser.corr4;
	merge sasuser.corr3 sasuser.means4 sasuser.means5;
	by param trial replicate;
	rc3star = (rxy + rxz * ryz * (largesz2 / smallsz2 - 1)) /
		  	  (sqrt ((1 + rxz ** 2 * (largesz2 / smallsz2 - 1)) * (1 + ryz ** 2 * (largesz2 / smallsz2 - 1))));
run;

proc means data = sasuser.corr4 noprint;
	output out = sasuser.means6 mean = rc3stardot;
	by param trial;
	var rc3star;
run;

data sasuser.means6;
	set sasuser.means6;
	keep param trial rc3stardot;
run;

data sasuser.corr5;
	merge sasuser.corr4 sasuser.means6;
	by param trial;
	diff2 = (rc3star - rc3stardot) ** 2;
run;

proc means data = sasuser.corr5 noprint;
	output out = sasuser.corr6 sum = summ;
	by param trial;
	var diff2;
run;

data sasuser.corr6;
	set sasuser.corr6;
	keep param trial seb;
	seb = sqrt (summ / (&nbootsamples. - 1));
run;

proc means data = sasuser.corr6 noprint;
	output out = sasuser.corr7 mean = sebmean;
	by param;
	var seb;
run;

data sasuser.corr7;
	set sasuser.corr7;
	keep param sebmean;
run;

/* Calculate BSI */

data sasuser.corr8;
	merge sasuser.corr2 sasuser.corr6;
	by param trial;
	bsi_l = rc3 - quantile ('normal', %sysevalf (1 - &alpha. / 2)) * seb;
	bsi_u = rc3 + quantile ('normal', %sysevalf (1 - &alpha. / 2)) * seb;
run;

/* Calculate BPI */

proc sort data = sasuser.corr4;
	by param trial rc3star;
run;

%let l = floor (%sysevalf (&nbootsamples. * &alpha. / 2));
%let u = ceil (%sysevalf (&nbootsamples. * (1 - &alpha. / 2)));

data sasuser.corr4;
	set sasuser.corr4;
	by param trial rc3star;
	if first.trial then rank = 1;
	else rank + 1;
	if rank = &l. then bpi_l = rc3star;
	if rank = &u. then bpi_u = rc3star;
	if bpi_l ~= . or bpi_u ~= .;
run;

data sasuser.corr4;
	set sasuser.corr4;
	by param trial rc3star;
	keep param trial bpi_l bpi_u;
	lagbpi_l = lag (bpi_l);
	if bpi_l = . then bpi_l = lagbpi_l;
	if last.trial;
run;

/* Calculate FSI */

data sasuser.var5;
	merge sasuser.corr2 sasuser.var1;
	by param trial;
	fsi_l = rc3 - quantile ('normal', %sysevalf (1 - &alpha. / 2)) * sef;
	fsi_u = rc3 + quantile ('normal', %sysevalf (1 - &alpha. / 2)) * sef;
run;

/* Accuracy of case III correction */

data sasuser.param1;
	set sasuser.data1;
	drop x y z accepted;
	by param;
	if first.param;
run;

data sasuser.means7;
	merge sasuser.means3 sasuser.param1;
	by param;
	biasr = 100 * ((rc3mean - rhoxy) / rhoxy);
	absbiasr = abs (biasr);
run;

proc means data = sasuser.means7 noprint;
	output out = sasuser.means8 mean = maper;
	var absbiasr;
run;

/* Accuracy of bootstrap standard error */

data sasuser.means9;
	merge sasuser.means3 sasuser.corr7;
	by param;
	biasseb = 100 * ((sebmean - sde) / sde);
	absbiasseb = abs (biasseb);
run;

proc means data = sasuser.means9 noprint;
	output out = sasuser.means10 mean = mapeseb;
	var absbiasseb;
run;

/* Accuracy of formula standard error */

data sasuser.var3;
	merge sasuser.means3 sasuser.var2;
	by param;
	biassef = 100 * ((sefmean - sde) / sde);
	absbiassef = abs (biassef);
run;

proc means data = sasuser.var3 noprint;
	output out = sasuser.var4 mean = mapesef;
	var absbiassef;
run;

/* Coverage probability of bootstrap confidence intervals */

data sasuser.corr9;
	merge sasuser.corr8 sasuser.param1;
	by param;
	if bsi_l <= rhoxy <= bsi_u then withinbsi = 1;
	else withinbsi = 0;
run;

proc means data = sasuser.corr9 noprint;
	output out = sasuser.corr10 mean = percwithinbsi;
	var withinbsi;
	by param;
run;

data sasuser.corr10;
	set sasuser.corr10;
	percwithinbsi = 100 * percwithinbsi;
run;

data sasuser.corr11;
	merge sasuser.corr4 sasuser.param1;
	by param;
	if bpi_l <= rhoxy <= bpi_u then withinbpi = 1;
	else withinbpi = 0;
run;

proc means data = sasuser.corr11 noprint;
	output out = sasuser.corr12 mean = percwithinbpi;
	var withinbpi;
	by param;
run;

data sasuser.corr12;
	set sasuser.corr12;
	percwithinbpi = 100 * percwithinbpi;
run;

/* Coverage probability of formula confidence intervals */

data sasuser.var6;
	merge sasuser.var5 sasuser.param1;
	by param;
	if fsi_l <= rhoxy <= fsi_u then withinfsi = 1;
	else withinfsi = 0;
run;

proc means data = sasuser.var6 noprint;
	output out = sasuser.var7 mean = percwithinfsi;
	var withinfsi;
	by param;
run;

data sasuser.var7;
	set sasuser.var7;
	percwithinfsi = 100 * percwithinfsi;
run;

/* Output */

data sasuser.output1;
	set sasuser.means7;
run;

proc print data = sasuser.output1 noobs;
	var param nr selr rhoxy rhoxz rhoyz rc3mean biasr;
run;

data sasuser.output2;
	set sasuser.means8;
run;

proc print data = sasuser.output2 noobs;
	var maper;
run;

data sasuser.output3;
	merge sasuser.means3 sasuser.var3 sasuser.means9 sasuser.param1;
	by param;
run;

proc print data = sasuser.output3 noobs;
	var param nr selr rhoxy rhoxz rhoyz sde sefmean sebmean biassef biasseb;
run;

data sasuser.output4;
	merge sasuser.means10 sasuser.var4;
	by _type_;
run;

proc print data = sasuser.output4 noobs;
	var mapesef mapeseb;
run;

data sasuser.output5;
	merge sasuser.corr10 sasuser.corr12 sasuser.var7 sasuser.param1;
	by param;
run;

proc print data = sasuser.output5 noobs;
	var param nr selr rhoxy rhoxz rhoyz percwithinfsi percwithinbsi percwithinbpi;
run;

%mend entire_program;

/*%entire_program (nr = 20 60 100, selr = 0.1 0.3 0.5, rhoxy = 0.2 0.5 0.8, rhoxz = 0.2 0.8,*/
/*				 rhoyz = 0.3 0.6, seed = 987654, ntrials = 1000, nbootsamples = 2000, alpha = 0.05,*/
/*				 cond = not (rhoxz eq 0.8 and rhoyz eq 0.6));*/
%entire_program (nr = 20, selr = 0.3, rhoxy = 0.5, rhoxz = 0.2 0.8, rhoyz = 0.3 0.6, seed = 987654,
				 ntrials = 10, nbootsamples = 50, alpha = 0.05, cond = not (rhoxz eq 0.8 and rhoyz eq 0.3));
/*%entire_program (nr = 100, selr = 0.1, rhoxy = 0.5, rhoxz = 0.2,*/
/*				 rhoyz = 0.3, seed = 987654, ntrials = 1000, nbootsamples = 2000, alpha = 0.05,*/
/*				 cond = 1);*/