//double DDEVF(void *Paramstuff,double *RandNumsPass,size_t dim,int pop,int maxy_t,double (*sim_results)[7])
double DDEVF(void *Paramstuff,gsl_rng *RandNumsPass,size_t dim,int pop,int maxy_t, double hatch, int q)
//double DDEVF2(struct STRUCTURE *Params,int pop,double sim_results[1+Params->MAXT[pop]/7][5])
{
// DDEVF sets up model and calls 0DE_SOLVER, returns Params->sim_results
//printf("DDEVF: population:%d\n",pop);
//printf("pop2=%d\n",pop2);
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
		//declarations of stuff for plotting
double PLOT=1.0;

double Fprob;
double p1,p2,p3;	        // simulation results (p_i is the probability of being in each of the three classes)
double Cpar, Cprob, Cprob2, Cprob3;					//CK//
double Opar, Oprob, Oprob2, Oprob3;					//CK//
double r1, r2, r3, c1, c2;	            //CK// simulation results (r1 is resting spore density, c1 is condidia density)
r1 = 0.0; r2=0.0; c1= 0.0; c2=0.0;

//double cover = Params->PARS[20];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double cover_C = Params->PARS[17];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double cover_R = Params->PARS[20];
double open_C = Params->PARS[24];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals
double open_R = Params->PARS[25];		//CK// effect of cage.  Testing to see if EXP bugs had higher or lower infection than ferals

int FlagDay=0;

//double junk3 = -300;
//printf("Resting spores day: %f\n", junk3);


//printf("maxt=%d dim=%d\n",Params->MAXT[pop],dim);getc(stdin);
//printf("determ: nuV=%f\t nuF=%e\n",Params->PARS[2],Params->PARS[3]);//getc(stdin);

//int MAXT3=(Params->EXPDATA[pop2][Params->MAXT2[pop2]][2]+1)*7;
int MAXT3=maxy_t;
//if(pop==6){MAXT3=77;}  ///can't figure out why plot 6 (UMBS 2012) is fucked up.  says it only has 2 weeks of data, but that's not true at all...

//printf("data set %d has %d days\n",pop,MAXT3);getc(stdin);

//printf("Final exp week: %d\n", MAXT3);		getc(stdin);

//int DIM = Params->PARS[9]+3;
int m = Params->PARS[9];
int n = Params->PARS[8];
int DIM = m+n+4+2+1;

int n1=(VFPass/exposetime)*Params->PARS[8];    //The number of the first group of exposed classes to virus
int n2=n-n1;                  //The number of the first group of exposed classes to virus

double t=h;		double t_next=h;	double t_0=h;	int i;				// time loop and index
double epsilon = pow(10,-6);
double y_ode[DIM];
double rand_nuR[MAXT3];
double rand_nuF[MAXT3];

double ave_R = Params->PARS[50+pop];
//double ave_R = Params->PARS[13];  //average R(0) for all the sites.  Fit like all other general params
double specific_muF = Params->PARS[6];   //general intercept for MAX TEMP decay function for conidia
double specific_nuF = Params->PARS[3];   //Site-specific infection rate for conidia
//double R_end = Params->PARS[16];   //CK//  Flat time for end of resting spore blooming
double rain_P = Params->PARS[21];  //fit param used to scale accumulating rain.
double rain_P2 = Params->PARS[26];  //fit param used to scale accumulating rain.
double rain_P3 = Params->PARS[29];  //fit param used to scale accumulating rain.
double RH_P = Params->PARS[22];   //CK//  Parameter for RH data
double temp_P = Params->PARS[23];   //CK//  Parameter for temperature data
double fourth_size=Params->PARS[28];	//CK// degree day when the bugs reach 4th instar

//double stop1	= Params->PARS[16];  //CK// param used for starting date
//double stop2	= Params->PARS[19];  //CK// param used for starting date
double DDstart	= Params->PARS[27];  //CK// param used for starting date
double DDstop	= Params->PARS[19];  //CK// param used for starting date

Params->size_C = 1.0;
Params->indexR = 0.0;
Params->indexV = 0.0;
double C_end=Params->PARS[16];	  //CK// fit param that turns off new conidia production once a specific size has been reached

double temp_now;  //CIK// used simplify decay functions
double total_rainfall=0.0;  //used to sum up rainfall
int rain_day;  //used to sum up rainfall
int beta;		//used to dictate how many days to go back when accumulating rain
int theta;    //used to determine lag period before calculating accumulated rainfall
//beta=Params->PARS[40+pop];
beta=Params->PARS[7];
//beta=1;
theta=1;
//theta=Params->PARS[13];
//printf("Beta: %d Theta: %d\n", beta, theta);		getc(stdin);


double DDtemp_now;  //CIK// used simplify decay functions


//double alpha = 0.01242568;  //param used to scale accumulating rain.  Will have it be fit later?

//double temp_P = Params->PARS[21];   //CK//  Parameter for temperature data
//double RH_P = Params->PARS[22];   //CK//  Parameter for RH data
//double RH_P = Params->PARS[21];   //CK//  Parameter for RH data



//double specific_nuF = Params->PARS[50+pop];   //Site-specific infection rate for conidia
//double specific_muF = Params->PARS[40+pop];   //Site-specific decay rate for conidia
//Params->muF = specific_muF;
//double size_S = Params->PARS[28];   //Scaling effect of size on susceptibility over time.
//double fourth_size=Params->PARS[16];	//CK// degree day when the bugs reach 4th instar
//double C_end=Params->PARS[19];	  //CK// fit param that turns off new conidia production once a specific size has been reached

double nuF2;
double nuR2;

double FIO_Cc;
double FIO_Cr;
double FIO_Oc;
double FIO_Or;

double DD10=0;    //accumulated degree days about 10 degrees C

//double specific_Rend = Params->PARS[40+pop];   //Site-specific decay rate for conidia


//printf("specific_nuF: %e\n", specific_muF);		getc(stdin);

// ----------------------------------- Generate Random Numbers -------------------------------------------- //
for (i=0;i<=MAXT3;i++)	{
	//rand_nuR[i]=gsl_ran_gaussian(RandNumsPass,Params->PARS[11]); //2nd entry is stdev. (or variance)
	//rand_nuF[i]=gsl_ran_gaussian(RandNumsPass,Params->PARS[12]);
	rand_nuR[i]=0;
	rand_nuF[i]=0;
	//rand_nuF[i]=gsl_cdf_gaussian_Pinv(RandNumsPass[i],Params->PARS[12]);
	//printf("i(%d)\t var1=%f\t var2=%f\t rand parts: nuR=%e\t nuF=%e\n",i,Params->PARS[11],Params->PARS[12],rand_nuR[i],rand_nuF[i]);
}//getc(stdin);
// ------------------------------------- Initial Conditions ---------------------------------------------- //
// THIS SHOULD BE A NEW FUNCTION (should be in main?? outside an extra loop??)//
	//printf("DDEVF check2 HERE!");		getc(stdin);

//CK//  Changing initial conditions for fungus only model.  Don't need to take out virus infected neonates
//CK//  Assume that all hatched bugs are susceptible to fungus
//CK//  Also going to convert to bugs per meter^2 for units of S(0), rather than egg masses per 1/40th of a hectare

/*
double temp_S= (double)Params->DATA[pop][0][0]+(double)Params->DATA[pop][1][0];
double temp_V= (double)Params->DATA[pop][0][1]+(double)Params->DATA[pop][1][1];

if (temp_S<1)	{			//CK// Need to deal with small initial populations and the 2 zero collections in Pop 1
	temp_S= (double)Params->DATA[pop][2][0]+(double)Params->DATA[pop][3][0];
}

if (temp_V<1)	{			// give positive viral ICs if data gives no viral infected
	temp_V = 1.0;
	temp_S = 2.0*temp_S;
}*/

//Params->INITS[0] = Params->PARS[30+pop]*temp_S/(temp_S+temp_V);			// initS
//Params->INITS[1] = Params->PARS[30+pop]*temp_V/(temp_S+temp_V);			// initV
Params->INITS[0] = Params->PARS[30+pop];			// initS
//Params->INITS[1] = Params->INITS[0]*0.01;						// initV  SET TO 0 FOR FUNGUS ONLY MODEL!!
//Params->INITS[3] = Params->PARS[40+pop];								// initR
Params->INITS[3] = ave_R;												// initR  //CK// changed to use average R(0), not site-specific

// END OF NEW FUNCTION //
double initS = Params->INITS[0];			// initS
double initV = VPass;			// initV, passed from VPass in head file
double initR = Params->INITS[3];			// initR

//printf("CHECKING S(0)!! base S(0)=%f\t temp_S=%f\t temp_V = %f\n",Params->PARS[30+pop],temp_S, temp_V);
//printf("Start DDVF pop:%d\t initR =%f\t initialS = %f\t initV = %f\n",pop,initR, initS, initV ); getc(stdin);

//printf("pop:%d\t tempS=%f tempV=%f\n",pop,temp_S,temp_V);
//printf("pop(%d):\t initS=%f\t initR=%f\n",pop,Params->PARS[30+pop],Params->PARS[50+pop]);	getc(stdin);
// printf("th_id=%d\t POP %d\t initS=%f\t initV=%f\t initR=%f\n",Params->th_id,pop,initS,initV,initR);
// printf("th_id=%d POP %d\t initS=%f params_inits=%f tempS=%f tempV=%f pop=%f\n",Params->th_id,pop,initS,Params->INITS[0],temp_S,temp_V,Params->PARS[30+pop]);
// printf("DATA population(%d): initial week: %d\t %d\t %d\n",pop,DATA[pop][0][0],DATA[pop][0][1],DATA[pop][0][2]);	getc(stdin);
// ----------------------------------------- Fixed Parameters ---------------------------------------------- //
int gstepsV		= (int) Params->PARS[8];	int gstepsF	= (int) Params->PARS[9];
//double ratio	= Params->PARS[10];
double ratio = 1;
//double neo_v	= Params->PARS[15];		// latent period of neonates (days)
double neo_v	= 7.0;			// latent period of neonates (days) FUNGUS ONLY MODEL!
//double start1	= Params->PARS[16];  //CK// param used for starting date
//double start2	= Params->PARS[19];  //CK// param used for starting date

//double stop1	= Params->PARS[26];  //CK// param used for starting date
//double stop1	= Params->PARS[26];  //CK// param used for starting date
//double stop2	= Params->PARS[27];  //CK// param used for starting date

//double R_end	= Params->R_END[pop];   //CK//  Change value for function of latitude
double R_end;   //CK//  Change value for function of latitude
double R_start;   //CK//  Change value for function of latitude

//Params->PARS[0] = 1/(Params->PARS[4]);		// THIS SHOULD BE IN PARHOOD.C AFTER PARAMS ARE SELECTED (CAREFUL FOR MCMC.C VERSION)

Params->PARS[0]=1.0;


// ------------------------------------- initialize model parameters --------------------------------------- //
int FlagWeek;	int FlagV=0;	int FlagR=0;	int FlagR_end=0;		// keep track of end of week
FlagV=1;  //At the beginning stage, part of the population is infected.
double ConiBefore=0.0;    //CK// Thing to store conidia the day before the feral collection
double RestBefore=0.0;    //CK// Thing to store resting spores the day before the feral collection

int day = 0; int week = 0;							// keeps track of day and week number

int line_ticker=0;   //CK// Ticker used to associate t in function with numbered days.
int line_ticker2;
int test_day=0;	//CK// used to find the line in the weather data that corresponds to the starting day of collections
int num_day =  hatch;  //CK// Starting day number
//printf("Params->DATA[pop2][0][4]=%d\n",Params->DATA[pop2][0][4]);
//printf("starting day number: %d\n", num_day);		getc(stdin);
//for (i=0;i<MAXT3;i++){
//    printf("Params->WDATA[2][%d][1]=%lf\n",i+1,Params->WDATA[2][i][1]);
//}


line_ticker = num_day;

line_ticker=line_ticker-1;

//printf("corresponding line in WDATA: %d\n", line_ticker);		getc(stdin);

double S,V,F,R;	double IV=0, IF=0, IVF=0;
double E_V[n2+1]; double E_F[gstepsF+1]; double E_VF[n1+1];
int num_weeks=MAXT3/7;

// -----------------------------------//CK// calculating ending blooming times //CK//--------------------------------------- //

DD10=0.0;		R_start = 0.0;
test_day = line_ticker;

while(DD10 <= DDstart){

	DDtemp_now = Params->CCDATA[q][test_day][3]-10.0;  //CK// begin calculation of accumulated Degree Days
	if(DDtemp_now<0.0){DDtemp_now=0.0;}
	DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
	R_start++;
	test_day++;
//printf("temp now: %f DD10: %f R_start: %f test_day: %d\n", DDtemp_now, DD10, R_start, test_day);		//getc(stdin);
}
//printf("temp now: %f DD10: %f R_start: %f test_day: %d\n", DDtemp_now, DD10, R_start, test_day);		//getc(stdin);


//R_start=0.0;

DD10=0.0;		R_end = 0.0;
test_day = line_ticker;

while(DD10 <= DDstop){

	DDtemp_now = Params->CCDATA[q][test_day][3]-10.0;  //CK// begin calculation of accumulated Degree Days
	if(DDtemp_now<0.0){DDtemp_now=0.0;}
	DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time
	R_end++;
	test_day++;
//printf("temp now: %f DD10: %f R_end: %f test_day: %d\n", DDtemp_now, DD10, R_end, test_day);		//getc(stdin);
}
//printf("temp now: %f DD10: %f R_end: %f test_day: %d\n", DDtemp_now, DD10, R_end, test_day);		//getc(stdin);

DD10=0.0;

/*
//Lat and long for all 6 sites in current order (KF, RC1, UM1, RC2, RC3, CCR)
//42.363523, -85.348499
//44.463390, -84.604086
//45.483875, -84.680951
//44.465764, -84.595857
//44.465764, -84.595857
//45.188782, -84.22861
//42.614460, -85.453187  <- YS 2011

double lats[8];

lats[0]=42.363523; lats[1]=44.463390; lats[2]=45.483875; lats[3]=44.465764; lats[4]=44.465764; lats[5]=45.483875; lats[6]=42.614460;

//R_start=start1*exp(start2*(lats[pop-1]-lats[0]));
R_start=0.0;

//R_end=200;  //CK// forgot to deal with the ending day when I repurposed the 16th param for start1.  Deal with later

//R_end=Params->PARS[26];
//R_end=specific_Rend;  //CK// forgot to deal with the ending day when I repurposed the 16th param for start1.  Deal with later

R_end=stop1*exp(stop2*(lats[pop-1]-lats[0]));
//R_end=200.0;

*/
//printf("Checking blooming params: pop=%d\t Rstart= %f\t lat=%f\t r_time= %f\t Rend= %f\n", pop,Params->DAY_F[pop],lats[pop-1], r_time, R_end);		getc(stdin);



//-----------------------------------//CK// Infections on day 0!!!! -------------------------------//
// At week 0, when the bugs hatch, no infections can occur because the bugs are too young.
//sim_results[0][0]=1.0;sim_results[0][1]=0.0;sim_results[0][2]=0.0;  //CK// Feral results for week 0.  No infections on day 0, apparently
//sim_results[0][3]=0.0;sim_results[0][4]=0.0;  //CK// No conidia on day 0
//sim_results[0][5]=initR; sim_results[0][6]=initR;
//printf("R_start: %f\n", R_start);

if(R_start<1.0){

	total_rainfall = 0.0;

	//printf("stochasticity: %f\t initR: %f\n", rand_nuF[0], initR);

	for (rain_day= (line_ticker - beta - theta - 1);rain_day <= line_ticker - theta -1;rain_day++){
		total_rainfall = Params->CCDATA[q][rain_day][0]+total_rainfall;
	//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
	}

	//Params->nuR = (alpha*total_rainfall)*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + rand_nuF[0]);
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall))*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2.0)*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #4

//Deterministic nuR!//
	//Params->nuR = total_rainfall;
	//Params->nuR = exp(rain_P*total_rainfall);
	//Params->nuR = rain_P2*exp(rain_P*total_rainfall + rand_nuR[0]);
	Params->nuR=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[0]);  //JL: Should be the same as DDEVF for MCMC?
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall));	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2.0);	//CK// Resting Spore Rain Response #4
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))*(total_rainfall/(rain_P + total_rainfall)));	//CK// Resting Spore Rain Response #4
//Deterministic nuR!//

	r2=Params->nuR*initR;

	//if(total_rainfall==0.0){Params->nuR=0.0;}

	//Params->nuR = (rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6])*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][2] + RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//printf("initR: %f\n", initR);
	//printf("rain: %f\t stochasticity: %f\t initR: %f\n", total_rainfall, rand_nuF[0], initR);

	//sim_results[0][5]= Params->nuR*initR*(size_S*fourth_size);
	//sim_results[0][5]= Params->nuR*initR;


	total_rainfall = 0.0;

//printf("Pop:%d\t Day:%d\t Week:%d\t Rain1: %f\t Stoch1: %f\t Spores1: %f\n", pop, day, week, total_rainfall, junk3, junk4);
//printf("Pop:%d\t Day:%d\t Week:%d\t Rain1: %f\t Stoch1: %f\t Spores1: %f\n", pop, day, week, total_rainfall, exp(rand_nuF[0]), initR);

	for (rain_day= (line_ticker2 - beta - theta - 1);rain_day <= line_ticker2 - theta -1;rain_day++){
		total_rainfall = Params->CCDATA[q][rain_day][0]+total_rainfall;
	//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
	}

	//Params->nuR = (alpha*total_rainfall)*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + rand_nuF[0]);
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall))*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2)*exp(rand_nuF[0]);	//CK// Resting Spore Rain Response #4

//Deterministic nuR!//
	//Params->nuR = total_rainfall;
	//Params->nuR = exp(rain_P*total_rainfall);
	//Params->nuR = rain_P2*exp(rain_P*total_rainfall + rand_nuR[0]);
	Params->nuR=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[0]);  //JL: Should be the same as DDEVF for MCMC?
	//Params->nuR = (total_rainfall/(rain_P + total_rainfall));	//CK// Resting Spore Rain Response #3
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2.0);	//CK// Resting Spore Rain Response #4
	//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))*(total_rainfall/(rain_P + total_rainfall)));	//CK// Resting Spore Rain Response #4
//Deterministic nuR!//

	//if(total_rainfall==0.0){Params->nuR=0.0;}

	//Params->nuR = (rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6])*exp(rand_nuF[0]);
	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//Params->nuR = exp(rain_P*total_rainfall + temp_P*Params->WDATA[pop][line_ticker - 1][2] + RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[0]);

	//sim_results[0][6]= Params->nuR*initR*(size_S*fourth_size);
	//sim_results[0][6]= Params->nuR*initR;

	c1 = 0.0;	c2 = 0.0;
	r1 = Params->nuR*initR;

	Cpar = cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2);	//covered cages
	Opar = open_C*((c1+c2)/2) + open_R*((r1+r2)/2);	//uncovered cages


	FIO_Cc = cover_C*((c1+c2)/2);
    FIO_Cr =  cover_R*((r1+r2)/2);
    FIO_Oc = open_C*((c1+c2)/2);
    FIO_Or = open_R*((r1+r2)/2);

	Cprob = exp(-Cpar);	Oprob = exp(-Opar);

//	printf("%d\t %d\t %f\t %f\t %f\t %f\t\n",pop,day,prob,1.0-prob, prob2, prob3); // getc(stdin);

	Cprob2 = 0.0;  //fraction of F infection due to conidia   COVERED CAGES!!!
	Cprob3 = 1-Cprob;  //fraction of F infection due to resting spores

	Oprob2 = ((((c1+c2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*((1 - Oprob));  //OPEN fraction of F infection due to conidia
	Oprob3 = ((((r1+r2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*((1 - Oprob));  //fraction of F infection due to resting spores

	if(Oprob == 1.0){Oprob2=0.0; Oprob3=0.0;}

	//printf("%d\t %d\t %f\t %f\t %f\t %f\n",pop, day, c1, c2, r1, r2);

	//printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,0,1.0,0.0,0.0, Cprob,1.0-Cprob,Cprob2,Cprob3, Oprob,1-Oprob,Oprob2,Oprob3); //getc(stdin);

	printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,0,0.0,0.0,0.0, Cprob,1.0-Cprob,Cprob2,Cprob3, Oprob,1-Oprob,Oprob2,Oprob3, FIO_Cc, FIO_Cr, FIO_Oc, FIO_Or,Params->INITS[0], 0.0, 0.0,0.0,0.0,0.0); //getc(stdin);
	r2= r1;
}
else{
	//sim_results[0][5]=0.0; sim_results[0][6]=0.0;
	//printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,0,1.0,0.0,0.0, 1.0,0.0,0.0,0.0, 1.0,0.0,0.0,0.0); //getc(stdin);
    //if(PLOT==1.0){printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,0,0.0,0.0,0.0, 1.0,0.0,0.0,0.0, 1.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,Params->INITS[0],0.0,0.0,0.0,0.0,0.0); //getc(stdin);
	//}

}

//printf("Pop:%d\t Day:%d\t Week:%d\t Rain1: %f\t Stoch1: %f\t Spores1: %f\n", pop, day, week, total_rainfall, exp(rand_nuF[0]), initR); //getc(stdin);

//printf("Pop:%d\t Day:%d\t Week:%d\t Rain: %f\n", pop, day, week, total_rainfall);
//printf("Pop:%d\t Day:%d\t Week:%d\t Conidia1:%f\t Conidia2:%f\t Spores1: %f\t Spores2: %f\n", pop, day, week, sim_results[week][3],sim_results[week][4], sim_results[week][5], sim_results[week][6]); //getc(stdin);

//R=0.0;  //CK//  Make sure R density starts at zero.  Wasn't sure if that was happening....

// ---------------------------- Calculated Parameters  and Population Sizes ------------------------------- //
double Vstart = ratio*initV;						// viral cadavers after infected neonates die

//printf("Starting Virus Conditions: initV=%f\t ratio=%f\t Vstart=%f\n", initV, ratio, Vstart);	//getc(stdin);

//double r_germ = Params->DAY_F[pop]-r_time;		// calculated day of resting spore germination
double r_germ = R_start;		//CK// nixing r_time because I made all germ dates start at beginning of the collections
if (r_germ<0)	r_germ=0;
//if		(initS>Vstart)	{	S = initS-Vstart;	}	// S(0) (number of healthy neonates)
//if		(initS>Vstart)	{	S = initS;	}	// S(0) (number of healthy neonates)
//else					{	S = 0;	printf("Error Rafa\n");			}
S=initS;
V=0; F=0.0; R=0.0;
for (i=1;i<=n2;i++)		{	E_V[i]=0;			}
for (i=1;i<=n1;i++)		{	E_VF[i]=0;		    }
for (i=1;i<=gstepsF;i++)		{	E_F[i]=0;			}
//double timing[6]={neo_v,r_germ,R_end,MAXT3};
double timing[6]={r_germ,R_end,MAXT3};
//printf("r_germ=%e,R_end=%e\n",r_germ,R_end);
// ----------------------------------------- initialize results ------------------------------------------- //
double P[4] = {0,0,0,0};				// model results (0,healthy,viral kill,fungal kill)
double FRnext=0, Vkill=0, Fkill=0, Vnext=0;				// total individuals infected (used in calculating fraction infected)
double Vcadaver=0, Fcadaver=0;
double VirFI[100], FungFI[100];			// fraction of infected individals (by week number)
VirFI[0] = 0.0; VirFI[1] = 0.0; FungFI[0] = 0.0; FungFI[1] = 0.0;		// fraction infected first two weeks
//double Pdelay[64]={0};
//int Vdelayday=12;

// -------------------- MAIN LOOP!! (calculate populations as time is increased) -------------------------- //
//printf("pop:%d\t t_0=%f\t neo_v=%f\t r_germ=%f R_end=%f week_end=%d maxt=%f\n",pop,t_0,timing[0],timing[1],timing[2],7*(week+1),timing[3]);
//printf("begin DDEVF main loop\n");		getc(stdin);
while (t_0<MAXT3+h)	{    //CK// change MAXT to MAXT2 to let it go to the end of the experimental datasets?
	//if (week>=num_weeks)	break;
	FlagWeek=0;



	//printf("t_0=%f\t day=%d\t week=%d\t num_weeks=%d max_time=%d\n",t_0,day,week,num_weeks,Params->MAXT[pop]);//getc(stdin);
	//printf("neo=%f\t rgerm=%f\t rend=%f\t day+1=%d\n",timing[0],timing[1],timing[2],day+1);	getc(stdin);
	// ------------------------------- Find Stoppage Event ----------------------------------------------- //
	// --------------------- end of day -------------------- //
	if (day+1<timing[0] && day+1<timing[1])	{

		FlagDay=1;

		t_next=day+1;
		day++;
		num_day++;
		line_ticker++;

		//printf("end of day: t=%f\n",t_next);

		if (day%7==0)	{
		  //printf("end of week!!: t=%f\n",t_next);
			week++;
			FlagWeek=1;
		}
		//getc(stdin);
	}
/*	// --------------------- infected neonates die -------------------- //
	else if (timing[0]<timing[1] && timing[0]<timing[2] && timing[0]<day+1)			{
		FlagV=1;
		t_next=timing[0];
		timing[0]=999.9;
		//printf("VIRAL INFECTED NEONATES DIE: t=%f\n",t_next);	//getc(stdin);
	}
*/	// --------------------- resting spores bloom -------------------- //
	else if (timing[0]<timing[1] && timing[0]<day+1)	{
		FlagR=1;
		t_next=timing[0];
		timing[0]=999.9;
		//printf("Resting spores bloom: t=%f\n",t_next);		//getc(stdin);
	}
	// --------------------- resting spores done -------------------- //
	else if (timing[1]<timing[0] && timing[1]<day+1)	{
		FlagR_end=1;
		t_next=timing[1];
		timing[0]=999.9;
		timing[1]=999.9;
		//printf("Resting spores done: t=%f\n",t_next);		//getc(stdin);
	}
/*	// --------------------- infected neonates die and end of day -------------------- //
	else if (abs(day+1-timing[0])<epsilon)	{
		FlagV=1;
		t_next=day+1;
		timing[0]=999.9;
		day++;
		num_day++;
		line_ticker++;

		if (day%7==0)	{
			week++;
			FlagWeek=1;
		}
		//printf("end of day and Resting Spores bloom: t=%f\n",t_next);	getc(stdin);
	}
*/	// --------------------- resting spores bloom and end of day -------------------- //
	else if (abs(day+1-timing[0])<epsilon)	{

		FlagDay=1;

		FlagR=1;
		t_next=day+1;
		timing[0]=999.9;
		day++;
		num_day++;
		line_ticker++;

		if (day%7==0)	{
			week++;
			FlagWeek=1;
		}
		//printf("end of day and Resting Spores bloom: t=%f\n",t_next);	//getc(stdin);
	}
	// --------------------- resting spores done and end of day -------------------- //
	else if (abs(day+1-timing[1])<epsilon)	{

		FlagDay=1;

		FlagR_end=1;
		t_next=day+1;
		timing[0]=999.9;
		timing[1]=999.9;
		day++;
		num_day++;
		line_ticker++;

		if (day%7==0)	{
			week++;
			FlagWeek=1;
		}
		//printf("end of day and Resting Spores done: t=%f\n",t_next);	//getc(stdin);
	}
	else {
		printf("ERROR: NO EVENT IS NEXT IN TIME???\n");
		getc(stdin);
	}
	// -------------------------- integrate until next stoppage event ---------------------------------- //
	while (t<t_next)	{
		y_ode[0]=S;	y_ode[m+n+1]=Fcadaver;	Params->POPS[3]=R;
		y_ode[m+n+3]=Vcadaver;
		//printf("pop:%d t=%f\t S=%4.3e\t V=%4.3e\t F=%4.3e\t R=%4.3e\n",pop,t,S,V,F,R);

		for (i=1;i<=gstepsF;i++)	{
			//y_ode[1+i]=E_V[i];
			y_ode[i]=E_F[i];
		}
		for (i=1;i<=n1;i++)	{
			//y_ode[1+i]=E_V[i];
			y_ode[gstepsF+i]=E_VF[i];
		}
		for (i=1;i<=n2;i++)	{
			//y_ode[1+i]=E_V[i];
			y_ode[gstepsF+n1+i]=E_V[i];
		}
		y_ode[m+n+2]=FRnext;
		y_ode[m+n+4]=Vkill;
		y_ode[m+n+5]=Fkill;
		y_ode[m+n+6]=Vnext;

		//for (i=0;i<=m+n+5;i++){
        //    printf("y_ode[%d]=%f\n",i,y_ode[i]);
		//}
		//printf("DDEVF:\n parm 2=%f parm 3=%f parm 4=%f parm 5=%f parm 6=%f\n",Params->PARS[2],Params->PARS[3],Params->PARS[4],Params->PARS[5],Params->PARS[6]);
		//getc(stdin);

		//Params->nuF = specific_nuF;
		//Params->nuR = initR;
		Params->nuV = Params->PARS[2];
		//Params->muF = specific_muF;

		DDtemp_now = Params->CCDATA[q][line_ticker - 1][3]-10.0;  //CK// begin calculation of accumulated Degree Days
		if(DDtemp_now<0.0){DDtemp_now=0.0;}
		DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time

		//if(DD10>=C_end){
        //        Params->size_C=0.0;
                //Params->indexR=1;
		//}	//CK// stops new conidia production once a fit size has been reached.
        //if(DD10>=fourth_size){
        //        Params->size_C=0.0;
        //        Params->indexR=1;
		//}

		//Params->nuF = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6]);
		//Params->nuF = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[(int)t]);
		//Params->nuF = (DD10/fourth_size)*specific_nuF*exp(RH_P*Params->WDATA[pop2][line_ticker - 1][6] * exp(rand_nuF[(int)t]));
		//nuF2 = specific_nuF*exp(RH_P*Params->WDATA[pop2][line_ticker - 1][6] + rand_nuF[(int)t]);
		//nuF2 = specific_nuF*exp(RH_P*Params->WDATA[pop2][line_ticker - 1][6] * exp(rand_nuF[(int)t]));  //JL: Make it the same as DDEVF for MCMC?
		//if(nuF2> pow(8.0,8.0)){nuF2= pow(8.0,8.0);}
		nuF2 = specific_nuF*exp(RH_P*Params->CCDATA[q][line_ticker - 1][1]) * exp(rand_nuF[(int)t]);    //JL: Make it the same as DDEVF for MCMC
		if(nuF2> pow(8.0,8.0)){nuF2= pow(8.0,8.0);}
        //printf("%e\t %e\t %e\n",specific_nuF,RH_P,rand_nuF[(int)t]);
		Params->nuF = (DD10/fourth_size)*nuF2;
        //printf("%e\t %e\n",DD10,fourth_size);


		temp_now = Params->CCDATA[q][line_ticker - 1][2];  //CK// putting max temp in smaller object
		//printf("DDtemp_now=%e, temp_now=%e\n",DDtemp_now,temp_now);

		//Params->muF = specific_muF;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		//Params->muF = specific_muF*temp_now;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		Params->muF = specific_muF*exp(temp_P*temp_now);	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		if(Params->muF > pow(10.0,10.0)){Params->muF= pow(10.0,10.0);}

		//total_rainfall = 0.0;
		//rain_day= line_ticker - beta - 1;
//printf("Start rain calc! day: %d\n", line_ticker-1);		getc(stdin);
		//for (rain_day= (line_ticker - beta - theta - 1);rain_day <= line_ticker - theta -1;rain_day++)	{
		//	total_rainfall = Params->WDATA[pop][rain_day][1]+total_rainfall;
//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);		getc(stdin);
		//}

		total_rainfall = 0.0;
		rain_day= line_ticker - beta - 1;

		//if(Params->POPS[3] > 0.0){
		for (rain_day= (line_ticker - beta - theta - 1);rain_day <= line_ticker - theta -1;rain_day++)	{

			if (rain_day<0){
                total_rainfall=0;
            }
            else{
                    total_rainfall = Params->CCDATA[q][rain_day][0]+total_rainfall;
            }
			//printf("line: %d\t accumulating rain: %lf\n", rain_day, total_rainfall);
		}
		//else{total_rainfall = 0.0;}

		//Params->nuR = total_rainfall;
		//Params->nuR = 1.0-exp(rain_P*total_rainfall);
		//Params->nuR = exp(rain_P*total_rainfall);
		//Params->nuR = rain_P2*exp(rain_P*total_rainfall + rand_nuR[(int)t]);
		//Params->nuR = (DD10/fourth_size)*rain_P2*exp(rain_P*total_rainfall + rand_nuR[(int)t]);
		//nuR2 = rain_P2*exp(rain_P*total_rainfall + rand_nuR[(int)t]);
		//printf("rainfall=%lf\n",total_rainfall);
		nuR2=(rain_P/(1+rain_P2*exp(-rain_P3*total_rainfall)) - rain_P/(rain_P2+1))*exp(rand_nuR[(int)t]);  //JL: Should be the same as DDEVF for MCMC?
		if(nuR2> pow(10.0,10.0)){nuR2= pow(10.0,10.0);}

		Params->nuR = (DD10/fourth_size)*nuR2;
        //printf("%e\t %e\t %e\t %e\n",rain_P,rain_P2,rain_P3,rand_nuR[(int)t]);
		//printf("%e\t %e\n",nuR2,DD10/fourth_size);
		//if(total_rainfall==0.0){Params->nuR=0.0;}
		if(Params->POPS[3] == 0.0){Params->nuR = 0.0; nuR2 = 0.0;}

//printf("day: %d\t line: %d\t accumulated rain: %lf\n", num_day, line_ticker-1, total_rainfall);		getc(stdin);

/*
//Deterministic Params START!//

		//Params->nuV = Params->PARS[2]*exp(rand_nuV[(int)t]);
		Params->nuV = Params->PARS[2];

		DDtemp_now = Params->WDATA[pop][line_ticker - 1][4]-10.0;  //CK// begin calculation of accumulated Degree Days
		if(DDtemp_now<0.0){DDtemp_now=0.0;}
		DD10 = DD10 + DDtemp_now;			//CK// summing degree days over time

		Params->nuR = total_rainfall*(size_S*DD10);
		nuR2 = total_rainfall*(size_S*fourth_size);	//CK//
		//Params->nuR = total_rainfall*(size_S*day);
		//nuR2 = total_rainfall*(size_S*28);
		//Params->nuR = exp(rain_P*total_rainfall);
		//Params->nuR = (total_rainfall/(rain_P + total_rainfall));	//CK// Resting Spore Rain Response #3
		//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))^2.0);	//CK// Resting Spore Rain Response #4
		//Params->nuR = ((total_rainfall/(rain_P + total_rainfall))*(total_rainfall/(rain_P + total_rainfall)));	//CK// Resting Spore Rain Response #4

//		Params->nuF = specific_nuF;
		//Params->nuF = Params->PARS[3];
		//Params->nuF = Params->PARS[3]*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6]);
		Params->nuF = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6])*(size_S*DD10);
		nuF2 = specific_nuF*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6])*(size_S*fourth_size);
		//Params->nuF = RH_P*exp(specific_nuF*Params->WDATA[pop][line_ticker - 1][6]);

		Params->size_C = Params->PARS[29]*DD10;   //Scaling effect of size on susceptibility over time.
		if(DD10>=C_end){Params->size_C=0.0;}	//CK// stops new conidia production once a fit size has been reached.

		temp_now = Params->WDATA[pop][line_ticker - 1][2];  //CK// putting max temp in smaller object

		//Params->muF = specific_muF;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		Params->muF = specific_muF*temp_now;	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		//Params->muF = specific_muF*exp(Params->PARS[26]*temp_now);	//CK// Conidia Decay Response #2.2  BEST SO FAR!!

		//Params->muF = specific_muF*temp_now;   //Conidia decay as a function of MAX DAILY TEMP
		//Params->muF = Params->PARS[26]*exp(specific_muF*temp_now);	//CK// Conidia Decay Response #2.1
		//Params->muF = specific_muF*exp(temp_now);	//CK// Conidia Decay Response #2.1
		//Params->muF = specific_muF*exp(Params->PARS[26]*temp_now);	//CK// Conidia Decay Response #2.2  BEST SO FAR!!
		//Params->muF = Params->PARS[26]*(temp_now/(specific_muF + temp_now));	//CK// Conidia Decay Response #3.1
		//Params->muF = specific_muF*(temp_now/(Params->PARS[26] + temp_now));	//CK// Conidia Decay Response #3.2
		//Params->muF = Params->PARS[26]*((temp_now/(specific_muF + temp_now))*(temp_now/(specific_muF + temp_now)));	//CK// Conidia Decay Response #4.1
		//Params->muF = specific_muF*((temp_now/(Params->PARS[26] + temp_now))*(temp_now/(Params->PARS[26] + temp_now)));	//CK// Conidia Decay Response #4.2

//Deterministic Params END!//
*/
		//Params->nuF = Params->PARS[3]*exp(rand_nuF[(int)t])*total_rainfall*rain_P;
//		//Params->nuF = Params->PARS[3]*exp(RH_P*Params->WDATA[pop][line_ticker - 1][6] + rand_nuF[(int)t]);
		//Params->nuF = Params->PARS[3]*exp(rand_nuF[(int)t]);

		//printf("call fast_ode: \t t=%f\t day=%d\t t_index=%d\t nuV_determ=%f\t nuV=%f\t rand_nuV=%f\n",
		//	t,day,(int)t,Params->PARS[2],Params->nuV,rand_nuV[(int)t]);	//getc(stdin);
		//printf("call fast_ode: \t t=%f\t day=%d\t t_index=%d\t nuF_determ=%f\t nuF=%f\t rand_nuF=%f\n",
		//	t,day,(int)t,Params->PARS[3],Params->nuF,rand_nuF[(int)t]);
		//printf("nuF: %f\n", Params->nuF);
		//printf("resting spores: %e\t nuR: %f\n", Params->POPS[3], Params->nuR );
		//if(Params->nuR > 0.0000000000001){getc(stdin);}

		//CK// Add in calculation of accumulated rain here.
		//Need to link first day of collections with numbered day for weather.  Will allow back calculation here
		//include total rain in Params, which is passed to ODE_Solver.  Make part of PARS.  Need to find end of PARS...


		t=ODE_Solver(t,t_next,Params,y_ode);

		S=y_ode[0];	Fcadaver=y_ode[m+n+1]; Vcadaver=y_ode[m+n+3];
		FRnext=y_ode[m+n+2];			// killed populations
		Vkill=y_ode[m+n+4];
		Fkill=y_ode[m+n+5];			// killed populations
		Vnext=y_ode[m+n+6];

		IF=0;
		IV=0;
		IVF=0;

//printf("pop: %d day: %d conidia: %f\n", pop, day, y_ode[7]*Params->nuF);


//CK// trying to track conidia and resting spore blooming

//double junk3 = -300;
//printf("Resting spores day: %f\n", junk3);


//printf("Pop:%d\t Day:%d\t Host:%f\t Virus:%f\t Conidia:%f\t Resting Spores: %f\t nuR: %f\n", pop, day, y_ode[0], y_ode[1], y_ode[7], Params->POPS[3], Params->nuR);
//if(pop == 1 && y_ode[7] >0.0){getc(stdin);}
//if(Params->POPS[3] >0.0){getc(stdin);}
//getc(stdin);

	if ((day+1)%7==0)	{
		ConiBefore=y_ode[m+1]*nuF2;  //CK// Saving the conidia 24 hours before feral collection
		RestBefore = Params->POPS[3]*nuR2;
		//ConiBefore=y_ode[m+1]*Params->nuF;  //CK// Saving the conidia 24 hours before feral collection
		//RestBefore = Params->POPS[3]*Params->nuR;
		//FlagConidia=2;
	}
		//getc(stdin);

		for (i=1;i<=gstepsF;i++)	{
			//E_V[i]=y_ode[1+i];
			E_F[i]=y_ode[i];
			//IV += E_V[i];
			IF += E_F[i];
			//printf("E_V(%d)=%e\t",i,E_V[i]);	//printf("IV=%f\n",IV);
		}//printf("\n");
		for (i=1;i<=n1;i++)	{
			//E_V[i]=y_ode[1+i];
			E_VF[i]=y_ode[gstepsF+i];
			//IV += E_V[i];
			IVF += E_VF[i];
			//printf("E_V(%d)=%e\t",i,E_V[i]);	//printf("IV=%f\n",IV);
		}//printf("\n");
		for (i=1;i<=n2;i++)	{
			//E_V[i]=y_ode[1+i];
			E_V[i]=y_ode[gstepsF+n1+i];
			//IV += E_V[i];
			IV += E_V[i];
			//printf("E_V(%d)=%e\t",i,E_V[i]);	//printf("IV=%f\n",IV);
		}//printf("\n");
		//if(IF>initS)		IF=initS;
		if (day==MAXT3-7){
            InfFungusWeekbefore=Fkill+IF;
            InfVirusWeekbefore=Vkill+IV;
            //Params->indexV=1;
		}

		if (day==MAXT3-19){
            InfFungusTwoWeekbefore=Fkill+IF;
            InfVirusTwoWeekbefore=Vkill+IV;
            Params->indexR=1;
            Params->size_C=0;
		}

		if (day==MAXT3-30){
            InfFungusMonthbefore=Fkill+IF;
            InfVirusMonthbefore=Vkill+IV;
            //Params->indexV=1;
            //Params->indexR=1;
            //if (Params->size_C==0){
            //    Params->indexR=1;
            //}
            //else {
            //    Params->size_C=0.5;
            //    Params->indexR=0.5;
            //}
		}

	}
	if (FlagV==1)			{	//printf("VIRAL INFECTED NEONATES DIE: t=%f\n",t_next);		//getc(stdin);
        Vcadaver=Vstart;
		FlagV=2;
	}

	//printf("t_0=%lf\n",t_0);
    //int tlag=t_0+Vdelayday;
	//if (tlag<64){
    //    Pdelay[tlag]=V;
    //    printf("Pdelay[%d]=%lf\n",tlag,Pdelay[tlag]);
    //}
	if (FlagR==1)			{	//printf("Resting spores bloom: t=%f\n",t_next);			//	getc(stdin);
		R=initR;
		FlagR=2;
	}

	else if (FlagR_end==1)	{	//printf("Resting spores done: t=%f\n",t_next);				// getc(stdin);
		R=0;
		FlagR_end=2;
	}


//Add daily plotting output here
	//if (FlagWeek==1 && PLOT==1.0)	{
    if (FlagDay==1 && PLOT==1.0)	{

		if ((S+IF)==0)	{
			Fprob=0;
		}
		else					{
			//VirFI[week]  = IV/(S+IV+IF);		// fraction infected at the end of each week
			Fprob = IF/(S+IF);
		}
		p1 = y_ode[m+n+1]*Params->nuF;
		p2 = Params->POPS[3]*Params->nuR;

		//p1 = 1-(Fprob);
		//p2 = 0;  			//CK//	formerly virus infected
		p3 = Fprob;

//		printf("%d\t %d\t %f\t %f\t %f\n",pop,i,p1,p2,p3); //getc(stdin);

		c1 = y_ode[m+n+1]*nuF2;
		r1 = Params->POPS[3]*nuR2;

		//c2=ConiBefore;
		//r2=RestBefore;
		//printf("c1=%e,c2=%e,r1=%e,r2=%e\n",c1,c2,r1,r2);

		Cpar = cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2);
		Opar = open_C*((c1+c2)/2) + open_R*((r1+r2)/2);	//uncovered cages

		//printf("cover_C=%e, cover_R=%e, open_C=%e, open_R=%e\n",cover_C, cover_R, open_C, open_R);
		//printf("c1=%e, c2=%e, r1=%e, r2=%e, nuF2=%e, nuR2=%e\n",c1,c2,r1,r2,nuF2,nuR2);

		FIO_Cc = cover_C*((c1+c2)/2);
		FIO_Cr =  cover_R*((r1+r2)/2);
		FIO_Oc = open_C*((c1+c2)/2);
		FIO_Or = open_R*((r1+r2)/2);

		Cprob = exp(-Cpar);	Oprob = exp(-Opar);

//		printf("%d\t %d\t %f\t %f\t %f\t %f\t\n",pop,day,prob,1.0-prob, prob2, prob3); // getc(stdin);

		Cprob2 = (( cover_C*((c1+c2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*((1 - Cprob));  //OPEN fraction of F infection due to conidia
		Cprob3 = ((cover_R*((r1+r2)/2))/(cover_C*((c1+c2)/2) + cover_R*((r1+r2)/2)))*((1 - Cprob));  //fraction of F infection due to resting spores

		//Cprob2 = 0.0;  //fraction of F infection due to conidia   COVERED CAGES!!!
		//Cprob3 = 1-Cprob;  //fraction of F infection due to resting spores

		Oprob2 = ((open_C*((c1+c2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*((1 - Oprob));  //OPEN fraction of F infection due to conidia
		Oprob3 = ((open_R*((r1+r2)/2))/(open_C*((c1+c2)/2) + open_R*((r1+r2)/2)))*((1 - Oprob));  //fraction of F infection due to resting spores

		//Oprob2 = ((((c1+c2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*((1 - Oprob));  //OPEN fraction of F infection due to conidia
		//Oprob3 = ((((r1+r2)/2))/(((c1+c2)/2) + ((r1+r2)/2)))*((1 - Oprob));  //fraction of F infection due to resting spores

		if(Cprob == 1.0){Cprob2=0.0; Cprob3=0.0;}
		if(Oprob == 1.0){Oprob2=0.0; Oprob3=0.0;}

		//printf("%d\t %d\t %f\t %f\t %f\t %f\n",pop, day, c1, c2, r1, r2);

		//printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,day,p1,p2,p3, Cprob,1.0-Cprob,Cprob2,Cprob3, Oprob,1-Oprob,Oprob2,Oprob3); //getc(stdin);
		//printf("IF=%e,S=%e\n",IF,y_ode[0]);
        //printf("nuF=%e,nuR=%e\n",Params->nuF,Params->nuR);
        //printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,day,p1,p2,p3, Cprob,1.0-Cprob,Cprob2,Cprob3, Oprob,1-Oprob,Oprob2,Oprob3, FIO_Cc, FIO_Cr, FIO_Oc, FIO_Or, y_ode[0], IF, y_ode[m+n+1],V,Fkill,Vkill); //getc(stdin);

        //printf("muF=%e\n",Params->muF);
        //printf("Vcadaver=%e\n", Vcadaver);
        //printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,day-1,initS,y_ode[0],Fkill,Vkill,Fcadaver,Vcadaver,IF,IV,y_ode[0]+Fkill+Vkill+IV+IF); //getc(stdin);

		c2=c1;	r2=r1;   //make today's C and R yesterday's C and R
        //printf("c2=%e,r2=%e\n",c2,r2);
		FlagDay=0;
		//FlagWeek=0;
	//printf("Rstart:%e\t Rend:%e\t R=%e\t nuF=%e\t nuR=%e\t muF=%e\n",R_start, R_end,R,Params->nuF,Params->nuR,Params->muF);
	//printf("%d\t %e\t %e\n",day-1,Params->size_C,Params->indexR);
	//printf("%d\t %d\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",pop,day-1,initS,y_ode[0],Fkill,FRnext,Vkill,Vnext,Fcadaver,Vcadaver,Fkill+IF+IVF,Vkill+IV,y_ode[0]+Fkill+Vkill+IV+IF+IVF);

	}

//Done with daily plotting output here


	// ---------------------------- end of the week updates ---------------------------------- //
/*
	if (FlagWeek==1)	{
		//printf("t=%f S=%e IV=%e IF=%e\n",t,S,IV,IF);
		//week++;							// add one to the week if 7 days are complete
		if ((S+IF)==0)	{
			VirFI[week]=0;	FungFI[week]=0;
		}
		else					{
			//VirFI[week]  = IV/(S+IV+IF);		// fraction infected at the end of each week
			FungFI[week] = IF/(S+IF);
		}
		P[1] = 1-(FungFI[week]);
		P[2] = 0;  			//CK//	formerly virus infected
		P[3] = FungFI[week];

		//printf("week=%d:\t VirFI=%f\t FungFI=%f\n",week,VirFI[week],FungFI[week]);
		sim_results[week][0]=P[1];	sim_results[week][1]=P[2];	sim_results[week][2]=P[3];
		//printf("pop=%d week %d results: P1=%f P2=%f P3=%f\n",pop,week,sim_results[week][0],sim_results[week][1],sim_results[week][2]);	//getc(stdin);

		//STORING THE F AND R AT EACH WEEK
		sim_results[week][3]=ConiBefore;  //CK// Saving the conidia 24 hours before feral collection

		sim_results[week][4]=y_ode[m+1]*Params->nuF;
		//sim_results[week][4]=y_ode[m+1]*nuF2;

		sim_results[week][5]=RestBefore;
		sim_results[week][6]=Params->POPS[3]*Params->nuR;
		//sim_results[week][6]=Params->POPS[3]*nuR2;


//printf("Pop:%d\t Day:%d\t Week:%d\t Conidia1:%f\t Conidia2:%f\t Spores1: %f\t Spores2: %f\n", pop, day, week, sim_results[week][3],sim_results[week][4], sim_results[week][5], sim_results[week][6]); getc(stdin);
//if( y_ode[7] >0.0){getc(stdin);}
//if(Params->POPS[3] >0.0){getc(stdin);}
//getc(stdin);

	}
*/
//printf("t_0=%lf\n",t_0);
t_0=t_next;
//printf("t_0=%lf\n",t_0);
}
SusEnd=y_ode[0];

if (initR==0){
   InfVirusEnd=Vkill+IV+IVF;
   InfFungusNext=0;
   InfFungusEnd=0;
   InfVirusNext=Vkill+IV+IVF;
}
else{
   InfVirusEnd=Vkill+IV;
   InfFungusNext=FRnext+IF+IVF;
   InfFungusEnd=Fkill+IF+IVF;
   InfVirusNext=Vkill+IV;
}
InfFungusAdj=InfFungusEnd-InfFungusWeekbefore;
InfVirusAdj=InfVirusEnd-InfVirusWeekbefore;

//InfFungusAdj=InfFungusEnd-InfFungusTwoWeekbefore;
//InfVirusAdj=InfVirusEnd-InfVirusTwoWeekbefore;

//InfFungusAdj=InfFungusEnd-InfFungusMonthbefore;
//InfVirusAdj=InfVirusEnd-InfVirusMonthbefore;

//printf("%e\t %e\n",InfFungusWeekbefore,InfFungusEnd);

//printf("After an epizootic: %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",SusEnd,Fkill,Vkill,Fcadaver,Vcadaver,IF,IV,InfFungusEnd+InfVirusEnd+SusEnd); //getc(stdin);


//int xx;
//for (xx=0;xx<64;xx++){
//    printf("Pdelay[%d]=%lf\n",xx,Pdelay[xx]);
//}
//printf("Resting spores day -1: %f\n", sim_results[0][5]);

//printf("end DDEVF\n");	getc(stdin);
return 0;
}
