
void GetString(int j,int i,char* strOut, unsigned int strSize)
{
// ---------------------------------------- Name for Output Files ------------------------------------------- //
char *Prefix=" ";
char *max_lhood_prefix="full_R+T+RH_GRID6";
char *mcmc_prefix="big8june_";
char *profile_prefix="intrinsic_";

char *FileName	= " ";
char *Path		= " ";
char *Name		= " ";
char *Type		= ".dat";						// type of output file
char profile[33];
// ---------------------------------------- create time suffixes -------------------------------------------- //
struct tm *ptr;									// struct defined in time.h;  used to hold the time and date
time_t lt;										// type defined in time.h; for storing the calendar time
lt = time(NULL);
ptr = localtime(&lt);							// converts calendar time to local time
char Date[30];									// string to hold date and time
strftime(Date, 30, "%m_%d_%H_%M_%S", ptr);		// adds to Date month, day, year, hour, minute, second from ptr
// ---------------------------------- create paths and allocate memory ------------------------------------- //
if	(j==0)	{										// MCMC
	Path="/home/jiaweiliu/likelihood_calculation_full/mcmc_results/";
	Prefix=mcmc_prefix;
	if		(i==1)	{	Name = "parms_";											}
	else if	(i==2)	{	Name = "pc_";												}
	else if (i==3)	{	Name = "acc_";												}
	else if (i==4)	{	Name = "lastpars_";											}
	else			{	printf("bad i value in filenames!!!\n");	getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name)+strlen(Type)+strlen(Date)+1), sizeof(char));
}
else if (j==1)	{									// Line Search MLE
	Prefix=max_lhood_prefix;
	if		(i==0)	{	Path="/home/jiaweiliu/likelihood_calculation_full/max_lhood/";			Name="max_lhood";			}
	else if (i==4)	{	Path="/home/jiaweiliu/likelihood_calculation_full/max_lhood/";			Name="L_";				}
	else			{	printf("bad i value in filenames!!!\n");		getc(stdin);	}
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+1), sizeof(char));
}
else if (j>1)	{									// Profile Likelihoods
	Path="/home/jiaweiliu/likelihood_calculation_full/profile_data/";
	Prefix=profile_prefix;
	sprintf(profile,"%d",j);
	Name="pro_";
	FileName = (char*)calloc((strlen(Path)+strlen(Prefix)+strlen(Name) + strlen(Type)+strlen(profile) +1), sizeof(char));
}
// ---------------------------------- add strings onto the file name ------------------------------------------ //
strcat(FileName,Path);
strcat(FileName,Prefix);
strcat(FileName,Name);

if		(j==0||(j==1 && i==4))	strcat(FileName,Date);		// add Date to output file name		(MCMC & max lhood)
else if (j>1)					strcat(FileName,profile);	// add profile number to file name	(profile likelihood)

strcat(FileName,Type);

if (j==0)			strncpy(strOut,FileName,strSize);
else if (j>=1)		strncpy(strOut,FileName,strSize);

return;
}



void output_file(void *Paramstuff,FILE *fp_results,double best_post_hood,int num_adj_pars,int pro)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int i,pop;

//fprintf(fp_results,"F only AVE T 5000 loops I: %d\t %f\t %f\t",pro,Params->PARS[pro],best_post_hood);
fprintf(fp_results,"%d\t %f\t",pro,best_post_hood);

//for (i=0;i<=50+DATA_SETS;i++)	{
for (i=0;i<=num_adj_pars;i++)	{
	fprintf(fp_results,"%e\t",Params->PARS[i]);
	//printf("%e\n",Params->PARS[i]);
}
//for (pop=1;pop<=DATA_SETS;pop++)	{
//	fprintf(fp_results,"%e,\t",Params->MLE_host[pop]);
//	//printf("%e\n",Params->PARS[i]);
//}
for (pop=1;pop<=DATA_SETS;pop++)		{
	fprintf(fp_results,"%e\t",Params->MLE_initR[pop]);
	//printf("%e\n",Params->PARS[i]);
}

fprintf(fp_results,"\n");


}
/*
void output_fileINFECTION(void *Paramstuff,FILE *fp_results,double INFECTED,double Rmean, double Rvar,int num_adj_pars,int pro)
{
STRUCTURE* Params;
Params = (STRUCTURE*) Paramstuff;
int i,pop;

fprintf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", Params.PARS[30+2], Params.PARS[50+1], INFECTED, Params.survivors, Params.total, Rmean, Rvar); //getc(stdin);

fprintf(fp_results,"\n");


}
*/




/*
void printf_last_parms(char fp_last_parms)
{
//char strFileName_last[99];
//GetString(0,4,strFileName_last,98);		printf("%s\n",strFileName_last);
FILE *fp_last_parms;

fp_last_parms = fopen(strFileName_last,"w");

for (i=1;i<=50+DATA_SETS;i++)	{
	fprintf(fp_last_parms,"%e\n",PARAMS[0].PARS[i]);
}
fclose(fp_last_parms);

return;
}
*/
