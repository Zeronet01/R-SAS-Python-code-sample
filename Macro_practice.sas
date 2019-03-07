
/*********************************************************************************************************************************
Macro name:		statistic

Purpose:		Calculate descriptive statistics by the variables of interest. 
                Variables in different types can be inputed simultaneously by inputing 0 and 1 in typelist according to the order of varlist.
Author:			Lingyi Tan

Creation Date: 	Oct.11

Revision Date:	Oct.13

SAS version:	9.3

Required 
Parameters:		data = 		data set name
				varlist = 	name of the interested variables
				typelist
					0 = continuous variable
					1 = categorical variable
				
				
Optional 
Parameters:     stat = 		specific statistic that would like to calculate for categorical data


Example: Calculate mean of contiguous variable ht2 and frequency of categorical variable wt9

		%statistic(data=bgd, varlist=ht2 wt9, typelist=0 1, stat=mean);
*********************************************************************************************************************************/


%let pathname=/folders/myfolders/sasuser.v94/note4／;

options symbolgen;
%macro statistic(data=,varlist=,typelist=,stat=);
  %let i=1;
    %let j=1;
  %do %until(%scan(&varlist, &i)=);
  %do %until(%scan(&varlist, &j)=);
     %let variable=%scan(&varlist,&i);
     %let type=%scan(&typelist,&j);
       %if &type=0 %then %do;
	    proc means data=&data &stat;
          var &variable;
          title "The &stat of &variable";
        run;
        %end;
        %else %if &type=1 %then %do;
          proc freq data=&data;
          table &variable/LIST;
          title "&variable Frequencies data";
          run;
          %end;
    %let i=%eval(&i+1);
    %let j=%eval(&j+1);
 %end;
 %end;
   
%mend statistic;


/*Demonstrate the macro statistic using the Berkeley Guidance data set*/

%let pathname=/folders/myfolders/sasuser.v94/note4／;

proc import datafile=" /folders/myfolders/sasuser.v94/note4/Berkeley Guidance Data.txt" out=bgd dbms=dlm replace;
 delimiter=',';
   getnames=yes;
label wt2="Weight at age 2 (kg)";
label ht2="Height at age 2 (cm)";
label wt9="Weight at age 9";
label ht9="Height at age 9";
label lg9="Leg circumference at age 9 (cm)";
label st9="A composite measure of strength at age 9 (high values=stronger)";
label wt18="Weight at age 18";
label ht18="Height at age 18";
label lg18="Leg circumference at age 18";
label SOMA="Somatotype, seven point scale, 1=slender, 7=fat";
run;

/*Suppose that we have ht2 as continuous variable and wt9 as categorical variable*/
%statistic(data=bgd, varlist=ht2 wt9, typelist=0 1, stat=max);
