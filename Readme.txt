BIRECT-V algorithm ( is a modification of the BIRECT (BIsecting of RECTangles).

How to run this program?

Go to Birect2.m file, and comment lines 1 and 4  if you want to run many or all test problems from (there are 54 test functions described in the file f.m). 

Then go to file all examples and run for example functions from 1 to 12. Click on the file resultats (results). The program will display:  The function number, the obtained value of f : f(xmin), and the number of function evaluations (f_eval), with the average and median values. 

You can also get the CPU time displayed by uncomment corresponding lines in the file (tous_exemples). If you want to run only one function then uncomment the lines 1 and 4 from Birect2.m and type the number of the function in exmple (example). 

You can use graph=1 (for two dimensional functions), for the graphic representation, and 0 if not. In total 6 figures and graphs are displayed.

The domain of certain functions have been modified (ex. fct. 13, 39, etc..) to get the required results (f_min and x_min). 

In the file po_opt.m, you can either run BIRECT-v or BIRECT-vl :
uncomment ind2=find(mn==D(I),1,'first'); and comment ind2=find(mn==D(I)); to run BIRECT-vl %%%
%%% do the opposite if you want to run BIRECT-v %%%
 %ind2=find(mn==D(I));
 ind2=find(mn==D(I),1,'first');
 ind2=I(ind2);

The file median-average is auxiliary and helps to get the median and the average from any results. 

%%%%%%%%%%%%%%%%%%%%%%%%%%
The code for BIRECT is available at : https://data.mendeley.com/datasets/t6vv9yknbc

Any suggestions or corrections are welcome.
email me to: lchiter@univ-setif.dz