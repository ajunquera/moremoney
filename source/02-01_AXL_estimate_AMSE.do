********************************************************************************
* ASSEGNO PER IL LAVORO - 02-01 Picking a polynomial degree for RDD
* Author: Álvaro F. Junquera
********************************************************************************

* Treatment D1

net install rdmse, from(https://raw.githubusercontent.com/peizhuan/rdmse/master) replace
use "C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script01/indi_ns_stata_D1.dta"

matrix res = J(1, 12, .) // creates the matrix to save results


 * Y1.1
rdmse post_interval6 scoringD1_0, deriv(0) c(0) p(1) h(0.0640470731480583) b(0.110211501689641) kernel(triangular)
matrix res[1,1] = r(amse_cl)

rdmse post_interval6 scoringD1_0, deriv(0) c(0) p(2) h(0.080732004637794) b(0.109571813649169) kernel(triangular)
matrix res[1,2] = r(amse_cl)

* Y1.2
rdmse post_interval712 scoringD1_0, deriv(0) c(0) p(1) h(0.0614476052733848) b(0.0936065400260332) kernel(triangular)
matrix res[1,3] = r(amse_cl)

rdmse post_interval712 scoringD1_0, deriv(0) c(0) p(2) h(0.095202145191588) b(0.127547259757979) kernel(triangular)
matrix res[1,4] = r(amse_cl)

* Y1.3
rdmse post_interval1318 scoringD1_0, deriv(0) c(0) p(1) h(0.0452136553673777) b(0.0729286568333182) kernel(triangular)
matrix res[1,5] = r(amse_cl)

rdmse post_interval1318 scoringD1_0, deriv(0) c(0) p(2) h(0.0807233617068251) b(0.110334341735329) kernel(triangular)
matrix res[1,6] = r(amse_cl)

* Y1.4
rdmse post_interval1924 scoringD1_0, deriv(0) c(0) p(1) h(0.0526641813478227) b(0.0787968615560281) kernel(triangular)
matrix res[1,7] = r(amse_cl)

rdmse post_interval1924 scoringD1_0, deriv(0) c(0) p(2) h(0.0869321287266897) b(0.11500118288894) kernel(triangular)
matrix res[1,8] = r(amse_cl)

* MECHANISMS (js hours)
rdmse jshours scoringD1_0, deriv(0) c(0) p(1) h(0.0457096546465008) b(0.0815086542260088) kernel(triangular)
matrix res[1,9] = r(amse_cl)

rdmse jshours scoringD1_0, deriv(0) c(0) p(2) h(0.0761072602485587) b(0.108717970703897) kernel(triangular)
matrix res[1,10] = r(amse_cl)

* MECHANISMS (tr hours)
rdmse attiv_form_ore_prev scoringD1_0, deriv(0) c(0) p(1) h(0.0880123802774152) b(0.132090581454506) kernel(triangular)
matrix res[1,11] = r(amse_cl)

rdmse attiv_form_ore_prev scoringD1_0, deriv(0) c(0) p(2) h(0.0672244652818837) b(0.0897278013176474) kernel(triangular)
matrix res[1,12] = r(amse_cl)

* Export to Excel
// Create a temporal dataset with the matrix
clear 
svmat res, names(col)

// Export the dataset as MS Excel
export excel using "C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script02/stata/resu_amse_1.xlsx", firstrow(variables) replace


* Treatment D2
use "C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script01/indi_ns_stata_D2.dta"

matrix resu = J(1, 12, .) // creates the matrix to save results


* Y1.1
rdmse post_interval6 scoringD2_0, deriv(0) c(0) p(1) h(0.0373618696007985) b(0.0714720564037321) kernel(triangular)
matrix resu[1,1] = r(amse_cl)

rdmse post_interval6 scoringD2_0, deriv(0) c(0) p(2) h(0.0613421179970907) b(0.0912042068600323) kernel(triangular)
matrix resu[1,2] = r(amse_cl)

* Y1.2
rdmse post_interval712 scoringD2_0, deriv(0) c(0) p(1) h(0.044258712319243) b(0.0850967616013454) kernel(triangular)
matrix resu[1,3] = r(amse_cl)

rdmse post_interval712 scoringD2_0, deriv(0) c(0) p(2) h(0.0556633846425218) b(0.0858110460502846) kernel(triangular)
matrix resu[1,4] = r(amse_cl)

* Y1.3
rdmse post_interval1318 scoringD2_0, deriv(0) c(0) p(1) h(0.0487790367230225) b(0.0865974616352447) kernel(triangular)
matrix resu[1,5] = r(amse_cl)

rdmse post_interval1318 scoringD2_0, deriv(0) c(0) p(2) h(0.0507201897416726) b(0.0795055390732065) kernel(triangular)
matrix resu[1,6] = r(amse_cl)

* Y1.4
rdmse post_interval1924 scoringD2_0, deriv(0) c(0) p(1) h(0.0438561291816327) b(0.0810442283600625) kernel(triangular)
matrix resu[1,7] = r(amse_cl)

rdmse post_interval1924 scoringD2_0, deriv(0) c(0) p(2) h(0.0467136666725013) b(0.0768616834451484) kernel(triangular)
matrix resu[1,8] = r(amse_cl)

* MECHANISMS (js hours)
rdmse jshours scoringD2_0, deriv(0) c(0) p(1) h(0.0744972275088399) b(0.115323400782257) kernel(triangular)
matrix resu[1,9] = r(amse_cl)

rdmse jshours scoringD2_0, deriv(0) c(0) p(2) h(0.075495074063369) b(0.102150590427237) kernel(triangular)
matrix resu[1,10] = r(amse_cl)

* MECHANISMS (tr hours)
rdmse attiv_form_ore_prev scoringD2_0, deriv(0) c(0) p(1) h(0.055321550828718) b(0.0855524173618022) kernel(triangular)
matrix resu[1,11] = r(amse_cl)

rdmse attiv_form_ore_prev scoringD2_0, deriv(0) c(0) p(2) h(0.097172950655242) b(0.135380873005315) kernel(triangular)
matrix resu[1,12] = r(amse_cl)

* Export to Excel
clear 
svmat resu, names(col)

export excel using "C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script02/stata/resu_amse_2.xlsx", firstrow(variables) replace

