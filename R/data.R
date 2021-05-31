#' PDAC blood pressure data [replicates study]
#'
#' Blood pressure, age and creatinine levels of 450 pregnant women from the
#' Pregnancy Day Assessment Clinic.
#'
#' This is a simulated dataset inspired by data that was originally published at the Dryad Digital Repository: <doi:10.5061/dryad.0bq15>
#'
#' @format A data frame with 450 rows and 6 variables:
#' \describe{
#'   \item{creatinine}{Serum creatinine (umol/L)}
#'   \item{age}{Age (years)}
#'   \item{sbp30}{Systolic blood pressure at 30 minutes (mm Hg)}
#'   \item{sbp60}{Systolic blood pressure at 60 minutes (mm Hg)}
#'   \item{sbp90}{Systolic blood pressure at 90 minutes (mm Hg)}
#'   \item{sbp120}{Systolic blood pressure at 120 minutes (mm Hg)}
#' }
#' @examples
#' data("bloodpressure", package = "mecor")
#' @references
#' Elizabeth Anne McCarthy, Thomas A Carins, Yolanda Hannigan, Nadia Bardien, Alexis Shub, and Susan P Walker. Data from: Effectiveness and safety of 1 vs 4h blood pressure profile with clinical and laboratory assessment for the exclusion of gestational hypertension and pre-eclampsia: a retrospective study in a university affiliated maternity hospital. Dryad (2015). <doi:10.5061/dryad.0bq15>.
"bloodpressure"
#' NEO visceral adipose tissue data [internal covariate-validation study]
#'
#' Insulin resistance, waist circumference, sex, age, total body fat and visceral adipose tissue of 650 individuals from the
#' Netherlands Epidemiology of Obesity (NEO) study. Visceral adipose tissue measurements were taken of approximately 40\% of the individuals, at random.
#'
#' This is a simulated data set inspired by the NEO data <doi:10.1007/s10654-013-9801-3>. A motivating example using the example data can be found here: <doi:10.1093/aje/kwab114>
#'
#' @format A data frame with 650 rows and 6 variables:
#' \describe{
#'   \item{ir_ln}{Natural logarithm of insulin resistance (fasting glucose (mmol/L) x fasting insulin (mU/L) / 22.5)}
#'   \item{wc}{Waist circumference (standardised, cm)}
#'   \item{sex}{Sex (0 = male, 1 = female)}
#'   \item{age}{Age (years)}
#'   \item{tbf}{Total body fat (standardised, \%)}
#'   \item{vat}{Visceral adipose tissue (standardised, cm^2)}
#' }
#' @examples
#' data("vat", package = "mecor")
#' @references
#' Renee de Mutsert, Martin den Heijer, Ton J Rabelink, Johannes WA Smit, Johannes A Romijn, Johan W Jukema, Albert de Roos, Christa M Cobbaert, Margreet Kloppenburg, Saskia le Cessie, Saskia Middeldorp, Frits R Rosendaal. The Netherlands epidemiology of obesity (NEO) study: Study design and data collection. European Journal of Epidemiology (2013). <doi:10.1007/s10654-013-9801-3>
#'
#' Linda Nab, Maarten van Smeden, Renee de Mutsert, Frits R Rosendaal, and Rolf HH Groenwold. Sampling strategies for internal validation samples for exposure measurement error correction: A study of visceral adipose tissue measures replaced by waist circumference measures. American Journal of Epidemiology (2021). <doi:10.1093/aje/kwab114>
"vat"
#' Visceral adipose tissue external data [external covariate-validation study]
#'
#' Waist circumference, visceral adipose tissue, sex, age, and total body fat of 100 individuals
#'
#' This is a simulated data set accompanying the dataset "vat", that is inspired by the NEO data <doi:10.1007/s10654-013-9801-3>. A motivating example using the example data can be found here: <doi:10.1093/aje/kwab114>
#'
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'   \item{wc}{Waist circumference (standardised, cm)}
#'   \item{vat}{Visceral adipose tissue (standardised, cm^2)}
#'   \item{sex}{Sex (0 = male, 1 = female)}
#'   \item{age}{Age (years)}
#'   \item{tbf}{Total body fat (standardised, \%)}
#' }
#' @examples
#' data("vat_ext", package = "mecor")
#' @references
#' Renee de Mutsert, Martin den Heijer, Ton J Rabelink, Johannes WA Smit, Johannes A Romijn, Johan W Jukema, Albert de Roos, Christa M Cobbaert, Margreet Kloppenburg, Saskia le Cessie, Saskia Middeldorp, Frits R Rosendaal. The Netherlands epidemiology of obesity (NEO) study: Study design and data collection. European Journal of Epidemiology (2013). <doi:10.1007/s10654-013-9801-3>
#'
#' Linda Nab, Maarten van Smeden, Renee de Mutsert, Frits R Rosendaal, and Rolf HH Groenwold. Sampling strategies for internal validation samples for exposure measurement error correction: A study of visceral adipose tissue measures replaced by waist circumference measures. American Journal of Epidemiology (2021). <doi:10.1093/aje/kwab114>
"vat_ext"
#' Low-dose iron supplements haemoglobin data [internal outcome-validation study]
#'
#' Capillary haemoglobin and venous haemoglobin levels of 400 subjects of a trial investigating the efficacy of low-dose iron supplements during pregnancy.
#' Venous haemoglobin levels were observed of approximately 25\% of the subjects included in the trial.
#'
#' This is a simulated data set inspired by a trial investigating low-dose iron supplements <doi:10.1093/ajcn/78.1.145>. A motivating example using the example data can be found here: <doi:10.1002/sim.8359>
#'
#' @format A data frame with 400 rows and 3 variables:
#' \describe{
#'   \item{capillary}{Haemoglobin levels measured in capillary blood (g/L)}
#'   \item{supplement}{Low-dose iron supplement (20 mg/d) (0 = no, 1 = yes)}
#'   \item{venous}{Haemoglobin levels measured in venous blood (g/L)}
#' }
#' @examples
#' data("haemoglobin", package = "mecor")
#' @references
#' Maria Makrides, Caroline A Crowther, Robert A Gibson, Rosalind S Gibson, and C Murray Skeaff. Efficacy and tolerability of low-dose iron supplements during pregnancy: a randomized controlled trial. The American Journal of Clinical Nutrition (2003). <doi:10.1093/ajcn/78.1.145>
#'
#' Linda Nab, Rolf HH Groenwold, Paco MJ Welsing, and Maarten van Smeden. Measurement error in continuous endpoints in randomised trials: Problems and solutions. Statistics in Medicine (2019). <doi:10.1002/sim.8359>
"haemoglobin"
#' Haemoglobin external data [external outcome-validation study]
#'
#' Capillary haemoglobin and venous haemoglobin levels of 100 individuals.
#'
#' This is a simulated data set accompanying the dataset "haemoglobin", that is inspired by a trial investigating low-dose iron supplements <doi:10.1093/ajcn/78.1.145>. A motivating example using the example data can be found here: <doi:10.1002/sim.8359>
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{capillary}{Haemoglobin levels measured in capillary blood (g/L)}
#'   \item{venous}{Haemoglobin levels measured in venous blood (g/L)}
#' }
#' @examples
#' data("haemoglobin_ext", package = "mecor")
#' @references
#' Maria Makrides, Caroline A Crowther, Robert A Gibson, Rosalind S Gibson, and C Murray Skeaff. Efficacy and tolerability of low-dose iron supplements during pregnancy: a randomized controlled trial. The American Journal of Clinical Nutrition (2003). <doi:10.1093/ajcn/78.1.145>
#'
#' Linda Nab, Rolf HH Groenwold, Paco MJ Welsing, and Maarten van Smeden. Measurement error in continuous endpoints in randomised trials: Problems and solutions. Statistics in Medicine (2019). <doi:10.1002/sim.8359>
"haemoglobin_ext"
#' TONE sodium data [outcome-calibration study]
#'
#' Self-reported sodium intake and urinary sodium in the TONE study, a randomized controlled trial designed to
#' investigate whether a reduction in sodium intake results in satisfactory blood pressure control.
#' Two replicate urinary sodium measures were available in 50\% of the subjects included in the trial.
#'
#' This is a simulated data set inspired by the TONE study <doi: 10.1016/1047-2797(94)00056-y>. A motivating example using the example data can be found here: <doi:10.1002/sim.7011>
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{recall}{Sodium intake measured by a 24h recall (mg)}
#'   \item{diet}{Usual diet or sodium-lowering diet (0 = usual, 1 = sodium-lowering)}
#'   \item{urinary1}{Sodium intake measured in urine (1st measure, mg)}
#'   \item{urinary2}{Sodium intake measured in urine (2nd measure, mg)}
#' }
#' @examples
#' data("sodium", package = "mecor")
#' @references
#' Lawrence J Appel, Mark Espeland, Paul K Whelton, Therese Dolecek, Shiriki Kumanyika, William B Applegate, Walter H Ettinger, John B Kostis, Alan C Wilson, Clifton Lacy, and Stephen T Miller. Trial of Nonpharmacologic Intervention in the Elderly (TONE). Design and rationale of a blood pressure control trial. Annals of Epidemiology (1995). <doi: 10.1016/1047-2797(94)00056-y>
#'
#' Ruth H Keogh, Raymond J Carroll, Janet A Tooze, Sharon I Kirkpatrick, Laurence S Freedman. Statistical issues related to dietary intake as the response variable in intervention trials. Statistics in Medicine (2016). <doi:10.1002/sim.7011>
"sodium"
#' @references
#' Maria Makrides, Caroline A Crowther, Robert A Gibson, Rosalind S Gibson, and C Murray Skeaff. Efficacy and tolerability of low-dose iron supplements during pregnancy: a randomized controlled trial. The American Journal of Clinical Nutrition (2003). <doi:10.1093/ajcn/78.1.145>
#'
#' Linda Nab, Rolf HH Groenwold, Paco MJ Welsing, and Maarten van Smeden. Measurement error in continuous endpoints in randomised trials: Problems and solutions. Statistics in Medicine (2019). <doi:10.1002/sim.8359>
"haemoglobin"
#' Simulated dataset for the \link[mecor]{ipwm} function
#'
#' A simulated dataset containing 5000 observations of the covariates L1-L10,
#' the true exposure A and true outcome Y, and the misclassified exposure B and
#' misclassified outcome Z.
#'
#' @format A data frame with 5000 rows and 14 variables:
#' \describe{
#'   \item{L1}{covariate, binary}
#'   \item{L2}{covariate, continuous}
#'   \item{L3}{covariate, binary}
#'   \item{L4}{covariate, continuous}
#'   \item{L5}{covariate, binary}
#'   \item{L6}{covariate, binary}
#'   \item{L7}{covariate, continuous}
#'   \item{L8}{covariate, binary}
#'   \item{L9}{covariate, binary}
#'   \item{L10}{covariate, continuous}
#'   \item{A}{exposure, binary}
#'   \item{Y}{outcome, binary}
#'   \item{B}{misclassified exposure, binary}
#'   \item{Z}{misclassified outcome, binary}
#' }
#' @examples
#' data("sim", package = "mecor")
"sim"

