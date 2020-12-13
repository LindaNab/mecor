## Resubmission
This is a resubmission. In this version I have:

* Described the methods implemented and added key references to the description field in DESCRIPTION
* Ommitted the use of '<<-' in ipwm
* Changed the broken doi in README

## Test environments
* local OS x86_64-apple-darwin17.0, R 4.0.2
* win-builder (devel)
* win-builder (release)

## R CMD check results

0 errors | 0 warnings | 1 notes

There was 1 NOTE:
Possibly mis-spelled words in DESCRIPTION:
  Buonaccorsi (8:522, 8:994)
  Crainiceanu (8:144, 8:1093)
  Fieller (8:966)
  JA (8:624)
  JW (8:376)
  Kipnis (8:247)
  RJ (8:114, 8:242, 8:614, 8:1063)
  Rosner (8:789, 8:877)
  Ruppert (8:118, 8:1067)
  Spiegelman (8:220, 8:799, 8:887)
  Stavola (8:380)
  Stefanski (8:129, 8:1078)
  Tooze (8:618)
  Willett (8:814, 8:902)
These are names and not misspelled.

Found the following (possibly) invalid DOIs:
  DOI: 10.1002/1097-0258(20010115)20:1<139::AID-SIM644
    From: DESCRIPTION
    Status: Not Found
    Message: 404
The doi contains a closing angle bracket, which causes the note to occur: https://doi.org/10.1002/1097-0258(20010115)20:1<139::AID-SIM644>3.0.CO;2-K
I tried different things but I don't know how I should solve this problem.

## revdepchecks
No reverse dependencies

