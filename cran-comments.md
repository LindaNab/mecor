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
  Buonaccorsi (8:543, 8:1015)
  Crainiceanu (8:144, 8:1114)
  Fieller (8:987)
  JA (8:645)
  JW (8:397)
  Kipnis (8:247)
  RJ (8:114, 8:242, 8:635, 8:1084)
  Rosner (8:810, 8:898)
  Ruppert (8:118, 8:1088)
  Spiegelman (8:220, 8:820, 8:908)
  Stavola (8:401)
  Stefanski (8:129, 8:1099)
  Tooze (8:639)
  URLencode (8:263)
  Willett (8:835, 8:923)

These are names and not misspelled.

The Description field contains
  URLencode("https://doi.org/10.1002/1097-0258(20010115)20:1<139::AID-SIM644>3.0.CO;2-K")
Please enclose URLs in angle brackets (<...>).
The Description field contains
  URLencode("https://doi.org/10.1002/1097-0258(20010115)20:1<139::AID-SIM644>3.0.CO;2-K")
Please write DOIs as <doi:10.prefix/suffix>.

I was asked to use URLencode() since the doi contains angle brackets. 

## revdepchecks
No reverse dependencies

