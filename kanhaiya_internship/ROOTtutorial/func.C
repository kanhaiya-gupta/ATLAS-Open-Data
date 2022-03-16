#include "myroot.h"
// example for a root macro. The standard extension of root macros is .C 
//
// root macros are implemented as functions
// the name of the function must match the filename
//
// macros are called from the root command line:  .x func.C
void func()  
{
   TF1* func= new TF1("func","sin(x)/x",0.,15.);
   func->Draw();
   c1->SaveAs("func.jpg"); // c1 is the name of ROOT's default canvas object
}


