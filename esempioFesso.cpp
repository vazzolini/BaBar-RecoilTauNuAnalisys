#include "TFile.h"
//#include <iostream.h>
void main()
{
  //  TFile f("a");
  TFile * pippo = new TFile("provola.root", "RECREATE");
  //  cout<< "mi chiamo " << pippo->GetName()<< endl;
  pippo->cd();
}
