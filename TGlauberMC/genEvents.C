#include "runglauber_v3.1.C"
#include "Hist2Txt.C"

void genEvents (int n1, int n2, const char* sys1, const char* sys2) {

  //runAndSmearNucleons (n1, n2, sys1, sys2, 42, 0.4, 0.4, false, Form ("outputs/outFile_%s%s_%i_%i.root", sys1, sys2, n1, n2));
  runAndSmearNucleons (n1, n2, sys1, sys2, 72, 0.4, 0.4, false, Form ("outputs/outFile_%s%s_%i_%i.root", sys1, sys2, n1, n2));

  Hist2Txt (n1, n2, sys1, sys2);

  return;

}
