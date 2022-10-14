#include "constant.h"

void compton()
{
    gROOT->SetBatch(1);
    double k = 2.334*eV;
    double me = 0.5*MeV;
    
    // double a = 1/(1 + 4*k*(E/me)/me);
    TF1 * a = new TF1("a", "1/(1+4*[0]*x/[1]/[1])", 1*MeV, 10*GeV);
    a->SetParameters(k, me);

    // double E = 1.165*GeV;
    double Ep = 0.95*GeV;
    TF1 * f1p = new TF1("f1p", "x*(1-[0])", 0, 1);
    f1p->SetParameter(0, a->Eval(Ep));
    TF1 * f2p = new TF1("f2p", "x*(1+[0])", 0, 1);
    f2p->SetParameter(0, a->Eval(Ep));
    TF1 * fp = new TF1("fp", "(1-f2p(x)) * (1 - 1/pow(1-f1p(x), 2)) / (pow(f1p(x), 2)/(1-f1p(x)) + 1 + pow((1-f2p(x))/(1-f1p(x)), 2)) ", 0, 1);
    fp->SetTitle("Analyzing power of Compton Scattering; #rho = E^{#gamma}/E^{#gamma}_{max}; A_{l}");
    fp->SetLineColor(kBlack);

    double Ec = 2.2*GeV;
    TF1 * f1c = new TF1("f1c", "x*(1-[0])", 0, 1);
    f1c->SetParameter(0, a->Eval(Ec));
    TF1 * f2c = new TF1("f2c", "x*(1+[0])", 0, 1);
    f2c->SetParameter(0, a->Eval(Ec));
    TF1 * fc = new TF1("fc", "(1-f2c(x)) * (1 - 1/pow(1-f1c(x), 2)) / (pow(f1c(x), 2)/(1-f1c(x)) + 1 + pow((1-f2c(x))/(1-f1c(x)), 2)) ", 0, 1);
    fc->SetTitle("Analyzing power of Compton Scattering; #rho = E^{#gamma}/E^{#gamma}_{max}; A_{l}");
    fc->GetYaxis()->SetTitleOffset(1.2);

    TCanvas c("c", "c", 800, 600);
    fc->Draw("L");
    fp->Draw("same");

    TLegend lg(0.2, 0.65, 0.4, 0.8);
    lg.AddEntry(fc, "E_{beam} = 2.2 GeV", "l");
    lg.AddEntry(fp, "E_{beam} = 0.95 GeV", "l");
    lg.Draw();

    TLine line(0, 0, 1, 0);
    line.SetLineStyle(2);
    line.Draw();

    c.SaveAs("compton_asym.png");
}
