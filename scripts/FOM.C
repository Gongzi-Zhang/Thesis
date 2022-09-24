#include "chuck.h"
#include "constant.h"

using namespace std;

/*
class PVES {
  private:
    double A = 208;
    double rho = 11.38*g/cm3; 
    double thickness = 0.55*mm;
    double E = 950*MeV;
    double theta = 5*deg; 
    double I = 70*uA;
    double dOmega = 0.0037;
    double acceptance = 0.4;
    double eta = 0.5;
    double rate_coef = (I/e_charge)*(rho/A)*thickness*NA*dOmega*acceptance*eta*1e-27/MHz;
    table *t0, *t1;
    double emin = 550*MeV;
    double emax = 2000*MeV;
    double thetamin = 2*deg;
    double thetamax = 15*deg;
  public:
    PVES(string n0, string n1) {
	t0 = new table(n0.c_str());
	t1 = new table(n1.c_str());
    };
    void set_A(const double val) {A = val;};
    void set_rho(const double val) {rho = val;};
    void set_thickness(const double val) {thickness = val;};
    void set_E(const double val) {E = val;};
    void set_theta(const double val) {theta = val;};
    void set_I(const double val) {I = val;};
    void set_dOmega(const double val) {dOmega = val;};
    void set_acceptance(const double val) {acceptance = val;};
    void set_eta(const double val) {eta = val;};
    void set_rate_coef() {rate_coef = (I/e_charge) * (rho/(g/cm3)/A) * (thickness/cm) * NA * dOmega * acceptance * (mbarn/cm2) /MHz;};
    void set_emin(const double val) {emin = val;};
    void set_emax(const double val) {emax = val;};
    void set_thetamin(const double val) {thetamin = val;};
    void set_thetamax(const double val) {thetamax = val;};
    void plot_xsec() {};
};
*/
void pb208()
{
    table pb("horpb.dat");
    table pb1("horpb_stretched.dat");
    const double A = 208;
    double rho = 11.38*g/cm3;
    double thickness = 0.55*mm;
    double E = 950*MeV;
    const double theta = 5*deg;
    const double I = 70*uA;
    const double dOmega = 0.0037;
    // const double acceptance = 0.4;
    const double eta = 0.5;

    rho /= (g/cm3);
    thickness /= cm;
    E /= MeV;
    const double rate_coef = (I/e_charge)*(rho/A)*thickness*NA*dOmega*eta*(mbarn/cm2)/MHz;

    double xsec, asym, asym1, rate, sen, fom;

    TGraph *g_asym_e = new TGraph();
    TGraph *g_rate_e = new TGraph();
    TGraph *g_sen_e = new TGraph();
    TGraph *g_fom_e = new TGraph();
    int i=0;
    for (double e=550; e<=2000; e+=50)
    {
	xsec = pb.interpolate(e, theta, 1);  // mb
	asym = pb.interpolate(e, theta, 0)/ppm;  
	asym1 = pb1.interpolate(e, theta, 0)/ppm;  
	rate = xsec*rate_coef;

	g_asym_e->SetPoint(i, e, asym);
	g_rate_e->SetPoint(i, e, rate);
	if (!isnan(asym1))
	{
	    sen = abs(asym1/asym-1)*100;
	    g_sen_e->SetPoint(i, e, sen);
	    fom = rate * asym*asym * sen*sen * 1e-3;
	    g_fom_e->SetPoint(i, e, fom);
	}
	i++;
    }

    TGraph *g_asym_theta = new TGraph();
    TGraph *g_rate_theta = new TGraph();
    TGraph *g_sen_theta = new TGraph();
    TGraph *g_fom_theta = new TGraph();
    i = 0;
    for (double t=2; t<=15; t+=0.2)
    {
	xsec = pb.interpolate(E, t, 1);  
	asym = pb.interpolate(E, t, 0)/ppm;
	asym1 = pb1.interpolate(E, t, 0)/ppm;  
	rate = rate_coef*xsec;

	g_asym_theta->SetPoint(i, t, asym);
	g_rate_theta->SetPoint(i, t, rate);
	if (!isnan(asym1))
	{
	    sen = abs(asym1/asym-1)*100;
	    g_sen_theta->SetPoint(i, t, sen);
	    fom = rate * asym*asym * sen*sen * 1e-3;
	    g_fom_theta->SetPoint(i, t, fom);
	}
	i++;
    }

    TCanvas *c = new TCanvas("c", "c", 1600, 600);
    c->Divide(2, 1);

    c->cd(1);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_asym_e->SetTitle("Pb208 (#theta = 5#circ); E (MeV); asym (ppm)");
    g_asym_e->GetXaxis()->CenterTitle(true);
    g_asym_e->GetYaxis()->CenterTitle(true);
    g_asym_e->SetMarkerStyle(20);
    g_asym_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_asym_theta->SetTitle("Pb208 (E = 950 MeV); #theta (deg); asym (ppm)");
    g_asym_theta->GetXaxis()->CenterTitle(true);
    g_asym_theta->GetYaxis()->CenterTitle(true);
    g_asym_theta->SetMarkerStyle(20);
    g_asym_theta->Draw("AP");

    c->SaveAs("Pb208_asym.png");

    c->Clear();
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_rate_e->SetTitle("Pb208 (#theta = 5#circ); E (MeV); rate (MHz)");
    g_rate_e->GetXaxis()->CenterTitle(true);
    g_rate_e->GetYaxis()->CenterTitle(true);
    g_rate_e->SetMarkerStyle(20);
    g_rate_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_rate_theta->SetTitle("Pb208 (E = 950 MeV); #theta (deg); rate (MHz)");
    g_rate_theta->GetXaxis()->CenterTitle(true);
    g_rate_theta->GetYaxis()->CenterTitle(true);
    g_rate_theta->SetMarkerStyle(20);
    g_rate_theta->Draw("AP");

    c->SaveAs("Pb208_rate.png");

    c->Clear();
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetGrid();
    g_sen_e->SetTitle("Pb208 (#theta = 5#circ); E (MeV); dA/A (%)");
    g_sen_e->GetXaxis()->CenterTitle(true);
    g_sen_e->GetYaxis()->CenterTitle(true);
    g_sen_e->SetMarkerStyle(20);
    g_sen_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    g_sen_theta->SetTitle("Pb208 (E = 950 MeV); #theta (deg); dA/A (%)");
    g_sen_theta->GetXaxis()->CenterTitle(true);
    g_sen_theta->GetYaxis()->CenterTitle(true);
    g_sen_theta->SetMarkerStyle(20);
    g_sen_theta->Draw("AP");

    c->SaveAs("Pb208_sen.png");

    c->Clear();
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetGrid();
    g_fom_e->SetTitle("Pb208 (#theta = 5#circ); E (MeV); FOM");
    g_fom_e->GetXaxis()->CenterTitle(true);
    g_fom_e->GetYaxis()->CenterTitle(true);
    g_fom_e->SetMarkerStyle(20);
    g_fom_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    g_fom_theta->SetTitle("Pb208 (E = 950 MeV); #theta (deg); FOM");
    g_fom_theta->GetXaxis()->CenterTitle(true);
    g_fom_theta->GetYaxis()->CenterTitle(true);
    g_fom_theta->SetMarkerStyle(20);
    g_fom_theta->Draw("AP");

    c->SaveAs("Pb208_fom.png");
}

void ca48()
{
    table ca("ca48_fsu.dat");
    table ca1("ca48_fsu_stretched.dat");
    const double A = 48;
    double rho = 1.855*g/cm3;
    double thickness = 6*mm;
    double E = 2200*MeV;
    const double theta = 4*deg;
    const double I = 150*uA;
    const double dOmega = 0.0037;
    // const double acceptance = 0.4;
    const double eta = 0.5;

    rho /= (g/cm3);
    thickness /= cm;
    E /= MeV;
    const double rate_coef = (I/e_charge)*(rho/A)*thickness*NA*dOmega*eta*(mbarn/cm2)/MHz;

    double xsec, asym, asym1, rate, sen, fom;

    TGraph *g_asym_e = new TGraph();
    TGraph *g_rate_e = new TGraph();
    TGraph *g_sen_e = new TGraph();
    TGraph *g_fom_e = new TGraph();
    int i=0;
    for (double e=800; e<=3000; e+=50)
    {
	xsec = ca.interpolate(e, theta, 1);  // mb
	asym = ca.interpolate(e, theta, 0)/ppm;  
	asym1 = ca1.interpolate(e, theta, 0)/ppm;  
	rate = xsec*rate_coef;

	g_asym_e->SetPoint(i, e, asym);
	g_rate_e->SetPoint(i, e, rate);
	if (!isnan(asym1))
	{
	    sen = abs(asym1/asym-1)*100;
	    g_sen_e->SetPoint(i, e, sen);
	    fom = rate * asym*asym * sen*sen * 1e-3;
	    g_fom_e->SetPoint(i, e, fom);
	}
	i++;
    }

    TGraph *g_asym_theta = new TGraph();
    TGraph *g_rate_theta = new TGraph();
    TGraph *g_sen_theta = new TGraph();
    TGraph *g_fom_theta = new TGraph();
    i = 0;
    for (double t=2; t<=15; t+=0.2)
    {
	xsec = ca.interpolate(E, t, 1);  
	asym = ca.interpolate(E, t, 0)/ppm;
	asym1 = ca1.interpolate(E, t, 0)/ppm;  
	rate = rate_coef*xsec;

	g_asym_theta->SetPoint(i, t, asym);
	g_rate_theta->SetPoint(i, t, rate);
	if (!isnan(asym1))
	{
	    sen = abs(asym1/asym-1)*100;
	    g_sen_theta->SetPoint(i, t, sen);
	    fom = rate * asym*asym * sen*sen * 1e-3;
	    g_fom_theta->SetPoint(i, t, fom);
	}
	i++;
    }

    TCanvas *c = new TCanvas("c", "c", 1600, 600);
    c->Divide(2, 1);

    c->cd(1);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_asym_e->SetTitle("Ca48 (#theta = 4#circ); E (MeV); asym (ppm)");
    g_asym_e->GetXaxis()->CenterTitle(true);
    g_asym_e->GetYaxis()->CenterTitle(true);
    g_asym_e->SetMarkerStyle(20);
    g_asym_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_asym_theta->SetTitle("Ca48 (E = 2200 MeV); #theta (deg); asym (ppm)");
    g_asym_theta->GetXaxis()->CenterTitle(true);
    g_asym_theta->GetYaxis()->CenterTitle(true);
    g_asym_theta->SetMarkerStyle(20);
    g_asym_theta->Draw("AP");

    c->SaveAs("Ca48_asym.png");

    c->Clear();
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_rate_e->SetTitle("Ca48 (#theta = 4#circ); E (MeV); rate (MHz)");
    g_rate_e->GetXaxis()->CenterTitle(true);
    g_rate_e->GetYaxis()->CenterTitle(true);
    g_rate_e->SetMarkerStyle(20);
    g_rate_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    gPad->SetLogy(1);
    g_rate_theta->SetTitle("Ca48 (E = 2200 MeV); #theta (deg); rate (MHz)");
    g_rate_theta->GetXaxis()->CenterTitle(true);
    g_rate_theta->GetYaxis()->CenterTitle(true);
    g_rate_theta->SetMarkerStyle(20);
    g_rate_theta->Draw("AP");

    c->SaveAs("Ca48_rate.png");

    c->Clear();
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetGrid();
    g_sen_e->SetTitle("Ca48 (#theta = 4#circ); E (MeV); dA/A (%)");
    g_sen_e->GetXaxis()->CenterTitle(true);
    g_sen_e->GetYaxis()->CenterTitle(true);
    g_sen_e->SetMarkerStyle(20);
    g_sen_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    g_sen_theta->SetTitle("Ca48 (E = 2200 MeV); #theta (deg); dA/A (%)");
    g_sen_theta->GetXaxis()->CenterTitle(true);
    g_sen_theta->GetYaxis()->CenterTitle(true);
    g_sen_theta->SetMarkerStyle(20);
    g_sen_theta->Draw("AP");

    c->SaveAs("Ca48_sen.png");

    c->Clear();
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetGrid();
    g_fom_e->SetTitle("Ca48 (#theta = 4#circ); E (MeV); FOM");
    g_fom_e->GetXaxis()->CenterTitle(true);
    g_fom_e->GetYaxis()->CenterTitle(true);
    g_fom_e->SetMarkerStyle(20);
    g_fom_e->Draw("AP");

    c->cd(2);
    gPad->SetGrid();
    g_fom_theta->SetTitle("Ca48 (E = 2200 MeV); #theta (deg); FOM");
    g_fom_theta->GetXaxis()->CenterTitle(true);
    g_fom_theta->GetYaxis()->CenterTitle(true);
    g_fom_theta->SetMarkerStyle(20);
    g_fom_theta->Draw("AP");

    c->SaveAs("Ca48_fom.png");
}
void FOM()
{
    pb208();
    ca48();
}
