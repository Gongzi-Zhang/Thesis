const double pi = 3.14159265358979323846;
const double e = 2.71828182845904523536;
const double MeV = 1e6;
const double GeV = 1e9;
const double fm = 1e-15;
const double E2L = 5.068/fm/GeV;    // 1 GeV = 5.068 fm⁻¹
const double barn = 1e-28;  // 1 b = 1e-28 m²
const TF1 * d2r = new TF1("d2r", "x/180*3.14159265358979323846", 0, 360);
void FF()
{
    const int npoints = 100;
    const double xmin = 15;
    const double xmax = 60;
    const double E = 757.5 * MeV;
    const double R = pow(48, 1/3.) * fm;
    TF1 * f = new TF1("f", "1e4/(137*137*4*[0]*[0])*pow(cos(d2r(x)/2), 2)/pow(sin(d2r(x)/2), 7)*pow(3/pow([1]*sin(d2r(x)/2), 3)*(sin([1]*sin(d2r(x)/2)) - [1]*cos([1]*sin(d2r(x)/2))), 2)", xmin, xmax);
    // TF1 * f = new TF1("f", "1/pow(cos(d2r(x)/2), 2)/pow(sin(d2r(x))/2, 7)*pow(3/pow([0]*sin(d2r(x)/2), 3)*(sin([0]*sin(d2r(x)/2)) - [0]*cos([0]*sin(d2r(x)/2))), 2)", xmin, xmax);
    f->SetParameters(E*E2L, 2*E*R*E2L);
    f->SetTitle("d#sigma/d#Omega (cm^{2}/str);deg");
    f->SetNpx(1000);
    f->SetLineColor(kBlack);
    // TF1 * f = new TF1("f", "1/pow(sin(x/2), 3)*(sin(5.5*sin(x/2)) - 5.5*cos(5.5*sin(x/2)))", xmin, xmax);
    TCanvas c("c", "c", 600, 800);
    c.SetLogy(1);
    f->Draw();
    c.SaveAs("FF.png");
}

void mott()
{
    const int npoints = 100;
    const double xmin = 15;
    const double xmax = 60;
    const double E = 757.5 * MeV;
    const double R = pow(48, 1/3.) * fm;
    TF1 * f = new TF1("f", "1e4/(137*137*4*[0]*[0])*pow(cos(d2r(x)/2), 2)/pow(sin(d2r(x)/2), 4)", xmin, xmax);
    f->SetParameter(0, E*E2L);
    f->SetTitle("d#sigma/d#Omega (cm^{2}/str);deg");
    f->SetNpx(1000);
    f->SetLineColor(kBlack);
    // TF1 * f = new TF1("f", "1/pow(sin(x/2), 3)*(sin(5.5*sin(x/2)) - 5.5*cos(5.5*sin(x/2)))", xmin, xmax);
    TCanvas c("c", "c", 600, 800);
    c.SetLogy(1);
    f->Draw();
    c.SaveAs("mott.png");
}

void fermi_dist()
{
    const double xmin = 0;
    const double xmax = 10;
    TF1 *f = new TF1("f", "[0]/(1+pow(2.7182818284590452353602874713, (x-[1])/[2]))", xmin, xmax);
    const double rho = 1;
    const double c = 6.5;
    const double a = 0.5;
    f->SetParameters(rho, c, a);
    TCanvas C("c", "c", 800, 600);
    f->Draw();
    f->SetTitle("Saxon-Woods Distribution;r (fm);#rho");
    f->GetXaxis()->SetNdivisions(10);
    f->GetYaxis()->SetNdivisions(3);
    TArrow *ar = new TArrow(0, 0.5, c, 0.5, 0.02, "<|>");
    ar->SetAngle(30);
    ar->SetLineWidth(2);
    ar->Draw();
    TText *t = new TText(c/2, 0.52, "c");
    t->Draw();
    C.SaveAs("fermi.png");
}
void plot_function()
{
    // mott();
    // FF();
    fermi_dist();
}
