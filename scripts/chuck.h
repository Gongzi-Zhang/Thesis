using namespace std;
class table {
    private:
	map<double, vector<double>> angles;
	map<double, vector<double>> xsections;
	map<double, vector<double>> asyms;
	vector<double> energies;

    public:
	table(const char *fn="c12_fsu.dat");
	double interpolate(double E, double angle, int value=0);
};

table::table(const char * fn) {
    ifstream fin(fn, ifstream::in);
    if (!fin.is_open()) {
	cerr << "Can't open file: " << fn << endl;
	exit(4);
    }

    angles.clear();
    xsections.clear();
    asyms.clear();
    energies.clear();

    char first[20], second[20], third[20];
    double E;
    while (fin >> first) {
	if (first[0] == 'E') {
	    int i=0;
	    while (first[i++] != '=');
	    E = atof(first+i);
	    cout << "INFO\tenergies = " << E << endl;
	    energies.push_back(E);
	    fin >> first;
	}
	fin >> second >> third;
	angles[E].push_back(atof(first));
	xsections[E].push_back(atof(second));
	asyms[E].push_back(atof(third));
    }
    fin.close();
}

double table::interpolate(double E, double angle, int value) { 
    // value: 0 for asym (default); 1 for xsection
    double asym, asym1, asym2;
    double xs, xs1, xs2;
    double angle1, angle2;
    if (E < energies[0] || E > energies[energies.size()-1]) {
	cerr << "ERROR--Energy out of range: " 
	    << energies[0] << "-" << energies[energies.size()-1] 
	    << "\tyou input energy is: " << E << endl;
	return 0./0.;
    }
    int i=1;
    double e1, e2;
    while (i < energies.size() && energies[i] < E) i++;
    e1 = energies[i-1];
    e2 = energies[i];

    if (angle < angles[e1][0] || angle > angles[e1][angles[e1].size()-1]) {
	cerr << "ERROR--Angle out of range: " 
	    << angles[e1][0] << "-" << angles[e1][angles[e1].size()-1] 
	    << " for energy=" << e1
	    << "\tyou input angle is: " << angle << endl;
	return 0./0.;
    }
    if (angle < angles[e2][0] || angle > angles[e2][angles[e2].size()-1]) {
	cerr << "ERROR--Angle out of range: " 
	    << angles[e2][0] << "-" << angles[e2][angles[e2].size()-1] 
	    << " for energy=" << e2
	    << "\tyou input angle is: " << angle << endl;
	return 0./0.;
    }

    i = 1;
    while (i<angles[e1].size() && angles[e1][i] < angle)    i++;
    angle1 = angles[e1][i-1];
    angle2 = angles[e1][i];
    asym1 = asyms[e1][i-1] + (asyms[e1][i]-asyms[e1][i-1])/(angle2-angle1)*(angle-angle1);
    xs1 = xsections[e1][i-1] + (xsections[e1][i]-xsections[e1][i-1])/(angle2-angle1)*(angle-angle1);

    i = 1;
    while (i<angles[e2].size() && angles[e2][i] < angle)    i++;
    angle1 = angles[e2][i-1];
    angle2 = angles[e2][i];
    asym2 = asyms[e2][i-1] + (asyms[e2][i]-asyms[e2][i-1])/(angle2-angle1)*(angle-angle1);
    xs2 = xsections[e2][i-1] + (xsections[e2][i]-xsections[e2][i-1])/(angle2-angle1)*(angle-angle1);

    asym = asym1 + (asym2-asym1)/(e2-e1)*(E-e1);
    xs = xs1 + (xs2-xs1)/(e2-e1)*(E-e1);

    // cout << Form("OUTPUT: energy=%.4f\tangle=%.2f\txsection=%.4f\tasym=%.4f ppm", E, angle, xs, asym/1e-6) << endl;
    return (value == 0) ? asym : xs;
}
