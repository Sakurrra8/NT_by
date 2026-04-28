#include "Global.h"

int main(int argc, char *argv[])

{
	Initialize(argc, argv);
	Read();
	Prepare();

	Moncar();
	Output();

	return 0;
}
// CalPec();
/*ofstream oup("doc/cross_CX_H_1e18.txt");
	double n___ = 1e18;
	for (double i = 0.1; i < 100; i += 0.1)
	{
		oup << i << "\t" << CCD96_H.cal(i, n___) << endl;
		// oup << i << "\t" << R3_2_3r_H4.cal(n___, i) << "\t" << R3_2_3d_H4.cal(n___, i) << "\t" << R3_2_3i_H4.cal(n___, i) << "\t" << R2_2_17r_H4.cal(n___, i) << "\t" << R2_2_17d_H4.cal(n___, i) << endl;
	}
	oup.close();

	ofstream oup("doc/cross-section_CX_n-n.txt");
	double n___ = 1e18;
	for (double i = 0.1; i < 100; i += 0.1)
	{
		oup << i << "\t" << R_H_H.cal(n___, i) << "\t" << R_H_H2.cal(n___, i) << "\t" << R_H2_H.cal(n___, i) << "\t" << R_H2_H2.cal(n___, i) << endl;
		// oup << i << "\t" << R3_2_3r_H4.cal(n___, i) << "\t" << R3_2_3d_H4.cal(n___, i) << "\t" << R3_2_3i_H4.cal(n___, i) << "\t" << R2_2_17r_H4.cal(n___, i) << "\t" << R2_2_17d_H4.cal(n___, i) << endl;
	}
	oup.close();*/