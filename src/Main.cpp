#include "Global.h"

int main(int argc, char *argv[])

{
	Initialize(argc, argv);
	Read();
	Prepare();
	// CalPec();
	Moncar();
	Output();
	return 0;
}
/*ofstream oup("doc/cross_mayue.txt");
	 for (double i = 0.1; i < 100; i += 0.1)
	 {
		 oup << i << "\t" << SCD12_H.cal(i, 1e18) << "\t" << ACD12_H.cal(i, 1e18) << "\t" << R2_2_12_H4.cal(1e18, i) << "\t" << R3_2_3_H3.cal(i, i) << endl;
	 }
	oup.close();*/