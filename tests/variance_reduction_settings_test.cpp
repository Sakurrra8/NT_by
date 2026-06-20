#include "Global.h"

#include <iostream>

int main()
{
	char program[] = "variance_reduction_settings_test";
	char setting[] = "Inputfile/settingfile/setting_Trimesh_D_5_vr_on.log";
	char *arguments[] = {program, setting};
	Initialize(2, arguments);

	if (K_Roulette != 1 || K_Splitting != 1 ||
		W_RouletteMin != 0.05 || W_RouletteTarget != 0.2 ||
		W_SplitMax != 2.0 || W_SplitTarget != 1.0 ||
		W_SplitMin != 0.05 || MaxSplit != 4 || MaxSplitDepth != 1 ||
		RegionImportance[0] != 1.0 || RegionImportance[1] != 2.0 ||
		RegionImportance[2] != 1.0 ||
		ImportanceMainPoloidalBegin != 25 || ImportanceMainPoloidalEnd != 72)
	{
		std::cerr << "variance-reduction settings were not parsed correctly\n";
		return 1;
	}

	std::cout << "variance reduction settings test passed\n";
	return 0;
}
