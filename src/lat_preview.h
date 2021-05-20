#ifndef LAT_PREVIEWxH
#define LAT_PREVIEWxH

#include "namics.h"
#include "input.h"

class Lat_preview {
public:
	Lat_preview(vector<Input*>,string);

virtual ~Lat_preview();

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;
	string name;
	vector<Input*> In;


	int gradients;
	string geometry;

	bool CheckInput(int start);
	string GetValue(string);

};
#endif
