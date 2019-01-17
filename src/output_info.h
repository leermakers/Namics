#pragma once

#include <string>

using namespace std;

class OutputInfo {
public:
	OutputInfo();

	void addProperty(const string &elem, const string &param, const string &value);
	void addFolderInfo(const string &param, const string &value);

	const string &getOutputPath() const;
	bool isOutputExists() const;

	static const std::string IN_CLASS_NAME;
	static const std::string DEFAULT_OUTPUT_PATH;

private:
	std::string output_path;

	std::string processOutputPath(std::string path);
};
