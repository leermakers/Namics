#pragma once

#include <string>

class OutputInfo {
public:
	OutputInfo();

	void addProperty(const std::string &elem, const std::string &param, const std::string &value);

	void addFolderInfo(const std::string &param, const std::string &value);

	const std::string &getOutputPath() const;

	static const std::string IN_CLASS_NAME;
	static const std::string DEFAULT_OUTPUT_PATH;

private:
	std::string output_path;

	std::string processOutputPath(std::string path);
};
