#include "file_reader.h"

constexpr uint8_t NUM_PRO_HEADER_TOKENS = 3;
constexpr uint8_t DENSITY_HEADER_TOKEN_INDEX = 2;
constexpr uint8_t BOUNDARIES = 2;
constexpr uint8_t SYSTEM_EDGE_OFFSET = 1;
constexpr uint8_t X_DIMENSION = 0;
constexpr uint8_t Y_DIMENSION = 1;
constexpr uint8_t Z_DIMENSION = 2;
constexpr uint8_t OFFSET_VTK_DIMENSIONS_TAG = 1;

std::map<Readable_filetype, std::string> Readable_file::extension_map  {
            {Readable_filetype::NONE, ""},
            {Readable_filetype::VTK_STRUCTURED_GRID, "vtk"},
            {Readable_filetype::PRO, "pro"}
};



Readable_filetype Readable_file::get_filetype()
{
    return m_filetype;
}

//Accepts filename,
Readable_file::Readable_file(const std::string filename_, Readable_filetype filetype_)
    : m_filename{filename_}, m_filetype{filetype_}
{
    try
    {
        check_filetype();
    }
    catch (std::string extension)
    {
        std::cerr << "Extension: " << extension << " not recognized." << std::endl;
        exit(error::ERROR_EXTENSION);
    }

    if (access(m_filename.c_str(), F_OK) == -1)
    {
        std::cerr << "Error opening file! Is the filename correct? Is there a vtk for each component, ending in [component number].vtk, starting from 1?" << std::endl;
        exit(error::ERROR_FILE_NOT_FOUND);
    }
}

void Readable_file::check_filetype()
{
    read_extension();

    if (Readable_file::extension_map[m_filetype] != m_extension)
    {
        cerr << "Extension " + m_extension + " doesn't correspond to this filetype!" << endl;
        throw error::ERROR_EXTENSION;
    }
}

void Readable_file::read_extension()
{
    if (m_filename.size() > 0)
    {
        m_extension = m_filename.substr(m_filename.find_last_of(".") + 1);
    }
}

void Lattice_geometry::set_jumps()
{
    switch (dimensions)
    {
    case 1:
        JX = 1;
        JY = 0;
        JZ = 0;
        break;
    case 2:
        JX = (MY);
        JY = 1;
        JZ = 0;
        break;
    case 3:
        JX = (MZ) * (MY);
        JY = (MZ);
        JZ = 1;
        break;
    }
}

void Lattice_geometry::assert_lattice_compatible(Lattice *Lat)
{
    assert(Lat->gradients == (int)dimensions);
    assert(Lat->MX + BOUNDARIES == (int)MX);
    assert(Lat->MY + BOUNDARIES == (int)MY);
    assert(Lat->MZ + BOUNDARIES == (int)MZ);
}

IReader::IReader(Readable_file file)
    : m_file{file.m_filename}
{
}

void IReader::assert_lattice_compatible(Lattice *Lat)
{
    file_lattice.assert_lattice_compatible(Lat);
}

IReader::~IReader() {}

std::vector<std::string> IReader::tokenize(std::string line, char delimiter)
{
    std::istringstream stream{line};

    std::vector<std::string> tokens;
    string token;

    //Read lines as tokens
    while (getline(stream, token, delimiter))
        tokens.emplace_back(token);

    return tokens;
}

void Pro_reader::check_delimiter(const std::string &line)
{

    if (line.find('\t') == string::npos)
    {
        cerr << "Wrong delimiter! Please use tabs.";
        throw ERROR_FILE_FORMAT;
    }
}

void Pro_reader::read_dimensions(const std::vector<std::string> &header_tokens)
{

    if (find(header_tokens.begin(), header_tokens.end(), "z") != header_tokens.end())
        file_lattice.dimensions = 3;

    else if (find(header_tokens.begin(), header_tokens.end(), "y") != header_tokens.end())
        file_lattice.dimensions = 2;

    else if (find(header_tokens.begin(), header_tokens.end(), "x") != header_tokens.end())
        file_lattice.dimensions = 1;

    else
        throw ERROR_FILE_FORMAT;
}

void Pro_reader::check_component_name_format(const std::string &header_token)
{

    std::vector<std::string> component_tokens = tokenize(header_token, ':');

    std::string ERROR = "No headers in the format mol:[molecule]:phi-[monomer].";

    if (component_tokens.size() != NUM_PRO_HEADER_TOKENS
        or component_tokens[0] != "mol"
        or component_tokens[DENSITY_HEADER_TOKEN_INDEX].substr(0, 4) != "phi-")
            throw ERROR;
}

std::vector<std::string> Pro_reader::parse_data(const size_t number_of_components, const size_t first_component_column)
{

    // Prepare vector of vectors for data
    m_data.resize(number_of_components);

    std::vector<std::string> tokens;
    std::string line;

        //Read the actual data, one line at a time
    while (getline(m_file, line)) {
        tokens = tokenize(line, '\t');
        for (size_t i = 0; i < number_of_components; ++i)
            m_data[i].emplace_back(
                //Convert string to float
                strtod( tokens[first_component_column + i].c_str() , NULL )
            );
    }
    
    //Return the last line that has been read.
    return tokens;
}

void Pro_reader::set_lattice_geometry(const std::vector<std::string> &last_line)
{
    switch (file_lattice.dimensions)
    {
    case 3:
        file_lattice.MZ = atof(last_line[Z_DIMENSION].c_str()) + SYSTEM_EDGE_OFFSET;
    case 2:
        file_lattice.MY = atof(last_line[Y_DIMENSION].c_str()) + SYSTEM_EDGE_OFFSET;
    case 1:
        file_lattice.MX = atof(last_line[X_DIMENSION].c_str()) + SYSTEM_EDGE_OFFSET;
        break;
    }

    //We need the correct jump sizes to access the indexed matrix
    file_lattice.set_jumps();
}

void Pro_reader::adjust_indexing()
{
    std::vector<std::vector<Real>> adjusted_data(m_data.size());

    for (vector<Real> &all_components : adjusted_data)
        all_components.resize(m_data[0].size());

    size_t n = 0;
    size_t z = 0;
    do
    {
        size_t y = 0;
        do
        {
            size_t x = 0;
            do
            {
                for (size_t c = 0; c < m_data.size(); ++c)
                    adjusted_data[c][x * file_lattice.JX + y * file_lattice.JY + z * file_lattice.JZ] = m_data[c][n];
                ++n;
                ++x;
            } while (x < file_lattice.MX);
            ++y;
        } while (y < file_lattice.MY);
        ++z;
    } while (z < file_lattice.MZ);

    m_data = adjusted_data;
}

Pro_reader::Pro_reader(Readable_file file)
    : IReader(file), m_data(0)
{
}

std::vector<std::vector<Real>> Pro_reader::get_file_as_vectors()
{
    //Find in which column the density profile starts
    //This depends on the fact that the first mon output is phi
    std::string header_line;

    //Read headers
    getline(m_file, header_line);

    check_delimiter(header_line);

    std::vector<std::string> headers = tokenize(header_line, '\t');
    try
    {

        read_dimensions(headers);

        // First component starts after the dimensions indicators x, y, or z. Remember, first index = 0.
        for (size_t i = file_lattice.dimensions; i < headers.size(); ++i)
            check_component_name_format(headers[i]);
    }
    catch (std::string ERROR)
    {

        cerr << ERROR << endl;
    }

    size_t number_of_components = headers.size() - file_lattice.dimensions;
    size_t first_component_column = IReader::file_lattice.dimensions;

    std::vector<std::string> last_line;

    last_line = parse_data(number_of_components, first_component_column);

    set_lattice_geometry(last_line);

    // Because .pro files are written in x-y-z order, whereas namics uses z-y-x for 3D
    adjust_indexing();

    return m_data;
}

void Vtk_structured_grid_reader::set_lattice_geometry(const std::vector<std::string> &tokens)
{
    switch (file_lattice.dimensions)
    {
    case 3:
        file_lattice.MZ = atof(tokens[OFFSET_VTK_DIMENSIONS_TAG + Z_DIMENSION].c_str()) + BOUNDARIES;
    case 2:
        file_lattice.MY = atof(tokens[OFFSET_VTK_DIMENSIONS_TAG + Y_DIMENSION].c_str()) + BOUNDARIES;
    case 1:
        file_lattice.MX = atof(tokens[OFFSET_VTK_DIMENSIONS_TAG + X_DIMENSION].c_str()) + BOUNDARIES;
        break;
    }

    file_lattice.set_jumps();
}

Vtk_structured_grid_reader::STATUS Vtk_structured_grid_reader::parse_next_data_block(std::vector<Real> &data)
{
    std::string line;

    while (line.find("LOOKUP_TABLE default") == string::npos)
    {
        getline(m_file, line);
    }

    while (getline(m_file, line))
    {
        if (std::regex_match(line, std::regex(R"(^[\d]+.[\d]+$)")))
            return STATUS::NEW_BLOCK_FOUND;
        else
            data.emplace_back(atof(line.c_str()));
    }

    return STATUS::END;
}

std::vector<Real> Vtk_structured_grid_reader::with_bounds(std::vector<Real> &input)
{
    size_t M_bounds = (file_lattice.MX) * (file_lattice.MY) * (file_lattice.MZ);
    vector<Real> output(M_bounds);

    size_t n = 0;
    size_t x = SYSTEM_EDGE_OFFSET;
    do
    {
        size_t y = SYSTEM_EDGE_OFFSET;
        do
        {
            size_t z = SYSTEM_EDGE_OFFSET;
            do
            {
                output[x * file_lattice.JX + y * file_lattice.JY + z * file_lattice.JZ] = input[n];
                ++n;
                ++z;
            } while (z < file_lattice.MZ - SYSTEM_EDGE_OFFSET);
            ++y;
        } while (y < file_lattice.MY - SYSTEM_EDGE_OFFSET);
        ++x;
    } while (x < file_lattice.MX - SYSTEM_EDGE_OFFSET);

    return output;
}

Vtk_structured_grid_reader::Vtk_structured_grid_reader(Readable_file file)
    : IReader(file)
{
}

std::vector<std::vector<Real>> Vtk_structured_grid_reader::get_file_as_vectors()
{
    std::string header_line;

    // Find headers with dimension information
    while (header_line.find("DIMENSIONS") ==std::string::npos)
    {
        getline(m_file, header_line);
    }

    std::vector<std::string> headers;
    headers = tokenize(header_line, ' ');

    if (headers.size() < OFFSET_VTK_DIMENSIONS_TAG+X_DIMENSION+1
        or headers.size() > OFFSET_VTK_DIMENSIONS_TAG+Z_DIMENSION+1)
            throw ERROR_FILE_FORMAT;
    else
        file_lattice.dimensions = headers.size() - OFFSET_VTK_DIMENSIONS_TAG;

    //JX, JY, JZ, needed for the with_bounds function.
    set_lattice_geometry(headers);

    std::vector<std::vector<Real>> output(0);

    std::vector<Real> data;

    Vtk_structured_grid_reader::STATUS status = STATUS::NEW_BLOCK_FOUND;

    while (status == STATUS::NEW_BLOCK_FOUND)
    {
        status = parse_next_data_block(data);
        //ASSUMPTION: VTK files a written without bounds, so add them
        data = with_bounds(data);
        output.emplace_back(data);
        data.clear();
    }

    if (status == STATUS::END)
    {
        data = with_bounds(data);
        output.emplace_back(data);
    } else
    {
        std::cerr << "No blocks found in file" << std::endl;
        throw ERROR_FILE_FORMAT;
    }

    return output;
}

Reader::Reader()
    : m_read_objects(0)
{
    m_file.precision(20);
}

//Callable for multiple files. returns number of objects read.
size_t Reader::read_objects_in(Readable_file file)
{

    switch (file.get_filetype())
    {
    case Readable_filetype::VTK_STRUCTURED_GRID:
        input_reader = make_unique<Vtk_structured_grid_reader>(file);
        break;
    case Readable_filetype::PRO:
        input_reader = make_unique<Pro_reader>(file);
        break;
    default:
        cerr << "Something went horribly wrong while constructing Readable_file. It seems that Readable_filetype is not set properly. Exiting." << endl;
        exit(0);
    }

    std::vector<std::vector<Real>> t_object;
    t_object = input_reader->get_file_as_vectors();

    m_read_objects.insert(m_read_objects.end(), t_object.begin(), t_object.end());

    return t_object.size();
}

void Reader::push_data_to_objects(std::vector<Lattice_object<Real>> &output)
{
    assert(output.size() == m_read_objects.size() && "Please resize your Lattice_object vector before passing!");
    for (size_t i = 0; i < m_read_objects.size(); ++i)
        output[i].m_data = m_read_objects[i];

    m_read_objects.clear();
}

void Reader::assert_lattice_compatible(Lattice *Lat)
{
    if (input_reader)
    {
        input_reader->assert_lattice_compatible(Lat);
    }
    else
    {
        std::cerr << "No input_reader found!" << endl;
    }
}

Reader::~Reader() {}
