#include "file_writer.h"
#include "../lattice.h"

// Writing a new writer? Check out the header for instructions

Register_class<IProfile_writer, Vtk_structured_grid_writer, Writable_filetype, Lattice*, Writable_file> Vtk_structured_grid_writer_factory(Writable_filetype::VTK_STRUCTURED_GRID);
Register_class<IProfile_writer, Vtk_structured_points_writer, Writable_filetype, Lattice*, Writable_file> Vtk_structured_points_writer_factory(Writable_filetype::VTK_STRUCTURED_POINTS);
Register_class<IProfile_writer, Pro_writer, Writable_filetype, Lattice*, Writable_file> Pro_writer_factory(Writable_filetype::PRO);

map<std::string, Writable_filetype> Profile_writer::output_options {
        {"vtk", Writable_filetype::VTK_STRUCTURED_POINTS},
        {"vtk_structured_grid", Writable_filetype::VTK_STRUCTURED_GRID},
        {"vtk_structured_points", Writable_filetype::VTK_STRUCTURED_POINTS},
        {"pro", Writable_filetype::PRO},
    };

Register_class<IParameter_writer, Kal_writer, Writable_filetype, Writable_file> Kal_writer_factory(Writable_filetype::KAL);
Register_class<IParameter_writer, Csv_parameter_writer, Writable_filetype, Writable_file> Csv_parameter_writer_factory(Writable_filetype::CSV);

map<std::string, Writable_filetype> Parameter_writer::output_options {
        {"kal", Writable_filetype::KAL},
        {"csv", Writable_filetype::CSV},
    };

typedef std::string Extension;

std::map<Writable_filetype, Extension> Writable_file::extension_map {
    {Writable_filetype::KAL, "kal"},
    {Writable_filetype::VTK_STRUCTURED_GRID, "vtk"},
    {Writable_filetype::VTK_STRUCTURED_POINTS, "vtk"},
    {Writable_filetype::CSV, "csv"},
    {Writable_filetype::PRO, "pro"}
};

Writable_file::Writable_file(const std::string filename_, Writable_filetype filetype_, int identifier_)
: m_identifier{identifier_}, m_filename{filename_}, m_filetype{filetype_}
{
    append_extension();
}

Writable_file::~Writable_file()
{

}

void Writable_file::increment_identifier()
{
    string dot = "\\.";
    string s_regex = "(\\_" + to_string(m_identifier) + dot + Writable_file::extension_map[m_filetype] + ')';
    std::regex ex(s_regex);

    string replacement = "_"+ to_string(++m_identifier) + "." + Writable_file::extension_map[m_filetype];

    m_filename = regex_replace(m_filename,ex, replacement);
}

Writable_filetype Writable_file::get_filetype()
{
    return m_filetype;
}

string Writable_file::get_filename()
{
    return m_filename;
}

void Writable_file::append_extension()
{
    m_filename.append("_" + to_string(m_identifier) + "." + Writable_file::extension_map[m_filetype]);
}

IProfile_writer::IProfile_writer(Lattice* geometry_, Writable_file file_)
: m_geometry{geometry_}, m_file{file_}, m_adapter{geometry_}
{
    m_filestream.precision(configuration.precision);
}

IProfile_writer::~IProfile_writer()
{

}

void IProfile_writer::bind_data(map< string, shared_ptr<IOutput_ptr>>& profiles_)
{
    m_profiles = profiles_;
}

void IProfile_writer::bind_subystem_loop(Boundary_mode mode_) {
    switch (mode_) {
        case Boundary_mode::WITH_BOUNDS:
            subsystem_loop = std::bind(&Lattice_accessor::system_plus_bounds, m_adapter, std::placeholders::_1);
            break;
        case Boundary_mode::WITHOUT_BOUNDS:
            subsystem_loop = std::bind(&Lattice_accessor::skip_bounds, m_adapter, std::placeholders::_1);
            break;
    }
}


Vtk_structured_grid_writer::Vtk_structured_grid_writer(Lattice* geometry_, Writable_file file_)
: IProfile_writer(geometry_, file_)
{
    configuration.boundary_mode = Boundary_mode::WITHOUT_BOUNDS;
}

Vtk_structured_grid_writer::~Vtk_structured_grid_writer()
{

}

void Vtk_structured_grid_writer::prepare_for_data()
{
    bind_subystem_loop(configuration.boundary_mode);

    m_filestream.open(m_file.get_filename(), std::ios_base::out);

	std::ostringstream vtk;

    int MX = m_geometry->MX;
	int MY = m_geometry->MY;
	int MZ = m_geometry->MZ;

	vtk << "# vtk DataFile Version 4.2 \n";
	vtk << "VTK output \n";
	vtk << "ASCII\n";
	vtk << "DATASET STRUCTURED_GRID \n";
	vtk << "DIMENSIONS " << MX << " " << MY << " " << MZ << "\n";
	vtk << "POINTS " << MX * MY * MZ << " int\n";

	for (int x = 1; x < MX + 1; ++x)
		for (int y = 1; y < MY + 1; ++y)
			for (int z = 1 ; z < MZ + 1 ; ++z )
				vtk << x << " " << y << " " << z << "\n";

	vtk << "POINT_DATA " << MX * MY * MZ << "\n";

	m_filestream << vtk.str();
	m_filestream.flush();

	m_filestream.close();
}

void Vtk_structured_grid_writer::write()
{
	m_filestream.open(m_file.get_filename(), std::ios_base::app);

    for (auto& profile : m_profiles) {

	    std::ostringstream vtk;

	    vtk << "SCALARS " << profile.first << " float\nLOOKUP_TABLE default \n";

	    subsystem_loop(
            [this, &vtk, profile] (size_t x, size_t y, size_t z)
            {
	    		vtk << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
            }
        );

	    m_filestream << vtk.str();
	    m_filestream.flush();

    }

	m_filestream.close();
    m_file.increment_identifier();
}

Vtk_structured_points_writer::Vtk_structured_points_writer(Lattice* geometry_, Writable_file file_)
: IProfile_writer(geometry_, file_)
{
    configuration.boundary_mode = Boundary_mode::WITHOUT_BOUNDS;
}

Vtk_structured_points_writer::~Vtk_structured_points_writer()
{

}

void Vtk_structured_points_writer::prepare_for_data()
{
    bind_subystem_loop(configuration.boundary_mode);

    m_filestream.open(m_file.get_filename(), std::ios_base::out);

	std::ostringstream vtk;

    int MX = m_geometry->MX;
	int MY = m_geometry->MY;
	int MZ = m_geometry->MZ;

    vtk << "# vtk DataFile Version 4.2 \n";
	vtk << "VTK output \n";
	vtk << "ASCII\n";
	vtk << "DATASET STRUCTURED_POINTS \n";
    vtk << "SPACING 1 1 1 \n";
    vtk << "ORIGIN 0 0 0";
    vtk << "POINT_DATA " << MX*MY*MZ << "\n";
	vtk << "DIMENSIONS " << MX << " " << MY << " " << MZ << "\n";

    m_filestream << vtk.str();
	m_filestream.flush();

	m_filestream.close();
}

void Vtk_structured_points_writer::write()
{
    m_filestream.open(m_file.get_filename(), std::ios_base::app);

    for (auto& profile : m_profiles) {

	    std::ostringstream vtk;

        vtk << "SCALARS " << profile.first << " float\nLOOKUP_TABLE default \n";

        #ifdef PAR_MESODYN
            stl::host_vector<Real> m_data;

            subsystem_loop(
                [this, &vtk, profile] (size_t x, size_t y, size_t z) {
	    		    vtk << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
                }
            );
        #else
        subsystem_loop(
            [this, &vtk, profile] (size_t x, size_t y, size_t z) {
	    		vtk << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
            }
        );
        #endif

        
        m_filestream << vtk.str();
	    m_filestream.flush();
    }

	m_filestream.close();
    m_file.increment_identifier();
}

Pro_writer::Pro_writer(Lattice* geometry_, Writable_file file_)
: IProfile_writer(geometry_, file_)
{
    configuration.boundary_mode = Boundary_mode::WITH_BOUNDS;
}

Pro_writer::~Pro_writer()
{

}

void Pro_writer::prepare_for_data()
{
    bind_subystem_loop(configuration.boundary_mode);
    
    m_filestream.open(m_file.get_filename(), std::ios_base::out);

	std::ostringstream pro;

    pro << "x";

    if (m_geometry->gradients > 1)
    {
        pro << "\ty";
    }   

    if (m_geometry->gradients > 2)
    {
        pro << "\tz";
    }

    pro << '\n';

    subsystem_loop(
        [this, &pro] (size_t x, size_t y, size_t z) {
			pro << x;
            if (m_geometry->gradients > 1)
            {
                pro << "\t" << y;
                if (m_geometry->gradients > 1)
                {
                    pro << "\t" << z;
                }

            }
            pro << '\n';
        }
    );

    m_filestream << pro.str();
	m_filestream.flush();

	m_filestream.close();
}

void Pro_writer::write()
{
    string temp_filename = m_file.get_filename()+".temp";

    for (auto& profile : m_profiles) {
        m_filestream.open(temp_filename, std::ios_base::out);
        std::ifstream in_file(m_file.get_filename());
        if (!in_file) {
            std::cerr << "Could not open input file. Did you prepare for data?\n";
            throw 1;
        }

	    std::ostringstream pro;

        std::string in;
        std::getline(in_file, in);
        m_filestream << in << "\t" << profile.first << '\n';

        subsystem_loop(
            [this, &pro, &in, &in_file, profile] (size_t x, size_t y, size_t z) {
                std::getline(in_file, in);
                m_filestream << in << "\t" << profile.second->data(x*m_geometry->JX+y*m_geometry->JY+z*m_geometry->JZ) << "\n";
                m_filestream << pro.str();
	            m_filestream.flush();
            }
        );

        in_file.close();
        m_filestream.close();

        rename(temp_filename.c_str(), m_file.get_filename().c_str());
    }

	m_filestream.close();
    m_file.increment_identifier();
}

IParameter_writer::IParameter_writer(Writable_file file_)
: m_file{file_}
{
    m_filestream.precision(configuration.precision);
}

void IParameter_writer::bind_data(map< string, shared_ptr<IOutput_ptr>>& param_) {
    m_params = param_;
}

Kal_writer::Kal_writer(Writable_file file_)
: IParameter_writer(file_)
{
    
}

Kal_writer::~Kal_writer() {

}

void Kal_writer::write() {
    m_filestream.open(m_file.get_filename(), ios_base::app);
    for (auto& description : m_selected_variables) {
        m_filestream << m_params[description]->data() << '\t';
    }
    m_filestream << '\n';
    m_filestream.close();
}

void Kal_writer::prepare_for_data(vector<string>& selected_variables_) {
    m_selected_variables = selected_variables_;

    m_filestream.open(m_file.get_filename());

    for (auto& description : m_selected_variables)
        m_filestream << description << '\t';
    m_filestream << '\n';

    m_filestream.close();
}

Csv_parameter_writer::Csv_parameter_writer(Writable_file file_)
: IParameter_writer(file_)
{
    
}

Csv_parameter_writer::~Csv_parameter_writer() {

}

void Csv_parameter_writer::write() {
    m_filestream.open(m_file.get_filename(), ios_base::app);
    for (auto& description : m_selected_variables) {
        m_filestream << m_params[description]->data() << ',';
    }
    m_filestream << '\n';
    m_filestream.close();
}

void Csv_parameter_writer::prepare_for_data(vector<string>& selected_variables_) {
    m_selected_variables = selected_variables_;

    m_filestream.open(m_file.get_filename());

    for (auto& description : m_selected_variables)
        m_filestream << description << ',';
    m_filestream << '\n';

    m_filestream.close();
}