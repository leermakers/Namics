#ifndef FILE_READER_H
#define FILE_READER_H

#include <cstdio>
#include <string>
#include <fstream>
#include <vector>
#include "lattice_object.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <memory>
#include <assert.h>
#include <stdlib.h>
#include <map>

typedef double Real;

/*  
 *  These follow a bridge pattern with prototype IReader, concrete classes filetype_reader and bridge Reader.
 *   
 *  WRITING A NEW PARSER? IT'S EASY!
 *      - Add your filetype to the filetype enumerator and extension_map below
 *      - Write your parser and (publicly) inherit and implement IReader
 *      - Add your class to Reader::read_objects_in()
 *      - Good to go!
 */ 

enum class filetype {
            NONE,
            VTK_STRUCTURED_GRID,
            PRO
        };

std::map<filetype, std::string> extension_map {
    {filetype::NONE, ""},
    {filetype::VTK_STRUCTURED_GRID, "vtk"},
    {filetype::PRO, "pro"}
};

//Provides necessary checks before reading file by reader.
class Readable_file {
    public:
        const std::string m_filename;
        
        enum error {
            ERROR_EXTENSION,
            ERROR_FILE_NOT_FOUND,
        };

        filetype get_filetype() {
            return m_filetype;
        }

        //Accepts filename, 
        Readable_file(const std::string filename_, filetype filetype_)
        : m_filename{filename_}, m_filetype{filetype_}
        { 
            try {
                check_filetype();
            } catch (std::string extension) {
                std::cerr << "Extension: " << extension << " not recognized." << std::endl;
                exit(error::ERROR_EXTENSION);
            }

            if ( access( m_filename.c_str(), F_OK ) == -1 ) {
                std::cerr << "Error opening file! Is the filename correct? Is there a vtk for each component, ending in [component number].vtk, starting from 1?" << std::endl;
                exit(error::ERROR_FILE_NOT_FOUND);
            }
        }

    private:

        filetype m_filetype;
        std::string m_extension;

        void check_filetype() {
            
            read_extension();

            if (extension_map[m_filetype] != m_extension) {
                cerr << "Extension " + m_extension + " doesn't correspond to this filetype!" << endl;
                throw error::ERROR_EXTENSION;
            }
        }

        void read_extension() {
            if (m_filename.size() > 0) {
                m_extension = m_filename.substr(m_filename.find_last_of(".")+1);
            }
        }
};


// We need this because we don't know what the lattice in the data file looks like.
// During reading we will infer all these elements and then compare them to the actual lattice
// that was constructed from the .in file by the Lattice class.
// Sometimes the jumps are needed for a particular operation, which is why we provide these as well.

struct Lattice_geometry {
    size_t dimensions{0};
    size_t MX{0};
    size_t MY{0};
    size_t MZ{0};

    size_t JX{0};
    size_t JY{0};
    size_t JZ{0};

    virtual void set_jumps() {
        switch (dimensions) {
            case 1:
            JX=1; JY=0; JZ=0;
            break;
            case 2:
            JX=(MY); JY=1; JZ=0;
            break;
            case 3:
            JX=(MZ)*(MY); JY=(MZ); JZ=1;
            break;
        }
    }

    void assert_lattice_compatible(Lattice* Lat) {
            assert(Lat->gradients == (int)dimensions);
            assert(Lat->MX+2 == (int)MX);
            assert(Lat->MY+2 == (int)MY);
            assert(Lat->MZ+2 == (int)MZ);
    }

};

//Interface for different filetypes
class IReader {
    public:
        //Accepts filename
        IReader(Readable_file file)
        : m_file{file.m_filename}
        { }

        virtual std::vector<std::vector<Real>> get_file_as_vectors() = 0;
        

        void assert_lattice_compatible(Lattice* Lat) {
            file_lattice.assert_lattice_compatible(Lat);
        }

        virtual ~IReader() { }

    protected:
        std::ifstream m_file;
        Lattice_geometry file_lattice;
        enum error {
            ERROR_FILE_FORMAT
        };

        virtual void set_lattice_geometry(std::vector<std::string>&) = 0;

        virtual std::vector<std::string> tokenize(std::string line, char delimiter) {
            std::istringstream stream{line};

            std::vector<std::string> tokens;
            string token;

            //Read lines as tokens
            while (getline(stream, token, delimiter))
                tokens.push_back(token);

            return tokens;
        }

    private:
        //Disable default constructor
        IReader();
};

class Pro_reader : public IReader {
    private:
       std::vector<std::vector<Real>> m_data;

        void check_delimiter(const std::string& line) {

            if (line.find('\t') == string::npos) {
                cerr << "Wrong delimiter! Please use tabs.";
                throw ERROR_FILE_FORMAT;
            }

        }

        void read_dimensions(const std::vector<std::string>& header_tokens) {

            if ( find(header_tokens.begin(), header_tokens.end(), "z") != header_tokens.end() )
                file_lattice.dimensions = 3;
            
            else if ( find(header_tokens.begin(), header_tokens.end(), "y") != header_tokens.end() )
                file_lattice.dimensions = 2;

            else if ( find(header_tokens.begin(), header_tokens.end(), "x") != header_tokens.end() )
                file_lattice.dimensions = 1;

            else
                throw ERROR_FILE_FORMAT;

        }

        void check_component_name_format(const std::string& header_token) {

            std::vector<std::string> component_tokens = tokenize(header_token, ':');

            std::string ERROR = "No headers in the format mol:[molecule]:phi-[monomer].";

            if ( component_tokens.size() != 3
              or component_tokens[0] != "mol"
              or component_tokens[2].substr(0, 4) != "phi-"
            ) throw ERROR;

        }

        std::vector<std::string> parse_data(const size_t number_of_components, const size_t first_component_column) {

            // Prepare vector of vectors for data
            m_data.resize(number_of_components);

            std::vector<std::string> tokens;
            std::string line;
              
            //Read the actual data, one line at a time
            while (getline(m_file, line)) {
                tokens = tokenize(line, '\t');

                for (size_t i = 0; i < number_of_components; ++i)
                    m_data[i].push_back(
                        //Convert string to float
                        strtod( tokens[first_component_column + i].c_str() , NULL )
                    );
            }

            //Return the last line that has been read.
            return tokens;
        }

        void set_lattice_geometry(std::vector<std::string>& last_line) override {
            switch (file_lattice.dimensions) {
    	    	case 3:
    	    		file_lattice.MZ = atof(last_line[2].c_str())+1;
    	    	case 2:
    	    		file_lattice.MY = atof(last_line[1].c_str())+1;
    	    	case 1:
    	    		file_lattice.MX = atof(last_line[0].c_str())+1;
    	    	break;
    	    }

            //We need the correct jump sizes to access the indexed matrix
            file_lattice.set_jumps();
        }

        void adjust_indexing() {
            std::vector< std::vector<Real> > adjusted_data(m_data.size());
	        for (vector<Real>& all_components : adjusted_data)
	        	all_components.resize(m_data[0].size());

	        size_t n = 0;
            size_t z = 0;
	        do { size_t y = 0;
	        	do { size_t x = 0;
	        		do {
	        			for (size_t c = 0 ; c < m_data.size() ; ++c)
	        				adjusted_data[c][x * file_lattice.JX + y * file_lattice.JY + z * file_lattice.JZ] = m_data[c][n];
	        			++n;
	        			++x; } while (x < file_lattice.MX);
	        		++y; } while (y < file_lattice.MY);
	        	++z; } while (z < file_lattice.MZ);
	        m_data = adjusted_data;
        }

    public:
        explicit Pro_reader(Readable_file file)
        : IReader(file), m_data(0)
        { }
        
        std::vector<std::vector<Real>> get_file_as_vectors() override {
            //Find in which column the density profile starts
            //This depends on the fact that the first mon output is phi
            std::string header_line;

            //Read headers
            getline(m_file, header_line);

            check_delimiter(header_line);

            std::vector<std::string> headers = tokenize(header_line, '\t');
            try {

                read_dimensions(headers); 

                // First component starts after the dimensions indicators x, y, or z. Remember, first index = 0.
                for (size_t i = file_lattice.dimensions; i < headers.size(); ++i)
                    check_component_name_format(headers[i]);

            } catch (std::string ERROR) {

                cerr << ERROR << endl;

            }

            size_t number_of_components = headers.size()-file_lattice.dimensions;
            size_t first_component_column = file_lattice.dimensions;

            std::vector<std::string> last_line;

            last_line = parse_data(number_of_components, first_component_column);

            set_lattice_geometry(last_line);

            // Because .pro files are written in x-y-z order, whereas namics uses z-y-x for 3D
            adjust_indexing();

            return m_data;
        }
};


class Vtk_structured_grid_reader : public IReader {
    private:
        void set_lattice_geometry(std::vector<std::string>& tokens) {
            switch (file_lattice.dimensions)
            {
              case 3:
                file_lattice.MZ = atof(tokens[3].c_str())+2;
              case 2:
                file_lattice.MY = atof(tokens[2].c_str())+2;
              case 1:
                file_lattice.MX = atof(tokens[1].c_str())+2;
                break;
            }

            file_lattice.set_jumps();
        }

        void parse_data(std::vector<Real>& data) {
            std::string line;

            while (line.find("LOOKUP_TABLE default") == string::npos ) {
              getline(m_file, line);
            }
            
            while( getline(m_file, line) ) {
              data.push_back(atof(line.c_str()));
            }
        }

        std::vector<Real> with_bounds(std::vector<Real>& input) {
            size_t M_bounds = (file_lattice.MX)*(file_lattice.MY)*(file_lattice.MZ);
	        vector<Real> output(M_bounds);

	        size_t n = 0;
	        size_t x = 1;
	        do {
	        	size_t y = 1;
	        	do {
	        		size_t z = 1;
	        		do {
	        			output[x * file_lattice.JX + y * file_lattice.JY + z * file_lattice.JZ] = input[n];
	        			++n;
	        			++z;
	        		} while (z < file_lattice.MZ-1);
	        		++y;
	        	} while (y < file_lattice.MY-1);
	        	++x;
	        } while (x < file_lattice.MX-1);

	        return output;
        }

    public:
        Vtk_structured_grid_reader(Readable_file file)
        : IReader(file)
        { }

        std::vector< std::vector<Real> > get_file_as_vectors() {
            
            std::string header_line;

            // Find headers with dimension information
            while (header_line.find("DIMENSIONS") == string::npos ) {
                getline(m_file, header_line);
            }
            
            std::vector<std::string> headers;
            headers = tokenize(header_line, ' ');

            if (headers.size() < 2 or headers.size() > 4)
                throw ERROR_FILE_FORMAT;
            else
                file_lattice.dimensions = headers.size()-1;

            //JX, JY, JZ, needed for the with_bounds function.
            set_lattice_geometry(headers);

            std::vector<Real> data;
            parse_data(data);

            //ASSUMPTION: VTK files a written without bounds, so add them
            data = with_bounds(data);

            std::vector< std::vector<Real> > output(0);
            output.push_back(data);
            
            return output;
        }
};

class Reader {
    public:
        Reader()
        : m_read_objects(0)
        {  }

        //Callable for multiple files. returns number of objects read.
        size_t read_objects_in(Readable_file file) {

            switch ( file.get_filetype() ) {
                case filetype::VTK_STRUCTURED_GRID:
                    input_reader = make_unique<Vtk_structured_grid_reader>(file);
                    break;
                case filetype::PRO:
                    input_reader = make_unique<Pro_reader>(file);
                    break;
                default:
                    cerr << "Something went horribly wrong while constructing Readable_file. It seems that filetype is not set properly. Exiting." << endl;
                    exit(0);
            }

            std::vector< std::vector<Real> > t_object;
            t_object = input_reader->get_file_as_vectors();

            m_read_objects.insert( m_read_objects.end(), t_object.begin(), t_object.end() );

            return t_object.size();
        }

        void push_data_to_objects(std::vector< Lattice_object<Real> >& output) {
            assert(output.size() == m_read_objects.size() && "Please resize your Lattice_object vector before passing!");
            for (size_t i = 0 ; i < m_read_objects.size() ; ++i)
                output[i].m_data = m_read_objects[i];
                
            m_read_objects.clear();
        }

        void assert_lattice_compatible(Lattice* Lat) {
            if (input_reader) {
                input_reader->assert_lattice_compatible(Lat);
            } else {
                std::cerr << "No input_reader found!" << endl;
            }
        }

        ~Reader() {
        }

    private:
        std::vector< std::vector<Real> > m_read_objects;
        std::ifstream m_file;
        unique_ptr<IReader> input_reader;
};


#endif