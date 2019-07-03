#ifndef FILE_WRITER_H
#define FILE_WRITER_H

//  ***** Writing new File writers *****
//
//  If you follow the method below, everything will automatically work with any existing code
//
//  Inherit IProfile_writer when writing arrays of size M (system size incl bounds)
//  Inherit IParameter_writer when writing invidual parameters, timesteps are supported.
//
//  You can basically overwrite two functions of the interface and you'll be done:
//      -   void prepare_for_data()
//          
//          Meant to do operations that have to be done before any time-stepping is being done.
//          In most writers use this to write headers or a pre-defined coordinate system.
//          This is called only once when the output class initializes.
//
//      -   void write();
//
//          Used to write the actual data to file, everything that cannot be done on a per-timestep basis
//          has to be done in prepare_for_data(). This is called for any timestep in a program.
//
//      -   NOTE: Files can then be opened using: m_filestream.open(m_file.get_filename());
//      -   NOTE: If you want to categorise your parameter data, a Category_map is available from IParameter_writer
//
//  When you're done, we need to make sure your new writer is registered with the factory classes so that
//  it is recognized when the user asks for your new writer:
//
//      -   Add an entry that describes your writer in the Writable_filetype enum below (use capitals)
//      -   In file_writer.cpp register your class with the factory:
//          Register_class< Interface, Implementation, Writable_filetype, Lattice* (only for Profile writers), Writable_file> Pro_writer_factory(Writable_filetype::YOURTYPE);
//          where Interface is either IProfile_writer or IParameter_writer and Implementation is your class that inherited either. Remove Lattice* if you've written a
//          parameter writer. Replace the function argument at the end by your newly added Writable_filetype from the previous step.
//      -   Add a key for the input file for your writer by adding your Writable_filetype and the corresponding key to output_options of
//          the corresponding namespace in file_writer.cpp.
//      -   Feed the writers the correct extension by adding your Writable_filetype to the extension_map in file_writer.cpp
//         
//  DONE! Everything should now automatically work with all the existing code.
//

#include "factory.h"
#include "../tools.h"
#include "stl_typedef.h"
#include "lattice_accessor.h"
#include <string>
#include <cstdio>
#include <vector>
#include <fstream>
#include <sstream>
#include <functional>
#include <iostream>
#include <regex>

#define Output_as_metadata IParameter_writer::CATEGORY::METADATA
#define Output_as_timespan IParameter_writer::CATEGORY::TIMESPAN
#define Output_as_constant IParameter_writer::CATEGORY::CONSTANT  

class Lattice;

enum class Writable_filetype
{
    KAL,
    CSV,
    VTK_STRUCTURED_GRID,
    VTK_STRUCTURED_POINTS,
    PRO
};

struct Header {
    std::string OUT_key;
	std::string OUT_name;
	std::string OUT_prop;
};

class IOutput_ptr {
    public:
        virtual string data(const size_t = 0) = 0;
};

template<typename T>
class Output_ptr : public IOutput_ptr
{
    public:
        Output_ptr(T* parameter_)
        : parameter{parameter_}, buffer{0}
        {
        }

        void set_buffer(const size_t size) {
            buffer.resize(size);
            #ifdef PAR_MESODYN
            TransferDataToHost(buffer.data(), const_cast<T*>(parameter), const_cast<size_t&>(size));
            #else
            for (size_t i = 0 ; i < size ; ++i)
                buffer[i] = parameter[i];
            #endif
        }

        void clear_buffer() {
            buffer.clear();
        }

        string data(const size_t offset = 0) override
        {
            if (buffer.size() > 0) {
                ostringstream out;
                out << buffer[offset];
                return out.str();
            } else {

                const T* data = new T;

            #ifdef PAR_MESODYN
                TransferDataToHost(const_cast<T*>(data), const_cast<T*>(parameter+offset), 1);
            #else
                data = parameter+offset;
            #endif

                ostringstream out;
                out << *data;
                return out.str();
            }
        }

        const T* parameter;

    private:
        std::vector<Real> buffer;
};

class Writable_file {
    public:
        Writable_file(const std::string, Writable_filetype, int = 0);
        ~Writable_file();
        Writable_filetype get_filetype();
        string get_filename();
        void increment_identifier();
        
    protected:
        int m_identifier;
        std::string m_filename;
        static std::map<Writable_filetype, std::string> extension_map;
        const Writable_filetype m_filetype;
        void append_extension();
};

class IProfile_writer
{
    public:
        IProfile_writer(Lattice*, Writable_file);
        ~IProfile_writer();

        virtual void write() = 0;

        enum class Boundary_mode
        {
            WITH_BOUNDS,
            WITHOUT_BOUNDS
        };

        struct Configuration {
            Boundary_mode boundary_mode;
            size_t precision = DEFAULT_PRECISION;
        } configuration;

        virtual void prepare_for_data() = 0;
        virtual void bind_data(map< string, shared_ptr<IOutput_ptr>>&);

        static constexpr uint8_t DEFAULT_PRECISION = 14;

    protected:
        Lattice* m_geometry;
        Writable_file m_file;
        std::map< string, shared_ptr<IOutput_ptr> > m_profiles;
        std::ofstream m_filestream;
        Lattice_accessor m_adapter;

        virtual void bind_subystem_loop(Boundary_mode);

        std::function<void(
            std::function<void(size_t, size_t, size_t)>
        )> subsystem_loop;
};

namespace Profile_writer {
    typedef Factory_template<IProfile_writer, Writable_filetype, Lattice*, Writable_file> Factory;
    extern map<std::string, Writable_filetype> output_options;
}

class Vtk_structured_grid_writer : public IProfile_writer
{
    public:
        Vtk_structured_grid_writer(Lattice*, Writable_file);
        ~Vtk_structured_grid_writer();

        virtual void write() override;
        virtual void prepare_for_data() override;
};

class Vtk_structured_points_writer : public IProfile_writer
{
    public:
        Vtk_structured_points_writer(Lattice*, Writable_file);
        ~Vtk_structured_points_writer();

        virtual void write() override;
        virtual void prepare_for_data() override;

};
class Pro_writer : public IProfile_writer
{
    public:
        Pro_writer(Lattice*, Writable_file);
        ~Pro_writer();

        virtual void write() override;
        virtual void prepare_for_data() override;
};


class IParameter_writer
{
    static constexpr uint8_t DEFAULT_PRECISION = 14;

    public:
        IParameter_writer(Writable_file);

        typedef std::map<string, shared_ptr<IOutput_ptr>> Prop_map;

        struct Configuration {
            size_t precision = DEFAULT_PRECISION;
            ios_base::open_mode write_mode = ios_base::app;
        } configuration;

        virtual void write() = 0;
        virtual void prepare_for_data(vector<string>&) = 0;
        virtual void bind_data(Prop_map&);
        

    protected:
        Prop_map m_params;
        vector<string> m_selected_variables;
        Writable_file m_file;
        std::ofstream m_filestream;
};

namespace Parameter_writer {
    typedef Factory_template<IParameter_writer, Writable_filetype, Writable_file> Factory;
    extern map<std::string, Writable_filetype> output_options;
}

class JSON_parameter_writer : public IParameter_writer
{
    public:
        JSON_parameter_writer(Writable_file);
        ~JSON_parameter_writer();


        void write();
        void prepare_for_data(vector<string>&);

        enum class STATE {
            NONE,
            IS_OPEN,
            IS_CLOSED,
        };

        STATE get_state();
        void finalize();
        
    private:
        STATE m_state;
        std::vector<std::string> m_constants;
        std::vector<std::string> m_selected_constants;
        std::vector<std::string> m_metadata;
        std::vector<std::string> m_selected_metadata;
        std::vector<std::string> m_timespan;
        std::vector<std::string> m_selected_timespan;
        void partition_selected_variables(std::vector<std::string>&, std::vector<std::string>&);
        void preprocess_categories();
        void write_list(std::vector<std::string>&);
        bool is_number(std::string&);
        void write_array_object(std::vector<std::string>&);
        void open_json();
        void close_json();
        void open_object(const std::string&);
        void close_object();
        void open_array(const std::string&);
        void close_array();
};

class Kal_writer : public IParameter_writer
{
    public:
        Kal_writer(Writable_file);
        ~Kal_writer();

        void write();
        void prepare_for_data(vector<string>&);
};

class Csv_parameter_writer : public IParameter_writer
{
    public:
        Csv_parameter_writer(Writable_file);
        ~Csv_parameter_writer();

        void write();
        void prepare_for_data(vector<string>&);
};

#endif