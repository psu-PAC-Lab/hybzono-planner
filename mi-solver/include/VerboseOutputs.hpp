#ifndef _MI_QPIPMPC_VERBOSEOUTPUTS_HPP_
#define _MI_QPIPMPC_VERBOSEOUTPUTS_HPP_

#include <sstream>
#include <mutex>
#include <string>
#include <memory>
#include <iomanip>
#include "Node.hpp"
#include "Data.hpp"
#include "LeavesQueue.hpp"
#include "DataVolatile.hpp"
#include "QP_consts.hpp"

namespace MI_QPIPMPC
{

class VerboseOutputs
{
    public:
        VerboseOutputs() = default;
        VerboseOutputs& operator=(const VerboseOutputs &other);

        void set_output_stream_fcn(void (*fcn)(const std::string &str));
        void print_headline();
        void print_progress(const std::unique_ptr<Node> &leaf, const LeavesQueue& leaves, const DataVolatile& data_volatile);
        void print_footer(const DataVolatile& data_volatile);
        std::string get_output_string() const;

        void write_to_output_stream(const std::string &str);

        void reset();

    protected:
        std::ostringstream os;
        void (*output_stream_fcn)(const std::string &str) = nullptr;
    
    private:
        mutable std::mutex mtx;

};

} // end namespace MI_QPIPMPC

#endif