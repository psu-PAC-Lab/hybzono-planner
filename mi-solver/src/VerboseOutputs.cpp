#include "VerboseOutputs.hpp"

using namespace MI_QPIPMPC;

// copy assignment
VerboseOutputs& VerboseOutputs::operator=(const VerboseOutputs &other)
{
    // lock the mutex
    std::lock_guard<std::mutex> lock(mtx);

    // copy the output stream
    os.str(other.os.str());

    // copy the function pointer
    output_stream_fcn = other.output_stream_fcn;

    return *this;
}

// print methods - portions copied from mi_osqp
void VerboseOutputs::print_headline()
{
    // lock mutex
    std::lock_guard<std::mutex> lock(mtx);

    // update the output stream
    os << std::endl << std::setw(15) << "Iter" << std::setw(15) << "Queue" << std::setw(15) <<
        "J_min" << std::setw(15) << "J_max" << std::setw(15) << "Gap [%]" << std::endl;

    // print if print function defined
    if (output_stream_fcn != nullptr)
    {
        // send to print function pointer
        output_stream_fcn(os.str());

        // reset output stream
        os.str("");
    }
}

void VerboseOutputs::print_progress(const std::unique_ptr<Node> &leaf, const LeavesQueue &leaves, const DataVolatile &data_volatile)
{
    // lock the mutex
    std::lock_guard<std::mutex> lock(mtx);

    // get data
    int iter_num = data_volatile.get_iter_num();
    int num_leaves = leaves.size();
    double lower_glob = data_volatile.get_lower_glob();
    double upper_glob = data_volatile.get_upper_glob();

    // update the output stream
    os << std::setw(15) << iter_num << std::setw(15) << num_leaves << std::setw(15) <<
        std::scientific << std::setprecision(3) << lower_glob << std::setw(15) << upper_glob <<
        std::setw(15) << ((upper_glob-lower_glob)/std::abs(lower_glob))*100 << std::endl;

    // print if print function defined
    if (output_stream_fcn != nullptr)
    {
        // send to print function pointer
        output_stream_fcn(os.str());

        // reset output stream
        os.str("");
    }
}

void VerboseOutputs::print_footer(const DataVolatile &data_volatile)
{
    // lock the mutex
    std::lock_guard<std::mutex> lock(mtx);

    // get data
    MIQP_Status status = data_volatile.get_status();
    double upper_glob = data_volatile.get_upper_glob();
    int qp_iter = data_volatile.get_qp_iter();
    double run_time = data_volatile.get_run_time();

    // send to stream
    os << std::endl;
    os << "Status: " << status << std::endl;
    switch (status)
    {
        case QP_SOLVED:
        {
            os << "Optimal solution found. Objective: " << upper_glob << std::endl;
            break;
        }
        case MIQP_INFEASIBLE:
        {
            os << "Solution is infeasible." << std::endl;
            break;
        }
        case MIQP_MAX_ITER_FEASIBLE:
        {
            os << "Maximum number of iterations reached or max time elapsed. Solution is feasible." << std::endl;
            os << "Objective: " << upper_glob << std::endl;
            break;
        }
        case MIQP_MAX_ITER_INFEASIBLE:
        {
            os << "Maximum number of iterations reached or max time elapsed. Solution is infeasible." << std::endl;
            break;
        }
    }
    os << "Total number of QP iterations: " << qp_iter << std::endl;
    os << "Total solve time: " << run_time << " seconds" << std::endl;

    // print if print function defined
    if (output_stream_fcn != nullptr)
    {
        // send to print function pointer
        output_stream_fcn(os.str());

        // reset output stream
        os.str("");
    }
}

void VerboseOutputs::set_output_stream_fcn(void (*fcn)(const std::string &str))
{
    // lock the mutex
    std::lock_guard<std::mutex> lock(mtx);

    // set the function pointer
    output_stream_fcn = fcn;
}

std::string VerboseOutputs::get_output_string() const
{
    std::lock_guard<std::mutex> lock(mtx);
    return os.str();
}

void VerboseOutputs::write_to_output_stream(const std::string &str)
{
    // lock the mutex
    std::lock_guard<std::mutex> lock(mtx);

    // update the output stream
    os << str;

    // print if print function defined
    if (output_stream_fcn != nullptr)
    {
        // send to print function pointer
        output_stream_fcn(os.str());

        // reset output stream
        os.str("");
    }
}

void VerboseOutputs::reset()
{
    // lock the mutex
    std::lock_guard<std::mutex> lock(mtx);

    // reset the output stream
    os.str("");
}