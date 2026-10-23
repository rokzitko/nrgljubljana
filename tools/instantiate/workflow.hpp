#ifndef NRG_INSTANTIATE_WORKFLOW_HPP
#define NRG_INSTANTIATE_WORKFLOW_HPP

#include <string>

namespace NRG::Workflow::detail {
// Called under the workflow preparation guard, with the calculation CWD active.
double instantiate(const std::string &parameter_file, const std::string &template_directory, unsigned nmax);
}

#endif
