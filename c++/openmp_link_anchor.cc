#include "openmp.hpp"

namespace {

#if NRG_MKL_REQUIRES_OPENMP
struct OpenMPRuntimeLinkAnchor {
  OpenMPRuntimeLinkAnchor() { NRG::detail::reference_linked_openmp_runtime(); }
};

[[maybe_unused]] const OpenMPRuntimeLinkAnchor openmp_runtime_link_anchor;
#endif

} // namespace
