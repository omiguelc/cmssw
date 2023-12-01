#include "L1Trigger/L1TMuonEndCapPhase2/interface/Tensors/NdArrayDesc.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DebugUtils.h"

using namespace emtf::phase2;

NdArrayDesc::NdArrayDesc(std::initializer_list<int>&& extents):
    num_dimensions_(extents.size()),
    num_elements_(0),
    extents_(new unsigned int[num_dimensions_]),
    strides_(new unsigned int[num_dimensions_])
{
    auto extents_it = extents.begin();

    const int last_index = num_dimensions_ - 1;

    for (int i = last_index; i > -1; --i) {
        extents_[i] = extents_it[i];

        if (i == last_index)
            strides_[i] = 1;
        else
            strides_[i] = extents_[i + 1] * strides_[i + 1];
    }

    num_elements_ = extents_[0] * strides_[0];
}

NdArrayDesc::~NdArrayDesc() {
    // Do Nothing
}

const unsigned int& NdArrayDesc::num_dimensions() const {
    return num_dimensions_;
}

const unsigned int& NdArrayDesc::num_elements() const {
    return num_elements_;
}

unsigned int NdArrayDesc::get_index_next_dim_contrib(const unsigned int& i_dim) const {
    emtf_assert(i_dim == num_dimensions_);  // Check that all the required indices were given

    return 0;
}
