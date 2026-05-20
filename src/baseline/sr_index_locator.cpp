#include "sr_index_locator.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <sr-index/config.h>
#include <sr-index/construct.h>
#include <sr-index/sr_index.h>

namespace lz77tax {

// ─── Implementación oculta (PIMPL) ───────────────────────────────────────────

struct SrIndexLocator::Impl {
    using SrIdx = sri::SrIndexValidArea<
        std::reference_wrapper<sri::GenericStorage>,
        sri::Alphabet<8>>;

    sri::GenericStorage          storage;
    std::unique_ptr<SrIdx>       idx;

    std::vector<std::size_t> locate(const std::string& pattern) const {
        using Sequence = sri::Alphabet<8>::string_type;
        Sequence seq(pattern.begin(), pattern.end());
        return idx->Locate(seq);
    }
};

// ─── Constructor / destructor / move ─────────────────────────────────────────

SrIndexLocator::SrIndexLocator()  = default;
SrIndexLocator::~SrIndexLocator() = default;

SrIndexLocator::SrIndexLocator(SrIndexLocator&&) noexcept            = default;
SrIndexLocator& SrIndexLocator::operator=(SrIndexLocator&&) noexcept = default;

// ─── load / build ─────────────────────────────────────────────────────────────

void SrIndexLocator::load(const std::string& data_name,
                           const std::string& data_dir,
                           std::size_t sr) {
    impl_ = std::make_unique<Impl>();
    impl_->idx = std::make_unique<Impl::SrIdx>(std::ref(impl_->storage), sr);

    sri::Config config(data_name, data_dir,
                       sri::SDSL_LIBDIVSUFSORT,
                       /*load=*/true, /*width=*/8);
    impl_->idx->load(config);
}

void SrIndexLocator::build(const std::string& data_path,
                            const std::string& cache_dir,
                            std::size_t sr) {
    impl_ = std::make_unique<Impl>();
    impl_->idx = std::make_unique<Impl::SrIdx>(std::ref(impl_->storage), sr);

    sri::Config config(data_path, cache_dir,
                       sri::SDSL_LIBDIVSUFSORT,
                       /*load=*/false, /*width=*/8);
    sri::construct(*impl_->idx, data_path, config);
}

// ─── Queries ──────────────────────────────────────────────────────────────────

std::pair<std::size_t, std::size_t>
SrIndexLocator::locate_extremal(const std::string& pattern) const {
    const auto positions = impl_->locate(pattern);
    if (positions.empty()) return {SIZE_MAX, 0};
    const auto [it_min, it_max] =
        std::minmax_element(positions.begin(), positions.end());
    return {*it_min, *it_max};
}

std::size_t SrIndexLocator::locate_count(const std::string& pattern) const {
    return impl_->locate(pattern).size();
}

std::size_t SrIndexLocator::size_in_bytes() const {
    std::ostringstream ss;
    impl_->idx->serialize(ss, nullptr, "");
    return static_cast<std::size_t>(ss.tellp());
}

bool SrIndexLocator::is_loaded() const {
    return impl_ != nullptr && impl_->idx != nullptr;
}

}  // namespace lz77tax
