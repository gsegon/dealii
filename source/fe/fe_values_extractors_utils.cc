// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/fe/fe_values_extractors_utils.h>
#include <deal.II/fe/fe_values_views.h>


DEAL_II_NAMESPACE_OPEN

namespace FEValuesExtractors
{
  namespace utils
  {
    namespace internal
    {
      template <int dim, int spacedim>
      std::vector<unsigned int>
      generate_component_order(
        const std::vector<FEValuesExtractors::AnyExtractor> &extractors,
        const unsigned int                                   fe_n_components)
      {
        // Initialize `component_order` with invalid unsigned ints representing
        // unassigned state.
        std::vector<unsigned int> component_order(
          fe_n_components, numbers::invalid_unsigned_int);

        // Extract the start component index and determine the number of
        // components for each extractor.
        unsigned int block_index = 0;
        for (const auto &extractor : extractors)
          {
            auto start_component_index = numbers::invalid_unsigned_int;
            auto n_components          = numbers::invalid_unsigned_int;

            if (std::holds_alternative<FEValuesExtractors::Scalar>(extractor))
              {
                start_component_index =
                  std::get<FEValuesExtractors::Scalar>(extractor).component;
                n_components = 1;
              }
            else if (std::holds_alternative<FEValuesExtractors::Vector>(
                       extractor))
              {
                start_component_index =
                  std::get<FEValuesExtractors::Vector>(extractor)
                    .first_vector_component;
                n_components = FEValuesViews::Vector<dim, spacedim>::
                  value_type::n_independent_components;
              }
            else if (std::holds_alternative<FEValuesExtractors::Tensor<2>>(
                       extractor))
              {
                start_component_index =
                  std::get<FEValuesExtractors::Tensor<2>>(extractor)
                    .first_tensor_component;
                n_components = FEValuesViews::Tensor<2, dim, spacedim>::
                  value_type::n_independent_components;
              }
            else if (std::holds_alternative<
                       FEValuesExtractors::SymmetricTensor<2>>(extractor))
              {
                start_component_index =
                  std::get<FEValuesExtractors::SymmetricTensor<2>>(extractor)
                    .first_tensor_component;
                n_components =
                  FEValuesViews::SymmetricTensor<2, dim, spacedim>::value_type::
                    n_independent_components;
              }
            else
              {
                AssertThrow(
                  false,
                  ExcNotImplemented(
                    "An unsupported ExtractorVariant was passed in the component_wise extractor_order argument."));
              }

            // Fill `component_order` vector with `n_components` starting at
            // `start_component_index`. Set the values to `block_index`.
            for (unsigned int i = start_component_index;
                 i < start_component_index + n_components;
                 ++i)
              {
                AssertThrow(i < component_order.size(),
                            ExcIndexRange(i, 0, component_order.size()));

                AssertThrow(
                  component_order[i] == numbers::invalid_unsigned_int,
                  ExcMessage(
                    "A component which has already been assigned a block "
                    "index is trying to be overwritten. This indicates that the "
                    "component_wise function is being called with an invalid set "
                    "of extractors in the extractor_order argument "
                    "that overlap in component indices."));

                component_order[i] = block_index;
              }
            // Increment block index
            block_index++;
          }
        return component_order;
      }
    } // namespace internal

  } // namespace utils

} // namespace FEValuesExtractors

/*------------------------------- Explicit Instantiations -------------*/

#include "fe/fe_values_extractors_utils.inst"

DEAL_II_NAMESPACE_CLOSE
