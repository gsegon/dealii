// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2012 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef deal_ii_fe_values_extractors_utils_h
#define deal_ii_fe_values_extractors_utils_h

#include <deal.II/fe/fe_values_extractors.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace FEValuesExtractors
{
  /**
   * Namespace for utility functions related to FEValuesExtractors.
   */
  namespace utils
  {

    namespace internal
    {

      /**
       * @internal
       * Utility function used to construct a `component_order` vector by
       * processing a vector of FEValuesExtractors in given order. Supported
       * extractors are listed in FEValueExtractors::AnyExtractor.
       */
      template <int dim, int spacedim>
      std::vector<unsigned int>
      generate_component_order(
        const std::vector<FEValuesExtractors::AnyExtractor> &extractors,
        const unsigned int                                   fe_n_components);

    } // namespace internal

  } // namespace utils

} // namespace FEValuesExtractors

DEAL_II_NAMESPACE_CLOSE

#endif // deal_ii_fe_values_extractors_utils_h
