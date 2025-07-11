// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/mapping_c1.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <memory>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
MappingC1<dim, spacedim>::MappingC1()
  : MappingQ<dim, spacedim>(3)
{
  Assert(dim > 1, ExcImpossibleInDim(dim));
  Assert(dim == spacedim, ExcNotImplemented());
}



template <>
void
MappingC1<1>::add_line_support_points(const Triangulation<1>::cell_iterator &,
                                      std::vector<Point<1>> &) const
{
  const unsigned int dim = 1;
  (void)dim;
  Assert(dim > 1, ExcImpossibleInDim(dim));
}



template <>
void
MappingC1<2>::add_line_support_points(
  const Triangulation<2>::cell_iterator &cell,
  std::vector<Point<2>>                 &a) const
{
  const unsigned int          dim = 2;
  const std::array<double, 2> interior_gl_points{
    {0.5 - 0.5 * std::sqrt(1.0 / 5.0), 0.5 + 0.5 * std::sqrt(1.0 / 5.0)}};

  // loop over each of the lines, and if it is at the boundary, then first get
  // the boundary description and second compute the points on it. if not at
  // the boundary, get the respective points from another function
  for (unsigned int line_no = 0; line_no < GeometryInfo<dim>::lines_per_cell;
       ++line_no)
    {
      const Triangulation<dim>::line_iterator line = cell->line(line_no);

      if (line->at_boundary())
        {
          // first get the normal vectors at the two vertices of this line
          // from the boundary description
          const Manifold<dim> &manifold = line->get_manifold();

          Manifold<dim>::FaceVertexNormals face_vertex_normals;
          manifold.get_normals_at_vertices(line, face_vertex_normals);

          // then transform them into interpolation points for a cubic
          // polynomial
          //
          // for this, note that if we describe the boundary curve as a
          // polynomial in tangential coordinate @p{t=0..1} (along the line)
          // and @p{s} in normal direction, then the cubic mapping is such
          // that @p{s = a*t**3 + b*t**2 + c*t + d}, and we want to determine
          // the interpolation points at @p{t=0.276} and @p{t=0.724}
          // (Gauss-Lobatto points). Since at @p{t=0,1} we want a vertex which
          // is actually at the boundary, we know that @p{d=0} and @p{a=-b-c},
          // which gives @p{s(0.276)} and @p{s(0.726)} in terms of @p{b,c}. As
          // side-conditions, we want that the derivatives at @p{t=0} and
          // @p{t=1}, i.e. at the vertices match those returned by the
          // boundary.
          //
          // The task is then first to determine the coefficients from the
          // tangentials. for that, first rotate the tangents of @p{s(t)} into
          // the global coordinate system. they are @p{A (1,c)} and @p{A
          // (1,-b-2c)} with @p{A} the rotation matrix, since the tangentials
          // in the coordinate system relative to the line are @p{(1,c)} and
          // @p{(1,-b-2c)} at the two vertices, respectively. We then have to
          // make sure by matching @p{b,c} that these tangentials are
          // orthogonal to the normals returned by the boundary object
          const Tensor<1, 2> coordinate_vector =
            line->vertex(1) - line->vertex(0);
          const double h = std::sqrt(coordinate_vector * coordinate_vector);
          Tensor<1, 2> coordinate_axis = coordinate_vector;
          coordinate_axis /= h;

          const double alpha =
            std::atan2(coordinate_axis[1], coordinate_axis[0]);
          const double c = -((face_vertex_normals[0][1] * std::sin(alpha) +
                              face_vertex_normals[0][0] * std::cos(alpha)) /
                             (face_vertex_normals[0][1] * std::cos(alpha) -
                              face_vertex_normals[0][0] * std::sin(alpha)));
          const double b = ((face_vertex_normals[1][1] * std::sin(alpha) +
                             face_vertex_normals[1][0] * std::cos(alpha)) /
                            (face_vertex_normals[1][1] * std::cos(alpha) -
                             face_vertex_normals[1][0] * std::sin(alpha))) -
                           2 * c;

          const double t1   = interior_gl_points[0];
          const double t2   = interior_gl_points[1];
          const double s_t1 = (((-b - c) * t1 + b) * t1 + c) * t1;
          const double s_t2 = (((-b - c) * t2 + b) * t2 + c) * t2;

          // next evaluate the so determined cubic polynomial at the points
          // 1/3 and 2/3, first in unit coordinates
          const Point<2> new_unit_points[2] = {Point<2>(t1, s_t1),
                                               Point<2>(t2, s_t2)};
          // then transform these points to real coordinates by rotating,
          // scaling and shifting
          for (const auto &new_unit_point : new_unit_points)
            {
              Point<2> real_point(std::cos(alpha) * new_unit_point[0] -
                                    std::sin(alpha) * new_unit_point[1],
                                  std::sin(alpha) * new_unit_point[0] +
                                    std::cos(alpha) * new_unit_point[1]);
              real_point *= h;
              real_point += line->vertex(0);
              a.push_back(real_point);
            }
        }
      else
        // not at boundary, so just use scaled Gauss-Lobatto points (i.e., use
        // plain straight lines).
        {
          // Note that the zeroth Gauss-Lobatto point is a boundary point, so
          // we push back mapped versions of the first and second.
          a.push_back((1.0 - interior_gl_points[0]) * line->vertex(0) +
                      (interior_gl_points[0] * line->vertex(1)));
          a.push_back((1.0 - interior_gl_points[1]) * line->vertex(0) +
                      (interior_gl_points[1] * line->vertex(1)));
        }
    }
}



template <int dim, int spacedim>
void
MappingC1<dim, spacedim>::add_line_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  std::vector<Point<spacedim>> &) const
{
  DEAL_II_NOT_IMPLEMENTED();
}



template <>
void
MappingC1<1>::add_quad_support_points(const Triangulation<1>::cell_iterator &,
                                      std::vector<Point<1>> &) const
{
  const unsigned int dim = 1;
  (void)dim;
  Assert(dim > 2, ExcImpossibleInDim(dim));
}



template <>
void
MappingC1<2>::add_quad_support_points(const Triangulation<2>::cell_iterator &,
                                      std::vector<Point<2>> &) const
{
  const unsigned int dim = 2;
  (void)dim;
  Assert(dim > 2, ExcImpossibleInDim(dim));
}



template <int dim, int spacedim>
void
MappingC1<dim, spacedim>::add_quad_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  std::vector<Point<spacedim>> &) const
{
  DEAL_II_NOT_IMPLEMENTED();
}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingC1<dim, spacedim>::clone() const
{
  return std::make_unique<MappingC1<dim, spacedim>>();
}



// explicit instantiations
#include "fe/mapping_c1.inst"


DEAL_II_NAMESPACE_CLOSE
