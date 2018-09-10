/*
  Copyright (C) 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_unit_tests_common_h
#define _aspect_unit_tests_common_h

#include <catch.hpp>

using Catch::Matchers::Contains;

/**
 * Compare the given two std::vector<double> entries with an epsilon (using Catch::Approx)
 */
inline void compare_vectors_approx(
  const std::vector<double> &computed,
  const std::vector<double> &expected)
{
  REQUIRE(computed.size() == expected.size());
  for (unsigned int i=0; i<computed.size(); ++i)
    {
      INFO("array index i=" << i << ": ");
      REQUIRE(computed[i] == Approx(expected[i]));
    }
}

#endif
