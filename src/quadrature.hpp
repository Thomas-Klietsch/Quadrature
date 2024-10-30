// Copyright (c) 2024 Thomas Klietsch, all rights reserved.
//
// Licensed under the GNU Lesser General Public License, version 3.0 or later
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or ( at your option ) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General
// Public License along with this program.If not, see < https://www.gnu.org/licenses/>. 

#pragma once

#include <algorithm>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>

// C++23
#if __STDCPP_FLOAT128_T__ == 1
#include <stdfloat>
using Real = std::float128_t;
#else
using Real = long double;
#endif

std::string RealToString(
	Real const& value,
	uint8_t const& decimals = 8)
{
	constexpr auto max_digits{ std::numeric_limits<Real>::digits10 + 1 };

	// Decimals = '0' outputs all available decimals,
	// unlike std::setprecision(0) which outputs none.
	// A prefix space is added if value is positive (or zero)
	std::stringstream stream{
		decimals ?
		std::format("{: .{}f}", value, std::min<uint8_t>(decimals, max_digits)) :
		std::format("{: }", value)
	};
	return stream.str();
};

// Currently no support in C++23 for std::cout of std::float128_t
std::ostream& operator<<(
	std::ostream& os,
	Real const& value)
{
	return os << RealToString(value, std::cout.precision());
};

constexpr Real pi = std::numbers::pi_v<Real>;

// Return 'Not a Number', without throwing an exception
constexpr Real NaN = std::numeric_limits<long double>::quiet_NaN();

// Numeric stability
// Smallest value such that 1+epsilon evaluates to 1
constexpr Real numeric_epsilon = std::numeric_limits<Real>::epsilon();
// Smallest interval a function/integral will be evaluated in
Real static const numeric_interval = std::numeric_limits<double>::epsilon();

// Adaptive numerical integration of a function from 'a' to 'b',
// limited by recursive depth
namespace Quadrature
{
	// Algorithm 103
	// Simpson's rule integrator
	// Guy F. Kuncir
	//
	// Returns NaN if f(a),f(b) or f(a/2 + b/2), is NaN,
	Real Simpson(
		std::function<Real(Real)> const& function,
		Real a,
		Real b,
		Real const& a_epsilon = 1e-10,
		uint8_t const& a_max_depth = 8)
	{
		if (b < a)
			std::swap(a, b);

		// Maximal intervals to evaluate, 2^22 ~= 4 million
		uint8_t const max_depth = std::min(a_max_depth, static_cast<uint8_t>(22));

		// Epsilon is divided in recursions, so ensure at least some recursions
		Real const epsilon = std::max(a_epsilon, 512 * numeric_epsilon);

		struct Data
		{
			Real x{ 0 };
			Real y{ 0 }; // f(x)
			Real area{ 0 }; // F(x)[a;b]
			Data() {};
			Data(Real const& x, Real const& y, Real const& area = 0)
				: x(x), y(y), area(area) {
			};
		};

		// Simpson's rule, three point area approximation
		// A3,j = (b-a)(g0 + 4g2 + g4)/(3*2^(n+1))
		auto evaluate = [](
			std::function<Real(Real)> function,
			Data const& start,
			Data const& end
			) -> Data
		{
			Data middle;
			middle.x = (start.x + end.x) / 2;
			middle.y = function(middle.x);
			middle.area = std::abs(end.x - start.x) * (start.y + 4 * middle.y + end.y) / 6;
			return middle;
		};

		auto recursive = [&evaluate, &max_depth](
			// Self reference, needed for recursion, C++23
			this auto const& meta,
			std::function<Real(Real)> const& function,
			Data const& start,
			Data const& middle, // A3,j = A[0]
			Data const& end,
			Real const& epsilon,
			uint8_t depth) -> Real
		{
			if ((epsilon < numeric_epsilon) || (std::abs(end.x - start.x) < numeric_interval))
				return middle.area;

			// ^ y
			// |
			// +--*-----*----*------*-----*---> x
			//    start left middle right end
			//
			// Simpson's rule, five point area approximation
			// A5,j = (b-a)(g0 + 4g1 + 2g2 + 4g3 + g4)/(3*2^(n+2))
			//      = A[1] + A[2] = left + right
			auto const left = evaluate(function, start, middle);
			auto const right = evaluate(function, middle, end);

			if (!std::isfinite(left.y) || !std::isfinite(right.y))
				return NaN;

			// | (A5,j-A3,j)/A5,j | <= epsilon / 2^n
			// Estimated error
			// Real const error = (left.area + right.area - middle.area) / (left.area + right.area);
			// if ((std::abs(error) < epsilon) || (++depth > max_depth))
			// 	return left.area + right.area;

			// J. N. Lyness
			// Notes on the Adaptive Simpson Quadrature Routine
			// Estimated error using modification 1 and 2
			Real const error = (left.area + right.area - middle.area) / 15;
			if ((std::abs(error) < epsilon) || (++depth > max_depth))
				return left.area + right.area + error;

			return meta(function, start, left, middle, epsilon / 2, depth) +
				meta(function, middle, right, end, epsilon / 2, depth);
		};

		Data const start(a, function(a));
		Data const end(b, function(b));
		Data const middle = evaluate(function, start, end);

		if (!std::isfinite(start.y) || !std::isfinite(end.y) || !std::isfinite(middle.y))
			return NaN;

		return recursive(function, start, middle, end, epsilon, 0);
	};

	// Adaptive Quadrature - Revisited
	// Walter Gander, Walter Gautschi
	Real Lobatto(
		std::function<Real(Real)> const& function,
		Real a,
		Real b,
		Real const& a_epsilon = 1e-10,
		uint8_t const& a_max_depth = 2)
	{
		Real static const node_lobatto = std::sqrt(Real(1) / Real(5));
		Real static const node_kronrod = std::sqrt(Real(2) / Real(3));

		if (b < a)
			std::swap(a, b);

		// Maximal intervals to evaluate, 7^8 ~= 6 million
		uint8_t const max_depth = std::min(a_max_depth, static_cast<uint8_t>(8));

		Real const epsilon = std::max(a_epsilon, numeric_epsilon);

		struct Data
		{
			Real x{ 0 };
			Real y{ 0 }; // f(x)
			Data() {};
			Data(Real const& x, Real const& y = 0)
				: x(x), y(y) {
			};
		};

		auto recursive = [&](
			// Self reference, needed for recursion, C++23
			this auto const& meta,
			std::function<Real(Real)> const& function,
			Data const& start, // point 1
			Data const& end, // point 7
			uint8_t depth) -> Real
		{
			Real const h = (end.x - start.x) / 2;

			Data p4((start.x + end.x) / 2); // Middle point

			Data p2(p4.x - node_kronrod * h);
			Data p3(p4.x - node_lobatto * h);
			Data p5(p4.x + node_lobatto * h);
			Data p6(p4.x + node_kronrod * h);

			p2.y = function(p2.x);
			p3.y = function(p3.x);
			p4.y = function(p4.x);
			p5.y = function(p5.x);
			p6.y = function(p6.x);

			// Seven point area approximation
			Real const area_kronrod = (h / 1470) *
				((start.y + end.y) * 77 + (p2.y + p6.y) * 432 + (p3.y + p5.y) * 625 + p4.y * 672);

			if (!std::isfinite(area_kronrod))
				return NaN;

			if ((std::abs(h) < numeric_interval) || (++depth > max_depth))
				return area_kronrod;

			// Four point area approximation
			Real const area_lobatto = (h / 6) * (start.y + end.y + (p3.y + p5.y) * 5);

			// Error estimate
			if (std::abs(area_kronrod - area_lobatto) < epsilon)
				return area_kronrod;

			return meta(function, start, p2, depth) +
				meta(function, p2, p3, depth) +
				meta(function, p3, p4, depth) +
				meta(function, p4, p5, depth) +
				meta(function, p5, p6, depth) +
				meta(function, p6, end, depth);

		};

		Data const start(a, function(a));
		Data const end(b, function(b));

		if (!std::isfinite(start.y) || !std::isfinite(end.y))
			return NaN;

		return recursive(function, start, end, 0);
	};

};
