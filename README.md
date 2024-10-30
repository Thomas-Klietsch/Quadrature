__Introduction__

A single header library for basic numeric integration, intended for 128bit floats.

If `std::float128_t` is not available/supported, `long double` will be used.

__Usage__

Define the equation either in a function or lambda:

    Real Function(Real const& x) {return std::sin(x);};

	auto lambda = [](Real const& x) -> Real {return std::sin(x);};

Perform the calculation with either Simpson's method or Lobatto, with an interval:

	std::cout << Quadrature::Simpson(Function, 0, pi) << "\n";
    std::cout << Quadrature::Lobatto(lambda, 0, pi) << "\n";

For increased accuracy of the quadrature, `epsilon` and recursive `max_depth` can be set.

	std::cout << Quadrature::Simpson(Function, 0, pi, 1e-16, 12) << "\n";

__Dependencies__

- C++23
