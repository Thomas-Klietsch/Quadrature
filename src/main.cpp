
#if __STDCPP_FLOAT128_T__ == 1
#include <stdfloat>
#endif

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

#include "./quadrature.hpp"

Real Function(Real const& x) {
	return std::sin(x);
};

int main(int argc, char* argv[])
{

	std::cout << std::setprecision(20) << std::endl;

	auto lambda = [](Real const& x) -> Real
	{
		return std::sin(x);
	};

	std::cout << "f(x)=sin(x), x=[0;pi]\n";
	std::cout << "Exact value: " << Real(-std::cos(pi)) - Real(-std::cos(0)) << "\n";
	std::cout << "Simpson:     " << Quadrature::Simpson(Function, 0, pi) << "\n";
	std::cout << "Lobatto:     " << Quadrature::Lobatto(lambda, 0, pi) << "\n";

	auto func_poly = [](Real const& x) -> Real
	{
		return 6 * x * x - 8 * x + 5;
	};

	std::cout << "\nf(x)=6x^2-8x+5, x=[1;4]\n";
	// 2 x^3 - 4 x^2 + 5 x + constant
	std::cout << "Exact value: " << Real(2 * 4 * 4 * 4 - 4 * 4 * 4 + 5 * 4) - Real(2 * 1 * 1 * 1 - 4 * 1 * 1 + 5 * 1) << "\n";
	std::cout << "Simpson:     " << Quadrature::Simpson(func_poly, 1, 4) << "\n";
	std::cout << "Lobatto:     " << Quadrature::Lobatto(func_poly, 1, 4) << "\n";

	auto func_log = [](Real const& x) -> Real
	{
		return std::log(x);
	};

	std::cout << "\nf(x)=ln(x), x=[1;2]\n";
	// x * ln(x) - x
	std::cout << "Exact value: " << (2 * std::log(Real(2)) - 2) - (1 * std::log(Real(1)) - 1) << "\n";
	std::cout << "Simpson:     " << Quadrature::Simpson(func_log, 1, 2) << "\n";
	std::cout << "Lobatto:     " << Quadrature::Lobatto(func_log, 1, 2) << "\n";

	auto func_sqrt = [](Real const& x) -> Real
	{
		return std::sqrt(x) + 1 / (3 * std::sqrt(x));
	};

	std::cout << "\nf(x)=sqrt(x)+1/(3*sqrt(x)), x=[4;9]\n";
	// 2/3 * (2 * sqrt(x) * (1 + x))
	std::cout << "Exact value: " << Real(40) / Real(3) << "\n";
	std::cout << "Simpson:     " << Quadrature::Simpson(func_sqrt, 4, 9) << "\n";
	std::cout << "Lobatto:     " << Quadrature::Lobatto(func_sqrt, 4, 9) << "\n";

	std::cout << "\nf(x)=x^i, x=[0;1]\n";
	// x^(1+i) / (1+i)
	for (uint8_t i{ 0 };i < 5;++i)
	{
		auto func_capture = [i](Real const& x) -> Real
		{
			return std::pow(x, i);
		};
		std::cout << "x^" << i + 0 << ": " << Quadrature::Lobatto(func_capture, 0, 1) << "\n";
	};

};
