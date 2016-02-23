#include "dist.h"

#include <cstdio>

int main(int argc, char** argv)
{
	DistMap<bool> res =
	/*bernoulli (0.5) >> [&](auto const& c1) {
		return bernoulli (0.5) >>[&](auto const & c2) {
			if (c1 || c2)
				return ret(c1);
			else
				return empty<bool>();
		};
	};*/
	bernoulli (0.5);
	auto l1 = [&] (auto c) {
	       return ret(c?2.0:0.0);
	};
	using F = decltype(l1);
       	Bind<bool,DistMap<bool>,double,F> res2 = res >> l1;
	using DT = Dist<double,decltype(res2)>;
	DistMap<double> res3(static_cast<DT const&>(res2));
	double prob = (res3).integrate([&](auto beh){
			return beh;
				});

	printf("%lf\n", prob);
}

