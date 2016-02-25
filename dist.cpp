#include "dist.h"

#include <cstdio>
#include <cstdlib>
	/*bernoulli (0.5) >> [&](auto const& c1) {
		return bernoulli (0.5) >>[&](auto const & c2) {
			if (c1 || c2)
				return ret(c1);
			else
				return empty<bool>();
		};
	};*/

int main(int argc, char** argv)
{

	/*
	auto f1 = [&](bool c1) {
		auto f2 = [&](bool c2){
			return DistMap<bool>(ret(c1&&c2));
		};
		using F2 = decltype(f2); // bool -> DistMap<bool>
		return Bind<bool,DistMap<bool>,bool,F2>(bernoulli(0.5),f2);
	};
	using F1 = decltype(f1); // bool -> Bind<bool,DistMap<bool>,bool,F2>

	auto finalDist = Bind<bool,DistMap<bool>,bool,F1>(bernoulli(0.5),f1);

	auto f3 = [&](bool beh){return beh?1.0:0.0; };
	using F3 = decltype(f3); // bool -> double

	//double prob = finalDist.integrate(f3);
	//double prob = Bind<bool,DistMap<bool>,bool,F1>(bernoulli(0.5),f1).integrate(f3);
	auto f4 = [&](bool const & c1){
		auto f22 = [&](bool c2){
			return DistMap<bool>(ret(c1&&c2));
		};
		using F22 = decltype(f22); // bool -> DistMap<bool>

		auto f5 = [&] (bool const & x) {
			return f22(x).integrate(f3);
		};
		return bernoulli(0.5).integrate(f5);
	};
	using T4 = decltype(f4); // bool -> Bind<bool,DistMap<bool>,bool,F2>.integrate(bool->double);
	double prob = bernoulli(0.5).integrate(f4);

	using U = decltype((l4(weights.begin()->first)));
		auto x =tmath::tzero<U>();
		for (auto const &e : weights) {
			auto v = fun(e.first);
			tmath::tmul(v,e.second,rank<2>());
			tmath::tadd(x,v);
			//tmath::tfma(x,v,e.second,rank<1>());
		}
		return x;


	printf("%lf\n", prob);
*/
	/*
	auto ff1 =
		[=](bool const & c1) {
			auto ff2 =
				[=](bool const &c2){
					DistMap<bool> rv = ret (c1 && c2);
					return rv;
				};
			using FF2 = decltype(ff2);
			Bind<bool,DistMap<bool>,bool,FF2> dret = bernoulli(0.5) >> ff2;
			return dret;
		}; 
	using FF1 = decltype(ff1);

	auto fffinta =  [=](bool const & c1){ auto ff4 = [](auto const &){ return ret(true); }; return Bind<bool,DistMap<bool>,bool,decltype(ff4)>(ret(true),ff4);};
	using FFF = decltype(fffinta);
	Bind<bool,DistMap<bool>,bool,FFF> d2 = bernoulli(0.5) >> fffinta; */
    /*auto d = bernoulli(0.5) >> [&](auto const &c1) {
	    return bernoulli(0.5) >> [&] (auto const &c2) {
		    //return ret(c1&&c2);
		    if (c1 || c2)
			    return ret(c1);
		    else
			    return empty<bool>();
	    };
    };
    //DistMap<bool> ciccio = d;
    auto normd = ~std::move(d);*/
	//double prob = (normd).integrate([](auto const &g){ return g?0.0:1.0; });
	//double prob2 =// d2.integrate([](bool const &x){bool y = x; return y;});
	/*bernoulli(0.5).integrate([&](auto const & x) {
			return bernoulli(0.5).integrate([&](auto const & x) {
				auto ff2 =
				[=](bool const &c2){
				DistMap<bool> rv = ret (c1 && c2);
				return rv;
				};
				using FF2 = decltype(ff2);
				DistMap<bool> dd = ff2(x);
				dd.integrate([](bool const &x){bool y = x; return y;});
				};
				});
*/

	srandomdev();

	auto prior = uniform<int>({-2,-1,0,1,2});
	for(int i = 0; i < 20000; i++)
	{
		long v = -2;
		for(int j = 0; j < 4; j++)
			v+=random() &1;
		
		DistMap<int> newPrior = ~(std::move(prior) >>[=](int p){
			double rejectProb = exp(-fabs((double)p - v)/10.0);
			//return ret(p);
			return bernoulli(rejectProb) >> [p](auto b){ return b?ret(p):empty(p);};
		});
		prior = newPrior;

	}
	
	for ( auto &p : prior.weights)
		printf("%d %lf\n", p.first, p.second);
	

}

void testBoh() {

	// la bind deve fare 1 -> true per ogni valore 
	DistMap<bool> b = bernoulli(0.1) >> [](auto b) { return ret(true); };
	for ( auto &p : b.weights)
		printf("%d %lf\n", p.first, p.second);

}

