#include <unordered_map>
#include <cmath>
typedef double tprob;

template<typename F>
struct f_trait : public f_trait<decltype(std::function(*(F*)(nullptr)))>{};

template<typename R,typename A>
struct f_trait<std::function<R(A)>> {
	
    typedef R result_type;
    enum { arity = 1 };

    template <size_t i>
    struct arg
    {
        typedef A type;// std::tuple_element<i, std::tuple<A...>>::type type;
    };
};


template<int I> struct rank : public rank <I-1>{};
template<> struct rank<0> {};

namespace tmath {
	template<class T>
	static T tzero() {
		return T(0);
	}
	
	template<class T>
	static auto tmul(T&& v,tprob scale,rank<2>) -> decltype(std::forward<T>(v*=scale)) {
		return std::forward<T>(v*=scale);
	}

	template<class T>
	static auto tmul(T&& v,tprob scale,rank<1>) -> decltype(std::forward<T>(v=v*scale)) {
		return std::forward<T>(v= v * scale);
	}

	template<class T>
	static auto tmul(T&& v,tprob scale,rank<0>) -> decltype(std::forward<T>(v=scale*v)) {
		return std::forward<T>(v=scale*v);
	}

	template<class T>
	static auto tmul(T&& v,tprob scale) -> decltype (std::forward<T>(tmul(v,scale,rank<2>()))) {
		return std::forward<T>(tmul(v,scale,rank<2>()));
	}

	template<class T, class U>
	static auto tadd(T&& lhs,U const &rhs,rank<1>) ->decltype(std::forward<T>(lhs+=rhs)) {
		return std::forward<T>(lhs+=rhs);
	}

	template<class T, class U>
	static auto tadd(T&& lhs,U const &rhs,rank<0>) -> decltype (std::forward<T>(lhs=lhs+rhs)){
		return std::forward<T>(lhs=lhs+rhs);
	}

	template<class T, class U>
	static auto tadd(T&& lhs,U const & rhs) -> decltype(std::forward<T>(tadd(std::forward<T>(lhs),rhs,rank<1>()))) {
		return std::forward<T>(tadd(std::forward<T>(lhs),rhs,rank<1>()));
	}

	template<class T, class U>
	static auto tfma(T&& acc,U const & rhs,tprob scale, rank<1>) -> decltype(std::fma(rhs,scale,acc),acc) {
		acc=std::fma(rhs,scale,acc);
		return acc;
	}

	template<class T, class U>
	static auto tfma(T&& acc,U && rhs,tprob scale, rank<0>) ->decltype(std::forward<T>(tadd(std::forward<T>(acc),std::forward<U>(tmul(std::forward<U>(rhs),scale))))) {
		return std::forward<T>(tadd(std::forward<T>(acc),std::forward<U>(tmul(std::forward<U>(rhs),scale))));
	}

	template<class T, class U>
	static auto tfma(T&& acc,U const & rhs,tprob scale) -> decltype(std::forward<T>(tfma(std::forward<T>(acc),rhs,scale,rank<1>()))) {
		return std::forward<T>(tfma(std::forward<T>(acc),rhs,scale,rank<1>()));
	}

	template<class T>
	static auto tdiv(T&& acc, tprob scale, rank<2>) -> decltype(std::forward<T>(acc /= scale)) {
		return std::forward<T>(acc /= scale);
	}

	template<class T>
	static auto tdiv(T&& acc, tprob scale, rank<1>) -> decltype(std::forward<T>( acc = acc / scale)) {
		return std::forward<T>( acc = acc / scale);
	}

	template<class T>
	static auto tdiv(T&& acc, tprob scale, rank<0>) -> decltype(std::forward<T>(tmul(acc,tprob(1.0)/scale))){
		return std::forward<T>(tmul(acc,tprob(1.0)/scale));
	}

	template<class T>
	static auto tdiv(T&& acc,tprob scale) -> decltype(std::forward<T>(tdiv(acc,scale,rank<2>()))) {
		return std::forward<T>(tdiv(acc,scale,rank<2>()));
	}
};

template<class T,class E>
struct Dist {
	using domain_type=T;

	template<class F>
	auto operator()(F cont) {
		return std::move(*this) >> [cont](auto x){return ret(cont(x));};
	}

	template<class F>
	auto operator[](F guard) {
		return std::move(*this) >> [guard](auto x){return guard(x)?ret(x):empty(x);};
	}

	tprob norm() const {
		return norm(rank<1>());
	}
	
	tprob norm(rank<1>) const {
		return static_cast<E const&>(*this).norm();
	}

	tprob norm(rank<0>) const {
		return integrate([](T){return tprob(1.0);});
	}

	template<class U>
	auto integrate(U&& fun) const {
		return static_cast<const E&>(*this).integrate(std::forward<U>(fun));
	}
};

template <class T>
struct DistMap:public Dist<T,DistMap<T>> {
	std::unordered_map<T,tprob> weights;
	tprob total_weight;
	tprob norm() const {
		return total_weight;
       	}
	template<class F>
	auto integrate(F&& fun) const {
		using U = typename f_trait<F>::result_type;
		auto x =tmath::tzero<U>();
		for (auto const &e : weights) {
			auto v = fun(e.first);
			tmath::tmul(v,e.second);
			tmath::tadd(x,v);
			//tmath::tfma(x,v,e.second,rank<1>());
		}
		return x;
	}

	DistMap<T>& operator*=(tprob scale) {
		for(auto& e: weights)
			tmath::tmul(e.second,scale);
		tmath::tmul(total_weight,scale);
		return *this;
	}

	DistMap<T>& operator+=(DistMap<T> const& rhs) {
		size_t target_size = weights.size() + rhs.weights.size();
		if(weights.max_load_factor() * weights.bucket_count() < target_size)
			weights.reserve(target_size);
		for(auto& v: rhs.weights)
		{
			add_util(v.first,v.second);
		}
		tmath::tadd(total_weight,rhs.total_weight);
		return *this;
	}
	template<class E>
	DistMap(const Dist<T,E> &d) {
		auto l = [](const T& x){ return ret(x);}; 
		auto r = d.integrate(l);
		weights = std::move(r.weights);
		total_weight = std::move(r.total_weight);
	}

	template<class E>
	DistMap& operator =(Dist<T,E> &&d) {
		auto l = [](const T& x){ return ret(x);}; 
		auto r = d.integrate(l);
		weights = std::move(r.weights);
		total_weight = std::move(r.total_weight);
		return *this;
	}

	//DistMap(const T& val):weights{{val,1.0}},total_weight{1.0} {}

	DistMap(DistMap&& mv) = default;//:weights(std::move(mv.weights)),total_weight(mv.total_weight) {}
	DistMap(DistMap const& cp) = default;//:weights(std::move(mv.weights)),total_weight(mv.total_weight) {}
	DistMap& operator = (DistMap&& mv) = default;
	DistMap& operator = (DistMap const & cp) = default;
	/*{
		weights = std::move(mv.weights);
		total_weight(mv.total_weight);
		return *this;
	}*/

	template<typename U>
	DistMap(std::piecewise_construct_t,U&& map):weights(std::forward<U>(map)),total_weight{tmath::tzero<tprob>()} {
		for(auto const & p : weights)
			tmath::tadd(total_weight,p.second);
	}

	DistMap(std::piecewise_construct_t,std::initializer_list<std::pair<const T,tprob>> const & init_list):weights(init_list),total_weight{tmath::tzero<tprob>()} {
		for(auto const & p : weights)
			tmath::tadd(total_weight,p.second);
	}

	//zero
	DistMap():total_weight(0) {}
//	DistMap(int x): total_weight(0){ if (x != 0) throw "meh"; }


	static auto ret(T const & x) { return DistMap<T>{std::piecewise_construct,{{x,1.0}} }; }
private:
	void add_util(T const & k,tprob weight) {
		auto it = weights.find(k);
		if (it!=weights.end())
			tmath::tadd(it->second,weight);
		else
			weights.emplace(std::piecewise_construct,std::forward_as_tuple(k),std::forward_as_tuple(weight));
	}

};

template<class T,class E>
struct Normalized: Dist<T,Normalized<T,E>> {
	DistMap<T> const  unnorm;
	
	Normalized(Dist<T,E> &&v): unnorm(std::move(v)){}

	template<typename F>
	auto integrate(F&& fun) const {
		auto v = unnorm.integrate(fun);
		return tmath::tdiv(v,unnorm.norm());
	}

	tprob norm() const { return 1.0; }
};

template<class T,class E>
auto operator ~(Dist<T,E> &&v) { return Normalized<T,E>(std::move(v)); }

template<class T,class E1,class E2>
struct Sum: Dist<T,Sum<T,E1,E2>> {
	DistMap<T> const lhs;
	DistMap<T> const rhs;
	Sum(Dist<T,E1> && tb, Dist<T,E2> && fb):lhs(std::move(tb)),rhs(std::move(fb)) {} 

	template<class F>
	auto integrate(F&& fun) const {

		auto v = lhs.integrate(fun);
		return std::move(tmath::tadd(v,rhs.integrate(fun)));
	}

	tprob norm() const { return lhs.norm() + rhs.norm(); }
};

template<class T,class E1,class E2>
auto operator +(Dist<T,E1> && lhs,Dist<T,E2> && rhs) { return Sum<T,E1,E2>(std::move(lhs),std::move(rhs)); }

template<class T,class E>
struct Mul: Dist<T,Mul<T,E>> {
	DistMap<T> const lhs;
	tprob scale;
	Mul(Dist<T,E> && tb, tprob fb):lhs(std::move(tb)),scale(fb) {} 

	template<class F>
	auto integrate(F&& fun) const {
		auto v = lhs.integrate(fun);
		return std::move(tmath::tmul(v,scale));
	}
	tprob norm() const { return lhs.norm() * scale; }

};
template<class T,class E>
auto operator *(Dist<T,E> &&lhs,tprob scale) { return Mul<T,E>(std::move(lhs),scale); }

template<class T,class E>
auto operator *(tprob scale,Dist<T,E> &&lhs) { return Mul<T,E>(std::move(lhs),scale); }

template<class T,class E1,class E2>
auto mix(tprob v, Dist<T,E1> lhs,Dist<T,E2> rhs) {
       	return ~((lhs * v) + (rhs * (tprob(1.0) - v)));
}

template<class T,class E,class U, class F>
struct Bind: public Dist<U,Bind<T,E,U,F>>
{
	DistMap<T> const d;
	F cont;
	Bind(Dist<T,E> && db,F tcont):d(std::move(db)),cont(tcont){}
	template<class G>
	auto integrate(G&& fun) const {
		return d.integrate([&](T const & x){
			       	return cont(x).integrate(fun);
				});
	}

};

template<class T,class E, class F>
auto bind(Dist<T,E> && source,F cont) {
	//typename f_trait<F>::result_type
	using U = decltype(cont(*(T*)nullptr).integrate([](auto const &y) {return y;}));
	//using U = decltype(source.integrate([&](T const & x){ return cont(x).integrate([](auto const & y){return y;});  }));
	return Bind<T,E,U,F>(std::move(source),cont);
}

template<class T,class E, class F>
auto operator>>(Dist<T,E> && source,F cont) {
	return bind(std::move(source),cont);
}


template<class T>
auto ret(T x) { return DistMap<T>::ret(x); }

using Bernoulli = DistMap<bool>;

Bernoulli bernoulli(tprob p) {
	tprob q = tprob(1);
	tmath::tfma(q,p,tprob(-1));
	return {std::piecewise_construct,{{true,p},{false,q}}};
}

template<class T,class It>
DistMap<T> uniform(It begin, It end) {
	ptrdiff_t count = end - begin;
	tprob scale = tprob(1);
	tmath::tdiv(scale,tprob(count));
	std::unordered_map<T,tprob> weights;
	weights.reserve(count);
	for(auto it = begin; it != end; it++)
	{
		weights.emplace(std::piecewise_construct,std::forward_as_tuple(*it),std::forward_as_tuple(scale));
	}
	return {std::piecewise_construct,weights};
}

template<class T>
DistMap<T> uniform(std::initializer_list<T> init) {
	return uniform<T,decltype(init.begin())>(init.begin(),init.end());
}
template<class T>
DistMap<T> empty(T= tmath::tzero<T>()) {
	return DistMap<T>();
}
