#include <vector>
#include <math.h>
#include <stdio.h>

#ifndef _CPP_
#include <emscripten/bind.h>
using namespace emscripten;
#endif

typedef double type1;

type1 tau_true=50, Ts_true=30;
//type1 mu_true = 100, lambda_true = 204.88;
type1 mu_true = 10, lambda_true = 8.1955;

type1 IG_pdf(type1 x, type1 mu, type1 lambda){
	return (x<=0) ? 0 : sqrt(lambda/2/M_PI/x/x/x)*exp(-lambda*(x-mu)*(x-mu)/2/mu/mu/x);
}

class Receiver {
	public:
		Receiver():tau_(tau_true), Ts_(Ts_true), k_(0), mu_(mu_true), lambda_(lambda_true){}

		void push_y(type1 y) {
			y_.push_back(y);
		}
		type1 get_y(int index){return y_[index];}
		int get_n(){return y_.size();}

		int get_k() const {return k_;}
		void increase_k(){++k_;}

		type1 get_tau() const {return tau_;}
		type1 get_Ts() const {return Ts_;}

		type1 mean(int start, int end) const {
			type1 ret = 0;
			for(int i=start; i<end; ++i)
				ret += y_[i];
			return ret/(end-start);
		}
		type1 moment2(int start, int end, type1 center) const {
			type1 ret = 0;
			for(int i=start; i<end; ++i){
				type1 tmp = y_[i] - center;
				ret += tmp*tmp;
			}
			return ret/(end-start);
		}
		type1 moment_1(int start, int end, type1 center) const {
			type1 ret = 0;
			for(int i=start; i<end; ++i){
				type1 tmp = y_[i] - center;
				ret += 1/tmp;
			}
			return ret/(end-start);
		}

		void LE_1(type1 mean_y0){
			tau_ = y_[0]-mean_y0;
		}

		void ini_tau(){
			int n = get_n();
			type1 y_bar = mean(0,n);
			type1 m2 = moment2(0,n,y_bar);
			type1 a = y_bar - y_[0];
			tau_ = y_[0]-a*a*a/2/m2/log(n);
			//fprintf(stderr,"ini_tau=%f\n",tau_);
		}


		type1 score(type1 s){
			//return lambda_/(2*mu_*mu_);
			//return 3/(2*s);
			//return lambda_/(2*s*s);
			return lambda_/(2*s*s) - 3/(2*s) - lambda_/(2*mu_*mu_);
		}

		void iter_tau(){
			int n = get_n(), iter_num = 0;
			type1 avg_s = 0, I0 = 1;
			do{
				type1 sum = 0, sum_2 = 0;
				tau_ -= avg_s/I0; //epsilon += so tau -=
				if(tau_ > y_[0]) tau_ = y_[0]-0.01;
				for(int i=0; i<n; ++i){
					type1 s = y_[i] - tau_;
					while(s > Ts_) {
						s -= Ts_;
					}
					type1 f0 = IG_pdf(s,mu_,lambda_);
					type1 p = f0/(f0+IG_pdf(s+Ts_,mu_,lambda_)+IG_pdf(s+Ts_+Ts_,mu_,lambda_));
					type1 sample_score = (score(s)*p + score(s+Ts_)*(1-p));
					sum += sample_score;
					//sum_2 += (sample_score*sample_score);
					//sum += score(s);
					//fprintf(stderr,"s=%f, score=%f\n", s, sample_score);
				}
				avg_s = sum/n;
				//I0 = sum_2/(n-1) - avg_s*avg_s*n/(n-1);
				//I0 = sum_2/n;
				I0 = lambda_/mu_/mu_/mu_ + 9/2/mu_/mu_ + 21/2/mu_/lambda_ + 21/2/lambda_/lambda_;
				iter_num ++;
			}while((avg_s > 0.0001 || avg_s < -0.0001) && iter_num < 10000);
			if(iter_num>100) fprintf(stderr,"iter_num=%d\n",iter_num);
			//fprintf(stderr,"I0=%f\n",I0);
			//fprintf(stderr,"avg_s=%f\n", avg_s);
			//if(tau_ > 55 || tau_ < 45) fprintf(stderr,"iter_tau=%f\n",tau_);
		}

	private:
		type1 tau_;
		type1 Ts_;
		int k_;
		type1 mu_;
		type1 lambda_;
		std::vector<type1> y_;
};

#ifndef _CPP_
// Binding code
EMSCRIPTEN_BINDINGS(my_class_example) {
	class_<Receiver>("Receiver")
		.constructor<>()
		.function("push_y", &Receiver::push_y)
		.function("get_y", &Receiver::get_y)
		.function("get_n", &Receiver::get_n)
		.property("k_", &Receiver::get_k)
		.function("increase_k",&Receiver::increase_k)
		.property("tau_", &Receiver::get_tau)
		.property("Ts_", &Receiver::get_Ts)
		.function("ini_tau", &Receiver::ini_tau)
		//.class_function("getStringFromInstance", &Receiver::getStringFromInstance)
		;
}
#endif

#ifdef _CPP_

#include <random>
#include <algorithm>    // std::sort

/*
type1 cdf_table[LARGE];
void IG_cdf_table(){
	type1 x = 0;
	type1 cdf = 0;
	for(int i=0; i<LARGE; ++i){
		x += SMALL;
		cdf += IG_pdf(x,mu_true,lambda_true);
		cdf_table[i] = cdf*SMALL;
	}
}
*/

//std::default_random_engine rd;
std::random_device rd;
std::normal_distribution<type1> normal(0.0,1.0);
std::uniform_real_distribution<type1> uni(0.0,1.0);

type1 get_IG_sample(){
	
	type1 mu = mu_true, lambda = lambda_true;
	type1 v = normal(rd);
	type1 yy = v*v;
	type1 xx = mu + (mu*mu*yy)/(2*lambda) - (mu/(2*lambda))*sqrt(4*mu*lambda*yy + mu*mu*yy*yy);
	type1 test = uni(rd); 
	type1 th = (mu)/(mu + xx);

	type1 ret = (test <= th)? xx : (mu*mu)/xx;
	

	/*
	type1 p = uni(rd);
	int b=0;
	while(b<LARGE && cdf_table[b++]<p);
	int a = b-1;
	type1 ret = 0;
	if(a==0) {
		type1 f_a = 0;
		type1 f_b = cdf_table[b-1];
		type1 w = (p-f_a)/(f_b-f_a);
		ret = a*w + b*(1-w);
	}
	else if(a==LARGE){
		type1 f_a = cdf_table[a-1];
		type1 f_b = 1;
		type1 w = (p-f_a)/(f_b-f_a);
		ret = a*w + b*(1-w);
	}
	else {
		type1 f_a = cdf_table[a-1];
		type1 f_b = cdf_table[b-1];
		type1 w = (p-f_a)/(f_b-f_a);
		ret = a*w + b*(1-w);
		//fprintf(stderr,"%f %d %f %d %f ",p,a,f_a,b,f_b);
	}
	ret = ret*SMALL;
	*/

	//fprintf(stderr,"IG_sample=%f\n",ret);
	return ret;
}

type1 offline(int q1){
	type1 delta = 0.001;
	type1 mean_y0 = 0;
	for(int j=1; j<=300*mu_true; ++j){
		type1 cdf = 0;
		for(int k=1; k<=j; ++k) {
			cdf += IG_pdf(delta*k, mu_true, lambda_true);
		}
		mean_y0 += pow(1-cdf*delta, q1);
	}
	//fprintf(stderr,"mean_y0=%f\n", mean_y0);
	return mean_y0*delta;
}

int main(int argc, char* argv[]){
	//fprintf(stderr,"hello world!\n");

	if(argc<3){
		fprintf(stderr,"Usage: %s <test_max> <K>\n", argv[0]);
		return 0;
	}

	//IG_cdf_table();

	type1 tau = tau_true, Ts=Ts_true;
	int q1 = 16, K = atoi(argv[2]);
	type1 mean_y0 = offline(q1);
	// 2_ary
	//std::uniform_int_distribution<int> uni_int(0,1);
	//int m_ary[] = {q1/2,q1*3/2};
	// 8_ary
	std::uniform_int_distribution<int> uni_int(0,7);
	int m_ary[] = {q1/8,q1*3/8,q1*5/8,q1*7/8,q1*9/8,q1*11/8,q1*13/8,q1*15/8};
	type1 err;

	int test_max = atoi(argv[1]);
	std::vector<type1> bias(2,0);
	std::vector<type1> mse_table(2,0);
	for(int test=0;test<test_max;test++){
		//fprintf(stderr,"test=%d\n",test);
		//tx
		std::vector<type1> x;
		for(int i=0; i<q1; ++i)
			x.push_back( tau + 0 + get_IG_sample());
			//x.push_back( tau + 0);
		for(int k=1; k<K; ++k){
			int q = m_ary[uni_int(rd)];
			for(int i=0; i<q; ++i)
				x.push_back( tau + k*Ts + get_IG_sample());
		}
		std::sort(x.begin(), x.end());

		//rx
		Receiver rx_c;
		for(int i=0; i<q1; ++i)
			rx_c.push_y(x[i]);
		//rx_c.ini_tau();
		rx_c.LE_1(mean_y0);
		err = (rx_c.get_tau() - tau);
		bias[0] += err;
		mse_table[0] += err*err;
		int N = x.size();
		for(int i=q1; i<N; ++i){
			rx_c.push_y(x[i]);
		}
		rx_c.iter_tau();
		err = (rx_c.get_tau() - tau);
		bias[1] += err;
		mse_table[1] += err*err;
	}
	fprintf(stderr, "q1=%d\n",q1);
	//printf("bias_ini_tau,\t%f\n", bias[0]/test_max);
	//printf("mse_ini_tau,\t%f\n", mse_table[0]/test_max);
	fprintf(stderr, "K\tmse_iter_tau\tbias_iter_tau\n");

	if(K==1){
		printf("1,\t%f,\t\n",mse_table[0]/test_max);
	}
	else{
		printf("%d,\t",K);
		printf("%f,\t", mse_table[1]/test_max);
		printf("%f,\t", bias[1]/test_max);
		printf("\n");
	}
	return 0;
}
#endif

