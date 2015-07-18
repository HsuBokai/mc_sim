#include <vector>
#include <math.h>
#include <stdio.h>

#ifndef _CPP_
#include <emscripten/bind.h>
using namespace emscripten;
#endif

typedef double type1;

type1 tau_true=0;
//type1 mu_true = 150, lambda_true = 51.22, Ts_true = 200;
type1 mu_true = 10, lambda_true = 8.1955, Ts_true = 30;;

type1 IG_pdf(type1 x, type1 mu, type1 lambda){
	return (x<=0) ? 0 : sqrt(lambda/2/M_PI/x/x/x)*exp(-lambda*(x-mu)*(x-mu)/2/mu/mu/x);
}

class Receiver {
	public:
		Receiver():tau_(0), Ts_(Ts_true), k_(0), mu_(1), lambda_(1){}

		void push_y(type1 y) {
			y_.push_back(y);
		}
		type1 get_y(int index){return y_[index];}
		int get_n(){return y_.size();}

		int get_k() const {return k_;}
		void increase_k(){++k_;}

		type1 get_tau() const {return tau_;}
		type1 get_Ts() const {return Ts_;}
		type1 get_mu() const {return mu_;}
		type1 get_lambda() const {return lambda_;}

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
			return ret/(end-start-1.0);
		}
		type1 moment_1(int start, int end, type1 center) const {
			type1 ret = 0;
			for(int i=start; i<end; ++i){
				type1 tmp = y_[i] - center;
				ret += 1.0/tmp;
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
			tau_ = y_[0]-(a*a*a)/(2.0*m2*log(n));
			//fprintf(stderr,"ini_tau=%f\n",tau_);
		}

		void ini_mu(){
			int n = get_n();
			type1 y_bar = mean(0,n);
			mu_ =  y_bar - tau_;
		}

		void ini_lambda(){
			int n = get_n();
			type1 m_1 = moment_1(0,n,tau_);
			lambda_ = 1.0/(m_1 - 1.0/mu_);
			//fprintf(stderr,"lambda=%f\n",lambda_);
		}


		type1 score(type1 s){
			return lambda_/(2*s*s) - 3/(2*s) - lambda_/(2*mu_*mu_);
		}
/*
		type1 epsilon_iter(type1 s){
			type1 ret = 0;
			type1 b = 1.0/lambda_+1.0/mu_-1.0/s;
			ret = score(s) + 3.0*lambda_*(1.0/mu_+1.0/lambda_)*b/2.0;
			return ret/(3.0/(2.0*mu_*lambda_) + 6.0/(lambda_*lambda_));
		}

		type1 lambda_iter(type1 s){
			type1 ret = 0;
			type1 b = 1.0/lambda_+1.0/mu_-1.0/s;
			ret = 3.0*lambda_*(1.0/mu_+1.0/lambda_)*score(s)+(9.0/(2.0*mu_*mu_) + 21.0/(2.0*mu_*lambda_) + 21.0/(2.0*lambda_*lambda_))*lambda_*lambda_*b;
			return ret/(3.0/(2.0*mu_*lambda_) + 6.0/(lambda_*lambda_));
		}
*/

		void iter_tau(){
			type1 init_tau = tau_;

			int n = get_n();
			type1 y_bar = mean(0,n);
			mu_ = y_bar - tau_ - Ts_*(k_-1)/2.0;

			int iter_num = 0;
			type1 max_iter_num = 1000.0;
			type1 smallest_s = 1000, best_tau = tau_;
			type1 smallest_s_lambda = 1000, best_lambda = lambda_;
			type1 avg_s = 0, I0 = 1.0;
			type1 m_1 = 0;
			type1 score_lambda = 0;
			//type1 m1 = 0;
			//type1 avg_epsilon_iter = 0, avg_lambda_iter = 0;
			do{
				//mu_ += avg_s/I0;
				//mu_ = mu_true;
				//mu_ += avg_epsilon_iter;

				//epsilon += so tau -=
				tau_ -= avg_s/I0;
				//tau_ = tau_true;
				//tau_ -= avg_epsilon_iter;
				//avg_epsilon_iter = 0;
				avg_s = 0;

				//lambda_ += score_lambda;
				//lambda_ = lambda_true;
				//lambda_ += avg_lambda_iter;
				//avg_lambda_iter = 0;
				score_lambda = 0;
				//m1 = 0;
				m_1 = 0;
				if(lambda_ < 0) lambda_ = iter_num/max_iter_num;
				if(mu_ < 0) mu_ = iter_num/max_iter_num;
				//if(tau_ > y_[0]) tau_ = y_[0]-0.01;
				if(tau_ > y_[0]) {
					//tau_ = y_[0]-Ts_/2.0-iter_num/max_iter_num;
					tau_ = init_tau;
				}
				if(tau_ < y_[0]-Ts_) {
					//tau_ = y_[0]-Ts_/2.0+iter_num/max_iter_num;
					tau_ = init_tau;
				}
				for(int i=0; i<n; ++i){
					type1 s = y_[i] - tau_;
					while(s > Ts_) {
						s -= Ts_;
					}
					type1 f0 = IG_pdf(s,mu_,lambda_);
					type1 f1 = IG_pdf(s+Ts_,mu_,lambda_);
					type1 f2 = IG_pdf(s+2*Ts_,mu_,lambda_);
					type1 p = f0/(f0+f1+f2);
					type1 sample_score = (score(s)*p + score(s+Ts_)*(1-p));
					avg_s += sample_score;
					m_1 += (p/s + (1.0-p)/(s+Ts_));
					//m1 += (s*f0+(s+Ts_)*f1+(s+2*Ts_)*f2)/(f0+f1+f2);
					//avg_epsilon_iter += (epsilon_iter(s)*p+epsilon_iter(s+Ts_)*(1-p));
					//avg_lambda_iter += (lambda_iter(s)*p+lambda_iter(s+Ts_)*(1-p));
				}
				m_1 /= n;
				//m1 /= n;
				score_lambda = (1.0/lambda_ + 1.0/mu_ - m_1)*lambda_*lambda_;
				avg_s = avg_s/n;
				I0 = lambda_/(mu_*mu_*mu_) + 9.0/(2.0*mu_*mu_) + 21.0/(2.0*mu_*lambda_) + 21.0/(2.0*lambda_*lambda_);
				//I0 = (3.0/(2.0*mu_*lambda_) + 6.0/(lambda_*lambda_));
				//I0 -= avg_s*avg_s;
				//avg_epsilon_iter = avg_epsilon_iter/n;
				//avg_lambda_iter = avg_lambda_iter/n;
				iter_num ++;
				//fprintf(stderr,"iter_num=%d\n", iter_num);
				mu_ = y_bar - tau_ - Ts_*(k_-1)/2.0;
				lambda_ = 1/(m_1 - 1/mu_);
				if(0 < avg_s && avg_s < smallest_s) {
					smallest_s = avg_s;
					best_tau = tau_;
				}
				if(0 < -avg_s && -avg_s < smallest_s){
					smallest_s = -avg_s;
					best_tau = tau_;
				}
				if(avg_s/I0 > Ts_ || avg_s/I0 < -Ts_) {
					//tau_ = best_tau;
					//avg_s = 0;
					//fprintf(stderr,"diverge! best_tau=%f\n", best_tau);
					//fprintf(stderr,"avg_epsilon_iter=%f\n",avg_epsilon_iter);
				}
				
				if(0 < score_lambda && score_lambda < smallest_s_lambda) {
					smallest_s_lambda = score_lambda;
					best_lambda = lambda_;
				}
				if(0 < -score_lambda && -score_lambda < smallest_s_lambda) {
					smallest_s_lambda = -score_lambda;
					best_lambda = lambda_;
				}
				if(score_lambda > 10000 || score_lambda < -10000) {
					lambda_ = best_lambda;
					score_lambda = 0;
					fprintf(stderr,"diverge! best_lambda=%f\n",best_lambda);
					//fprintf(stderr,"avg_lambda_iter=%f\n",avg_lambda_iter);
				}
				
				if(iter_num > max_iter_num){
					fprintf(stderr,"iter_num over %f\n", max_iter_num);
					tau_ = best_tau;
					break;
				}
				//fprintf(stderr,"I0=%f\n",I0);
				//fprintf(stderr,"avg_s/I0=%f\n", avg_s/I0);
				//fprintf(stderr,"iter_tau=%f\n",tau_);
				//fprintf(stderr, "m_1=%f\n",m_1);
				//fprintf(stderr, "mu_=%f\n",mu_);
				//fprintf(stderr,"iter_lambda=%f\n",lambda_);
			//}while(avg_s > 0.00001 || avg_s < -0.00001 || score_lambda > 0.01 || score_lambda < -0.01);
			}while(avg_s > 0.00001 || avg_s < -0.00001);
			//fprintf(stderr,"I0=%f\n",I0);
			//fprintf(stderr,"score_lambda=%f\n", score_lambda);
			//fprintf(stderr,"avg_s/I0=%f\n", avg_s/I0);
			//if(tau_ > 55 or tau_ < 45) fprintf(stderr,"iter_tau=%f\n",tau_);
			//fprintf(stderr,"mu_=%f\n",mu_);
			//fprintf(stderr,"lambda_=%f\n",lambda_);
		}

		void est_Ts(int q1){
			int start = 0;
			type1 A = mean(0,q1);
			type1 B = 0;
			for(int i=2; i<=k_; ++i){
				start += q1;
				type1 y_bar_new = mean(start, start+q1);
				A = A*(i-1)/i + y_bar_new/i;
				B = B + (i-1)*y_bar_new;
			}
			Ts_ = 12*B/(k_*(k_-1)*(k_+1)) - 6*A/(k_+1);
			//fprintf(stderr,"Ts_ = %f\n",Ts_);
		}

#define LARGE 10*mu_true
#define DELTA 0.01

		void pitman(){
			int n = get_n();
			type1 denominator = 0, nominator = 0;
			for(type1 u = 0; u<LARGE; u+= DELTA){
				//fprintf(stderr,"u=%f\n",u);
				long double value = 1;
				for(int i=0; i<n; ++i) {
					long double ig_large = IG_pdf(u+y_[i]-y_[0], mu_true, lambda_true)*n/50;
					value *= ig_large;
					//fprintf(stderr,"value=%f\n",value);
				}
				denominator += value;
				nominator += value*u;
				//fprintf(stderr,"deno_ = %f\n",denominator);
				//fprintf(stderr,"no_ = %f\n",nominator);
			}
			//denominator *= DELTA;
			//nominator *= DELTA;
			tau_ = y_[0] - nominator/denominator;
			//fprintf(stderr,"y0=%f\n",y_[0]);
			//fprintf(stderr,"tau_ = %f\n",tau_);
		}

		type1 pitman_scale(){
			int n = get_n();
			long double denominator = 0, nominator = 0;
			for(type1 u = 0; u<LARGE; u+= DELTA){
				//fprintf(stderr,"u=%f\n",u);
				long double value = 1;
				for(int i=0; i<n; ++i) {
					long double ig_large = exp(u)*y_[i]*IG_pdf(exp(u)*y_[i]/y_[0], mu_true, lambda_true);
					value *= ig_large;
					//fprintf(stderr,"value=%f\n",value);
				}
				//value *= exp(u*n);
				denominator += value;
				nominator += value*u;
				//fprintf(stderr,"deno_ = %f\n",denominator);
				//fprintf(stderr,"no_ = %f\n",nominator);
			}
			//denominator *= DELTA;
			//nominator *= DELTA;
			type1 ret = log(y_[0]) - nominator/denominator;
			//fprintf(stderr,"log(y0)=%f\n",log(y_[0]));
			//fprintf(stderr,"ret = %f\n",ret);
			return exp(ret);
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
		.property("mu_", &Receiver::get_mu)
		.property("lambda_", &Receiver::get_lambda)
		.function("ini_tau", &Receiver::ini_tau)
		.function("ini_mu", &Receiver::ini_mu)
		.function("ini_lambda", &Receiver::ini_lambda)
		.function("iter_tau", &Receiver::iter_tau)
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
	type1 xx = mu + (mu*mu*yy)/(2.0*lambda) - (mu/(2.0*lambda))*sqrt(4.0*mu*lambda*yy + mu*mu*yy*yy);
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


type1 pitman_scale(int q1){
	type1 err;
	Receiver rx_c;
	//tx
	std::vector<type1> x;
	for(int i=0; i<q1; ++i)
		x.push_back( get_IG_sample() );
	std::sort(x.begin(), x.end());

	//rx
	for(int i=0; i<q1; ++i)
		rx_c.push_y(x[i]);
	err = (rx_c.pitman_scale() - 1);
	return err;
}


type1 pitman_location(int q1){
	type1 err;
	Receiver rx_c;
	//tx
	std::vector<type1> x;
	for(int i=0; i<q1; ++i)
		x.push_back( tau_true + 0 + get_IG_sample() );
	std::sort(x.begin(), x.end());

	//rx
	for(int i=0; i<q1; ++i)
		rx_c.push_y(x[i]);
	rx_c.pitman();
	err = (rx_c.get_tau() - tau_true);
	return err;
}


void run_pitman(int test_max, int q1){
	std::vector<type1> bias(1,0);
	std::vector<type1> mse_table(1,0);
	for(int test=0;test<test_max;test++){
		//fprintf(stderr,"test=%d\n",test);

		//type1 err = pitman_location(q1);
		type1 err = pitman_scale(q1);

		bias[0] += err;
		mse_table[0] += err*err;
		//fprintf(stderr,"%d,\t%f\n", q1,  mse_table[0]/(test+1));
	}
	fprintf(stderr,"pitman_bias,\t%f\n", bias[0]/test_max);
	fprintf(stderr,"q1,\tpitman_mse,\t\n");
	printf("%d,\t%f\n", q1,  mse_table[0]/(test_max));
}




int main(int argc, char* argv[]){
	//fprintf(stderr,"hello world!\n");

	if(argc<4){
		fprintf(stderr,"Usage: %s <test_max> <K> <n1>\n", argv[0]);
		return 0;
	}

	//IG_cdf_table();

	type1 tau = tau_true, Ts=Ts_true;
	int test_max = atoi(argv[1]);
	int K = atoi(argv[2]), q1 = atoi(argv[3]);
	//type1 mean_y0 = offline(q1);
	// 2_ary
	std::uniform_int_distribution<int> uni_int(0,1);
	//int m_ary[] = {q1/2,q1*3/2};
	int m_ary[] = {q1,q1};
	// 8_ary
	//std::uniform_int_distribution<int> uni_int(0,7);
	//int m_ary[] = {q1/8,q1*3/8,q1*5/8,q1*7/8,q1*9/8,q1*11/8,q1*13/8,q1*15/8};

	run_pitman(test_max,q1);
/*
	type1 err;
	std::vector<type1> bias(3,0);
	std::vector<type1> mse_table(3,0);
	for(int test=0;test<test_max;test++){
		//fprintf(stderr,"test=%d\n",test);
		Receiver rx_c;
		//tx
		std::vector<type1> x;
		for(int i=0; i<q1; ++i)
			x.push_back( tau + 0 + get_IG_sample());
			//x.push_back( tau + 0);
		rx_c.increase_k();
		for(int k=1; k<K; ++k){
			rx_c.increase_k();
			int q = m_ary[uni_int(rd)];
			for(int i=0; i<q; ++i)
				x.push_back( tau + k*Ts + get_IG_sample());
		}
		// transmit more one symbol
		int N = x.size();
		for(int i=0; i<q1; ++i)
			x.push_back( tau + K*Ts + get_IG_sample());
		// transmit more one symbol
		std::sort(x.begin(), x.end());

		//rx
		for(int i=0; i<q1; ++i)
			rx_c.push_y(x[i]);
		rx_c.ini_tau();
		rx_c.ini_mu();
		rx_c.ini_lambda();
		//rx_c.LE_1(mean_y0);
		err = (rx_c.get_tau() - tau);
		bias[0] += err;
		mse_table[0] += err*err;
		for(int i=q1; i<N; ++i){
			rx_c.push_y(x[i]);
		}
		rx_c.est_Ts(q1);
		rx_c.iter_tau();
		err = (rx_c.get_tau() - tau);
		bias[1] += err;
		mse_table[1] += err*err;
		//err = 0;

		err += (K-1)*(rx_c.get_Ts()-Ts);
		bias[2] += err;
		mse_table[2] += err*err;
	}
	fprintf(stderr, "q1=%d\n",q1);
	fprintf(stderr,"bias_ini_tau,\t%f\n", bias[0]/test_max);
	fprintf(stderr,"mse_ini_tau,\t%f\n", mse_table[0]/test_max);
	fprintf(stderr, "K\tmse_iter_tau\tbias_iter_tau\tmse_s_k\tbias_s_k\n");

	if(K==1){
		printf("1,\t%f,\t\n",mse_table[0]/test_max);
	}
	else{
		printf("%d,\t",K);
		printf("%f,\t", mse_table[1]/test_max);
		printf("%f,\t", bias[1]/test_max);
		printf("%f,\t", mse_table[2]/test_max);
		printf("%f,\t", bias[2]/test_max);
		printf("\n");
	}
*/
	return 0;
}
#endif

