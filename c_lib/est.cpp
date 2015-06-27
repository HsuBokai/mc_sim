
#include <vector.h>
#include <math.h>
#include <emscripten/bind.h>

using namespace emscripten;

class Receiver {
	public:
		Receiver():tau_(50), Ts_(200), k_(0){}

		void push_y(float y) {
			y_.push_back(y);
		}
		float get_y(int index){return y_[index];}
		int get_n(){return y_.size();}


		int get_k() const {return k_;}
		void increase_k(){++k_;}

		float get_tau() const {return tau_;}
		float get_Ts() const {return Ts_;}

	private:
		float tau_;
		float Ts_;
		//float mu_;
		//float lambda_;
		int k_;
		std::vector<float> y_;
};

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
		//.class_function("getStringFromInstance", &Receiver::getStringFromInstance)
		;
}

