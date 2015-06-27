

function printf(str, value){
	console.log(str);
	console.log(value);
}


var Receiver = {
	createNew: function(tx){
		var rx = {};
		rx.rx_c;
		rx.init = function(){
			rx.rx_c = new Module.Receiver();
			for(var i = 0; i<tx.K; ++i){
				var id = "rx_time_"+parseInt(i);
				$("#simu").append("<div id='"+id+"'/>");
				$("#"+id).attr("class","rx_time");
			}
		}
		rx.catchMolecules = function(){
			for(var i=0; i<tx.K; ++i){
				tx.symbols[i].catchMolecules(rx);
			}
		}
		rx.move = function(){
			for(var k = 0; k<tx.K; ++k){
				var id = "rx_time_"+parseInt(k);
				tx.s1.move_sth(id);
			}
		}
		rx.update = function(){
			var max_k = rx.rx_c.k_;
			for(var k = 0; k<Math.min(tx.K, max_k); ++k){
				var id = "rx_time_"+parseInt(k);
				var est = rx.rx_c.tau_ + (max_k-1-k)*rx.rx_c.Ts_;
				$("#"+id).position({my:"center", at:"center-"+parseInt(tx.s1.get_count()-est), of:"#space"});
			}
		}
		rx.init_est = function(){
			/*
			var n_1 = tx.n_i[0];
			var y_bar = rx.get_mean(0,n_1);
			var m_2 = rx.get_m2(y_bar, n_1);
			var y_1 = rx.observ[0];
			//rx.tau_hat = y_1 - Math.pow(y_bar-y_1, 3)/2/m_2/Math.log(n_1);
			rx.tau_hat = y_1 - Math.pow(y_bar-y_1, 3)/m_2/Math.log(n_1);
			//here is some problem!!!!!!!!!!!! why need divid two
			//rx.tau_hat = rx.observ[0] - 200;
			rx.m_y = y_bar;

			//mu init est
			rx.mu = y_bar - rx.tau_hat;
			printf("mu_init",rx.mu);
			//lambda init est
			var m_n1 = rx.get_n1_moment(0,n_1);
			printf("neg_1_moment",m_n1)
			rx.lambda = 1/(m_n1-1/rx.mu);
			printf("lambda_init",rx.lambda);

			for(var i=0; i<rx.observ.length; ++i){
				//console.log(rx.observ[i]);
			}
			*/
			//console.log(y_bar);
			//console.log(m_2);
			//console.log(n_1);
			//console.log(y_1);
			printf("tau_init",rx.rx_c.tau_);
			rx.update();
		}
		rx.iter_est = function(){
			/*
			var k = rx.rx_c.get_k();
			var y_bar = rx.get_mean((k-1)*tx.n_i[0], tx.n_i[0]);
			rx.m_y = rx.m_y*(k-1)/k + y_bar/k;
			rx.m_xy = rx.m_xy + (k-1)*y_bar;
			rx.Ts_hat = (rx.m_xy*2/(k-1)/k - rx.m_y)*6/(k+1);
			*/

			//console.log(rx.Ts_hat);
			printf("tau_init",rx.rx_c.tau_);
			printf("Ts",rx.rx_c.Ts_);
			rx.update();
		}
		rx.caught = function(time){
			rx.rx_c.push_y(time);
			n = rx.rx_c.get_n();
			console.log(rx.rx_c.get_y(n-1));
			if(n == tx.n_i[0]){
				rx.rx_c.increase_k();
				rx.init_est();
			}
			if(n>tx.n_i[0] && n%tx.n_i[0]==0){
				rx.rx_c.increase_k();
				rx.iter_est();
			}
		}
		rx.init();
		return rx;
	}
}


