
var Receiver = {
	createNew: function(tx){
		var rx = {};
		rx.tau_hat = 0;
		rx.Ts_hat = 0;
		rx.observ = [];
		rx.sum_n_i = [];

		rx.N_ = 0;
		rx.m_y = 0;
		rx.m_xy = 0;



		rx.init = function(){
			for(var i = 0; i<tx.K; ++i){
				var id = "rx_time_"+parseInt(i);
				$("#simu").append("<div id='"+id+"'/>");
				$("#"+id).attr("class","rx_time");
			}
			var ss = 0;
			for(var i=0; i<tx.K; ++i){
				ss = ss + tx.n_i[i];
				rx.sum_n_i.push(ss);
			}
		}
		rx.catchMolecules = function(){
			for(var i=0; i<tx.K; ++i){
				tx.symbols[i].catchMolecules(rx);
			}
		}
		rx.update = function(max_k){
			for(var k = 0; k<max_k; ++k){
				var id = "rx_time_"+parseInt(k);
				var est = rx.tau_hat + k*rx.Ts_hat;
				//$("#"+id).position({my:"left", at:"left+"+parseInt(est), of:"#simu"});
			}
		}
		rx.get_mean = function(start, num){
			var ret = 0;
			for(var i=0; i<num; ++i)
				ret = ret + rx.observ[start + i];
			return ret/num;
		}
		rx.get_2_moment = function(start, num){
			var ret = 0;
			for(var i=0; i<num; ++i)
				ret = ret + Math.pow(rx.observ[start+i],2);
			return ret/num;
		}
		rx.init_est = function(){
			var n_1 = tx.n_i[0];
			var y_bar = rx.get_mean(0,n_1);
			var m_2 = rx.get_2_moment(0, n_1) - Math.pow(y_bar,2);
			var y_1 = rx.observ[0];
			//rx.tau_hat = y_1 - Math.pow(y_bar-y_1, 3)/2/m_2/Math.log(n_1);
			rx.tau_hat = y_1 - Math.pow(y_bar-y_1, 3)/m_2/Math.log(n_1);
			//rx.tau_hat = tx.tau;
			//rx.tau_hat = rx.observ[0] - 200;
			rx.m_y = y_bar;

			for(var i=0; i<rx.observ.length; ++i){
				console.log(rx.observ[i]);
			}
			console.log(y_bar);
			console.log(m_2);
			console.log(n_1);
			console.log(y_1);
			console.log(rx.tau_hat);
			rx.update(1);
		}
		rx.iter_est = function(i){
			var k = i+1;
			var y_bar = rx.get_mean(i*tx.n_i[0], tx.n_i[0]);
			rx.m_y = rx.m_y*i/k + y_bar/k;
			rx.m_xy = rx.m_xy + i*y_bar;
			rx.Ts_hat = (rx.m_xy*2/(k-1)/k - rx.m_y)*6/(k+1);

			console.log(rx.Ts_hat);
			rx.update(k);
		}
		rx.caught = function(time){
			rx.observ.push(time);
			/*
			for(var i=0; i<rx.observ.length; ++i){
				console.log(rx.observ[i]);
			}
			*/
			if( rx.observ.length == tx.n_i[0])
				rx.init_est();
			for(var i=1; i<tx.K; ++i){
				if(rx.observ.length == rx.sum_n_i[i])
					rx.iter_est(i);
			}
		}
		rx.init();
		return rx;
	}
}


