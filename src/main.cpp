#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include <algorithm>
#include "spline.h"
#include <map>
#include <float.h>

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

/*
 *      Finite state machine class
 */

class FSM {

    string state;
    vector<string> states;
    map<string, vector<string>> next_states;

    int lane;
//    int final_lane;
    int intended_lane;

    const double v_max = 49.5;
    
    vector<double> get_car(vector<vector<double>> sensor_fusion, double car_s, int from_lane, bool from_front, int prev_size);

public:
    FSM(string init_state, int init_lane);

    vector<vector<double>> get_wp(vector<vector<double>> sensor_fusion, double car_s, double &v_target, int prev_size, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y);

    void update_state(vector<vector<double>> sensor_fusion, double car_s, double car_d, int prev_size);
};

/*
 *      FSM implementation
 */

FSM::FSM(string init_state, int init_lane) {
    states = {"kl", "plcl", "plcr", "lcl" "lcr"};

    next_states["kl"] = {"kl", "plcl", "plcr"};
    next_states["plcl"] = {"plcl", "kl", "lcl"};
    next_states["plcr"] = {"plcr", "kl", "lcr"};
    next_states["lcl"] = {"lcl", "kl"};
    next_states["lcr"] = {"lcr", "kl"};

    vector<string>::iterator it = find(states.begin(), states.end(), init_state);

    if (it != states.end()) {
        state = init_state;
    } else {
        state = "kl";
    }

    lane = init_lane;
}


vector<double> FSM::get_car(vector<vector<double>> sensor_fusion, double car_s, int from_lane, bool from_front, int prev_size) {

    int car_i = -1; 
    double dist = DBL_MAX;
    if (!from_front) dist = DBL_MIN;

    double speed = 0;

    for (int i = 0; i < sensor_fusion.size(); ++i) {
        double d = sensor_fusion[i][6];

        if (d < (2 + 4 * from_lane + 2) && d > (2 + 4 * from_lane - 2)) {

            double vx = sensor_fusion[i][3];
            double vy = sensor_fusion[i][4];
            double check_speed = sqrt(vx*vx + vy*vy);
            double check_car_s = sensor_fusion[i][5];

            check_car_s += (double)prev_size * 0.02 * check_speed;
            double test_dist = check_car_s - car_s;

            if (from_front) {
                if (test_dist > 0 && test_dist < dist) {
                    car_i = i;
                    dist = test_dist;
                    speed = check_speed;
                }
            } else {
                if (test_dist < 0 && test_dist > dist) {
                    car_i = i;
                    dist = test_dist;
                    speed = check_speed;
                }
            }
        }
    }

    vector<double> car;

    if (car_i != -1) {
        car = sensor_fusion[car_i];
        car.push_back(dist);
        car.push_back(speed);
    }

    return car;

}

vector<vector<double>> FSM::get_wp(vector<vector<double>> sensor_fusion, double car_s, double &v_target, int prev_size, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y) {

    vector<vector<double>> wps;
    int i;

    double d1 = 25;
    double d2 = 50;

    if (state == "kl") {
        for (i = 0; i < 3; ++i) {
            vector<double> wp = getXY(car_s + 30 * (i + 1), 2 + 4 * lane, maps_s, maps_x, maps_y);
            wps.push_back(wp);
        }
        
        v_target = v_max;

        vector<double> front_car = get_car(sensor_fusion, car_s, lane, true, prev_size);

        if (front_car.size() > 0){
            double dist = front_car[7];
            double speed = front_car[8];
            if (dist < d1) v_target = 2.24 * speed * 0.75;
            else if (dist < d2) v_target = 2.24 * speed;
        } 

    } else if (state == "plcl" || state == "plcr") {        
        v_target = v_max;

        vector<double> front_car = get_car(sensor_fusion, car_s, lane, true, prev_size);

        if (front_car.size() > 0){
            double dist = front_car[7];
            double speed = front_car[8];
            if (dist < d1) v_target = 2.24 * speed * 0.75;
            else if (dist < d2) v_target = 2.24 * speed;
        }

        vector<double> front_side_car = get_car(sensor_fusion, car_s, intended_lane, true, prev_size);

        if (front_side_car.size() > 0) {
            double dist = front_side_car[7];
            double speed = front_side_car[8];
            if (dist < d1) v_target = min(v_target, 2.24 * speed * 0.75);
            else if (dist < d2) v_target = min(v_target, 2.24 * speed);
        }

        for (i = 0; i < 3; ++i) {
            vector<double> wp = getXY(car_s + 30 * 3 * v_target / v_max * (i + 1), 2 + 4 * lane, maps_s, maps_x, maps_y);
            wps.push_back(wp);
        }

    } else if (state == "lcl" || state == "lcr") {
        for (i = 0; i < 3; ++i) {
            vector<double> wp = getXY(car_s + 30 * (i + 1), 2 + 4 * intended_lane, maps_s, maps_x, maps_y);
            wps.push_back(wp);
        }
        
        v_target = v_max;

        vector<double> front_car = get_car(sensor_fusion, car_s, lane, true, prev_size);

        if (front_car.size() > 0){
            double dist = front_car[7];
            double speed = front_car[8];
            if (dist < d1) v_target = 2.24 * speed * 0.75;
            else if (dist < d2) v_target = 2.24 * speed;
        }

        if (intended_lane != lane) { 
            vector<double> front_side_car = get_car(sensor_fusion, car_s, intended_lane, true, prev_size);

            if (front_side_car.size() > 0) {
                double dist = front_side_car[7];
                double speed = front_side_car[8];
                if (dist < d1) v_target = min(v_target, 2.24 * speed * 0.75);
                else if (dist < d2) v_target = min(v_target, 2.24 * speed);
            }
        }

    } 

    return wps;
}

void FSM::update_state(vector<vector<double>> sensor_fusion, double car_s, double car_d, int prev_size) {
    vector<string> possible_states = next_states[state];

    double speed_weight = 0.5;
    double safety_weight = 70.0;
    int optimal_i = 0;
    double cost = DBL_MAX;

    double ds1 = 70;
    double ds2 = 50;

    for (int i = 0; i < possible_states.size(); ++i) {
        double state_cost = DBL_MAX;
        string possible_state = possible_states[i];
        
        if (possible_state == "kl") {

            vector<double> front_car = get_car(sensor_fusion, car_s, lane, true, prev_size);
            if (front_car.size() > 0) {
                double dist = front_car[7];
                double speed = front_car[8];
                if (dist < 60) {
                    state_cost = speed_weight * (v_max - speed);
                } else {
                    state_cost = 0;
                }
            } else {
                state_cost = 0;
            }

            if ((state == "lcl" && car_d < 2 + 4 * lane - 1) || (state == "lcr" && car_d > 2 + 4 * lane + 1)) {
                state_cost = DBL_MAX;
            }

        } else if (possible_state == "plcl" && lane > 0) {

            vector<double> front_left_car = get_car(sensor_fusion, car_s, lane - 1, true, prev_size);

            vector<double> rear_left_car = get_car(sensor_fusion, car_s, lane - 1, false, prev_size);

            state_cost = 1;

            if (front_left_car.size() > 0) {
                double dist = front_left_car[7];
                double speed = front_left_car[8];
                if (dist < ds1) {
                    state_cost += speed_weight * (v_max - speed);
                }
                if (dist < ds2) {
                    state_cost += safety_weight / max(dist, 0.01);
                }
            }

            if (rear_left_car.size() > 0) {
                double dist = rear_left_car[7];
                double speed = rear_left_car[8];
                if (abs(dist) < ds2) {
                    state_cost += safety_weight / (2 * max(abs(dist), 0.01));
                }
            }

        } else if (possible_state == "plcr" && lane < 2) {

            vector<double> front_right_car = get_car(sensor_fusion, car_s, lane + 1, true, prev_size);

            vector<double> rear_right_car = get_car(sensor_fusion, car_s, lane + 1, false, prev_size);

            state_cost = 1;

            if (front_right_car.size() > 0) {
                double dist = front_right_car[7];
                double speed = front_right_car[8];
                if (dist < ds1) {
                    state_cost += speed_weight * (v_max - speed);
                }
                if (dist < ds2) {
                    state_cost += safety_weight / max(dist, 0.01);
                }
            }

            if (rear_right_car.size() > 0) {
                double dist = rear_right_car[7];
                double speed = rear_right_car[8];
                if (abs(dist) < ds2) {
                    state_cost += safety_weight / (2 * max(abs(dist), 0.01));
                }
            }

        } else if (possible_state == "lcl" || possible_state == "lcr") {

            vector<double> front_side_car = get_car(sensor_fusion, car_s, intended_lane, true, prev_size);

            vector<double> rear_side_car = get_car(sensor_fusion, car_s, intended_lane, false, prev_size);

            state_cost = 0.2;

            if (front_side_car.size() > 0) {
                double dist = front_side_car[7];
                double speed = front_side_car[8];
                if (dist < ds1) {
                    state_cost += speed_weight * (v_max - speed);
                }
                if (dist < ds2) {
                    state_cost += safety_weight / max(dist, 0.01);
                }
            }

            if (rear_side_car.size() > 0) {
                double dist = rear_side_car[7];
                double speed = rear_side_car[8];
                if (abs(dist) < ds2) {
                    state_cost += safety_weight / max(abs(dist), 0.01);
                }
            }

        }

        if (state_cost < cost) {
            optimal_i = i;
            cost = state_cost;
        }
    }

    string old_state = state;
    state = possible_states[optimal_i];

    if (old_state != state) cout << "Changed to state: " << state << endl;

    if (state == "plcl") {
        intended_lane = lane - 1;
    } else if (state == "plcr") {
        intended_lane = lane + 1;
    } else if ((state == "lcl" && car_d < 4 * lane) || (state == "lcr" && car_d > 4 * lane + 4)) {
        lane = intended_lane;
    }

}

/*
 *      End of FSM implementation
 */


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;

  double v_ref = 0;
  const double v_max = 49.5;

  FSM fsm ("kl", lane);

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &v_ref, &v_max, &fsm](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
			double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.

          	auto sensor_fusion = j[1]["sensor_fusion"];

			int prev_size = min(previous_path_x.size(), previous_path_y.size());

          	double v_target;

			json msgJson;

			vector<double> xs;
			vector<double> ys;

			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);

			if (prev_size < 2) {
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);

				xs.push_back(prev_car_x);
				xs.push_back(car_x);

				ys.push_back(prev_car_y);
				ys.push_back(car_y);
			} else {
				ref_x = previous_path_x[prev_size - 1];
				ref_y = previous_path_y[prev_size - 1];

				double prev_ref_x = previous_path_x[prev_size - 2];
				double prev_ref_y = previous_path_y[prev_size - 2];
				ref_yaw = atan2(ref_y - prev_ref_y, ref_x - prev_ref_x);

				xs.push_back(prev_ref_x);
				xs.push_back(ref_x);

				ys.push_back(prev_ref_y);
				ys.push_back(ref_y);
			}
			
			int i;

            fsm.update_state(sensor_fusion, car_s, car_d, prev_size);

            vector<vector<double>> wps = fsm.get_wp(sensor_fusion, car_s, v_target, prev_size, map_waypoints_s, map_waypoints_x, map_waypoints_y);

			for (i = 0; i < 3; i++) {
				vector<double> wp = wps[i];

				xs.push_back(wp[0]);
				ys.push_back(wp[1]);
			}

			for (i = 0; i < xs.size(); i++) {
				double shift_x = xs[i] - ref_x;
				double shift_y = ys[i] - ref_y;

				xs[i] = shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw);
				ys[i] = shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw);

			}

			tk::spline s;

			s.set_points(xs, ys);

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds


			for (i = 0; i < prev_size; ++i) {
		  		next_x_vals.push_back(previous_path_x[i]);
		  		next_y_vals.push_back(previous_path_y[i]);
			}

			double target_x = 30.0;
			double target_y = s(target_x);
			double target_dist = sqrt(target_x * target_x + target_y * target_y);

            double dv_max = 0.4; //mph
            double dv_min = 0.05;
            double v_err = 0.5;
			
            double x_add = 0.0;

   			for(i = 0; i < 50 - prev_size; i++) {
                double dv = min(dv_max, (v_target - v_ref)/5);
                dv = max(dv, dv_min);
                if (v_ref < v_target - v_err) {
                    v_ref += dv;
                } else if (v_ref > v_target) {
                    v_ref -= dv;
                }
                //cout << "Target: " << v_target << ", ref: " << v_ref << endl;
			    double N = target_dist / (0.02 * v_ref / 2.24);

				double x_point = x_add + target_x / N;
				double y_point = s(x_point);
				x_add = x_point;
				
				double xp = x_point;
				double yp = y_point;

				x_point = xp * cos(ref_yaw) - yp * sin(ref_yaw) + ref_x;
				y_point = xp * sin(ref_yaw) + yp * cos(ref_yaw) + ref_y;
//				cout << "i: " << i + prev_size << ", x: " << x_point << ", y: " << y_point << endl;
				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
    		}

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
