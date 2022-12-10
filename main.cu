//This is mainly based on the sequence given by https://conference.sdo.esoc.esa.int/proceedings/sdc7/paper/14/SDC7-paper14.pdf
//The timestep calculation is based on Gravitational N-Body Simulations, SVERRE J. AARSETH

#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <fstream>
#include <exception>
#include <sstream>
#include <thrust/copy.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


using namespace std;

__constant__ double G0=6.6743015e-11, PI0=3.1415926535898;
const double G=6.6743015e-11, PI=3.1415926535898;
thrust::host_vector <double> add_init_position, add_init_velocity, body_values(7, 0.), output_vec, r_rel, v_rel, acting_values(7, 0.), body_final,
					r_f,v_f,accel_output,a_0,a_0_dot,a_p,a_p_dot,
					a_i,a_i_dot,a_j,a_j_dot,r_i,v_i,r_j,v_j,a_ij,a_dot_ij,a_t_dot_ij,a_d_dot_ij,a_d_dot,a_t_dot,
					burn_vector, burn_vector_next, burn_ori, burn_ori_rate;
thrust::host_vector <string> load_results;
double add_mass, r_dot_v_relative, abs_r_rel, a_comp, a_dot_comp_1, a_dot_comp_2, abs_v_rel, comp_a, comp_b, comp_c, body_timestep, abs_a_t_dot, abs_a_d_dot, abs_a_dot, abs_a;
string line, add_id, body_name;
double body_id, acting_id, burn_body, burn_id;
int num_bodies = 0;
int i, itts, time_counter;
int burn_count = 0;
// thrust::pair<string, thrust::host_vector<double> > burn_values, burn_values_next;
thrust::host_vector<double> burn_values, burn_values_next;

map<string,thrust::host_vector<double> > bodies, bodies_next;
map<string, thrust::pair<string, thrust::host_vector<double> > > burns;

thrust::host_vector<string> body_names;

thrust::host_vector<double> burn_data;
thrust::host_vector<double> body_data;
thrust::host_vector<double> body_data_next;

// thrust::device_vector<double> body_ddata;

double accuracy, timestep, time0, next_timestep;
int output_rate;
string info, debug_info;
void SetAccuracy(int new_accuracy){
	accuracy=new_accuracy;
}
void SetStartTime(int start_time){
	time0=start_time;
}
void SetOutToFile(const char * file_name, int rate){
	freopen(file_name,"w",stdout);
	output_rate = rate;
}
void AddBody (string id, double mass, thrust::host_vector<double> init_position, thrust::host_vector<double> init_velocity){
    bodies[id].clear();//so adding bodies with the same id twice doesn't break it
	bodies[id].push_back(mass);

	for(i=0;i<init_position.size();i++){
		bodies[id].push_back(init_position[i]);
	}

	for(i=0;i<init_velocity.size();i++){
		bodies[id].push_back(init_velocity[i]);
	}
}

void AddBody2(string id, double mass, thrust::host_vector<double> init_position, thrust::host_vector<double> init_velocity){
	body_names.push_back(id);
	body_data.push_back(mass);

	for(i=0;i<init_position.size();i++){
		body_data.push_back(init_position[i]);
	}

	for(i=0;i<init_velocity.size();i++){
		body_data.push_back(init_velocity[i]);
	}
}

void AddBurn (double burn_id, double body_id, double start_time, double end_time, double acceleration, thrust::host_vector<double> orientation, thrust::host_vector<double> orientation_rate){
	//Gives burn_id:(body_id:[start,end,accel,orie_x,orie_y,orie_z,orie_ra_x...])
	// burn_vector.clear();
	burn_data.push_back(burn_id);
	burn_data.push_back(body_id);
	burn_data.push_back(start_time);
	burn_data.push_back(end_time);
	burn_data.push_back(acceleration);
	for(i=0;i<3;i++){
		// burn_vector.push_back(orientation[i]);
		burn_data.push_back(orientation[i]);
	}
	for(i=0;i<3;i++){
		// burn_vector.push_back(orientation_rate[i]);
		burn_data.push_back(orientation_rate[i]);
	}
	// burns[burn_id] = thrust::make_pair(body_id,burn_vector);
}


void LoadFile(string filename){
	//Format is #body,id,mass,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z
	body_data.clear();
	ifstream file(filename);
	if (file.is_open()) {
		string line;
		int j = 0;
		while (getline(file, line)) {
			if(line.compare(0,1,"#")==0){
				if(line.compare(1, 4,"Body")==0){
					load_results.clear();
					stringstream s_stream(line);
					while(s_stream.good()) {
					    string substr;
					    getline(s_stream, substr, ',');
					    load_results.push_back(substr);
					}
					add_init_position.clear();
					add_init_velocity.clear();

					add_id = load_results[1];
					add_mass = stold(load_results[2]);

					for(i=0;i<3;i++){
						add_init_position.push_back(stold(load_results[3+i]));
						add_init_velocity.push_back(stold(load_results[6+i]));
					}
					// AddBody(add_id, add_mass, add_init_position, add_init_velocity);
					AddBody2(add_id, add_mass, add_init_position, add_init_velocity);
					num_bodies++;
				}
				if(line.compare(1, 4,"Burn")==0){
					load_results.clear();
					stringstream s_stream(line);
					while(s_stream.good()) {
					    string substr;
					    getline(s_stream, substr, ',');
					    load_results.push_back(substr);
					}
					burn_ori.clear();burn_ori_rate.clear();
					for(i=0;i<3;i++){
						burn_ori.push_back(stold(load_results[5+i]));
						burn_ori_rate.push_back(stold(load_results[8+i]));
					}
					// burn_id = load_results[1];
					cout<<load_results[1];
					AddBurn(burn_count,num_bodies-1,stold(load_results[2]),stold(load_results[3]),stold(load_results[4]),burn_ori,burn_ori_rate);
					burn_count++;
				}
				if(line.compare(1, 8,"Addition")==0){
					cout<<"This functionality has not yet been added";
					throw exception();
				}
				++j;
			}
		}
		body_data_next = body_data;
		file.close();
	}
}


#if __CUDA_ARCH__ < 600
__device__ double atomicAdd1(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif


__device__ void CalcAAndDotD(double* output_vec, double body_id, double* r, double* v, double* body_data, 
	double* burn_data, size_t bodySize, size_t burnSize, double time0){
	// vector<double> a, a_dot, acting_values, r_rel, v_rel;
	double acting_values[7], r_rel[3], v_rel[3], a[3], a_dot[3];
	double acting_id, r_dot_v_relative, abs_r_rel, a_comp, a_dot_comp_1, a_dot_comp_2;
	int i;

	for(i = 0;i < 3;i++){
		a[i] = 0;
		a_dot[i] = 0;
	}
	acting_id = threadIdx.x;
	int q = 7 * acting_id;
	if(acting_id != body_id){
		// acting_values = body_itterator.second;
		for (i = 0; i < 7; i++)
			acting_values[i] = body_data[q+i];
		for(i=0;i<3;i++){
			r_rel[i] = r[i] - acting_values[1+i];
		}
		for(i=0;i<3;i++){
			v_rel[i] = v[i] - acting_values[4+i];
		}
		for(i=0;i<3;i++){
			r_dot_v_relative = r_rel[i] * v_rel[i];
		}

		abs_r_rel = sqrt(pow(r_rel[0],2)+pow(r_rel[1],2)+pow(r_rel[2],2));

		a_comp = -G0*acting_values[0]/pow(abs_r_rel,3);
		a_dot_comp_1 = 3*G0*acting_values[0]*r_dot_v_relative/pow(abs_r_rel,5);
		a_dot_comp_2 = -G0*acting_values[0]/pow(abs_r_rel,3);

		for(i=0;i<3;i++){
			a[i] += a_comp*r_rel[i];
			a_dot[i] += a_dot_comp_1*r_rel[i]+a_dot_comp_2*v_rel[i];
		}
	}

	for(i=0;i<3;i++){
		// output_vec[i] = a[i];
		if (a[i])
			atomicAdd1(&output_vec[i], a[i]);
	}
	for(i=0;i<3;i++){
		// output_vec[i+3] = a_dot[i];
		if (a_dot[i])
			atomicAdd1(&output_vec[i + 3], a_dot[i]);
	}
	__syncthreads();
	if (threadIdx.x)
		return;
	// master thread shall work on burn, if any
	double burn_body;
	double burn_values[11], burn_vector[10];

	for (int q = 0; q < burnSize; q += 11){
		for (i = 0; i < 11; i++){
			// burn_values.push_back(burn_data[q + i]);
			burn_values[i] = burn_data[q + i];
		}
		// burn_values = burn_itt.second;
		burn_body = burn_values[1];
		if(burn_body == body_id){
			for (i = 1; i < 11; i++){
				// burn_vector.push_back(burn_values[i]);
				burn_vector[i - 1] = burn_values[i];
			}
			// burn_vector = burn_values.second;
			if(time0>burn_vector[0]&&time0<burn_vector[1]){
				for(i=0;i<3;i++){
					output_vec[i] = burn_vector[2]*burn_vector[3+i];
				}
			}
		}
	}
}


__global__ void perform(double* body_data0, double* burn_data, size_t bodySize, size_t burnSize, double time0, double timestep, 
	double* body_data_next) {
	int i, q;
	// create shared array for every block instead of using global data
	__shared__ double body_data[100];
	for (i =0; i < 7; i++)
		body_data[threadIdx.x * 7 + i] = body_data0[threadIdx.x * 7 + i];
	__syncthreads();

	double r_0[3], v_0[3];
	double body_values[7], body_id;
	q = blockIdx.x * 7;
	for (i = 0; i < 7; i++)
		body_values[i] = body_data[q+i];
	body_id = q / 7;

	for(i=0;i<3;i++){
		r_0[i] = body_values[1+i];
		v_0[i] = body_values[4+i];
	}

	__shared__ double accel_output[6];
	if (!threadIdx.x)
		for (i = 0; i < 6; i++)
			accel_output[i] = 0;
	__syncthreads();
	CalcAAndDotD(accel_output, body_id, r_0, v_0, body_data, burn_data, bodySize, burnSize, time0);
	// barrier not really needed here but whatever
	__syncthreads();
	__shared__ double r_p[3], v_p[3];
	double a_0[3], a_0_dot[3];
	if (!threadIdx.x){
		for(i=0;i<3;i++){
			a_0[i] = accel_output[i];
			a_0_dot[i] = accel_output[3+i];
		}
		
		for(i=0;i<3;i++){
			r_p[i] = r_0[i]+v_0[i]*timestep+0.5*a_0[i]*pow(timestep,2)+(1/6)*a_0_dot[i]*pow(timestep,3);
			v_p[i] = v_0[i]+a_0[i]*timestep+0.5*a_0_dot[i]*pow(timestep,2);
		}
	}	

	double a_p[3], a_p_dot[3];
	for(int itts=0;itts<2;itts++){
		__syncthreads();
		CalcAAndDotD(accel_output, body_id, r_p, v_p, body_data, burn_data, bodySize, burnSize, time0);
		// barrier not really needed here but whatever
		__syncthreads();
		if (!threadIdx.x){
			for(i=0;i<3;i++){
				a_p[i] = accel_output[i];
				a_p_dot[i] = accel_output[3+i];
			}

			for(i=0;i<3;i++){
				v_p[i] = v_0[i]+0.5*(a_0[i]+a_p[i])*timestep+(1/12)*(a_0_dot[i]-a_p_dot[i]*pow(timestep,2));
				r_p[i] = r_0[i]+0.5*(v_p[i]+v_0[i])*timestep+(1/12)*(a_0[i]-a_p[i])*pow(timestep,2);
			}
		}
	}
	// only master thread is premitted to write to global data
	if (threadIdx.x) return;
	double body_final[7];
	body_final[0] = body_values[0];

	for(i=0;i<3;i++){
		body_final[i + 1] = r_p[i];
	}
	for(i=0;i<3;i++){
		body_final[i + 4] = v_p[i];
	}
	// begin
	for (i = 0; i < 7; i++){
		body_data_next[q+i] = body_final[i];
	}
}


thrust::device_vector<double> body_ddata, burn_ddata, body_ddata_next;
double* body_ddata_ptr, *burn_ddata_ptr, *body_ddata_next_ptr;
void initialize(){
	body_ddata = body_data;
	body_ddata_ptr = thrust::raw_pointer_cast(body_ddata.data());
	burn_ddata = burn_data;
	burn_ddata_ptr = thrust::raw_pointer_cast(burn_ddata.data());
	body_ddata_next = body_data;
	body_ddata_next_ptr = thrust::raw_pointer_cast(body_ddata_next.data());
}


void swapPtr(double **r, double **s)
{
    double *pSwap = *r;
    *r = *s;
    *s = pSwap;
}


void Step() {
	next_timestep = 9999999999;
	// begin
	perform<<<num_bodies, num_bodies>>>(body_ddata_ptr, burn_ddata_ptr, body_data.size(), burn_data.size(), time0, timestep, 
		body_ddata_next_ptr);
	// body_ddata = body_ddata_next;
	swapPtr(&body_ddata_ptr, &body_ddata_next_ptr);
	time0 = time0 + timestep;
}


void Output(){
	if(time_counter==output_rate){
		cout<<"#"+to_string(time0)+"\n";
		// begin
		for (int q = 0; q < body_data.size(); q+=7){
			for (i = 0; i < 7; i++)
				// copy from gpu to cpu
				body_values[i] = body_ddata[q+i];
			body_id = q / 7;
			body_name = body_names[body_id];
			cout<<body_name;
			for(i = 0; i < body_values.size(); i++){
				cout << "," + to_string(body_values[i]);
			}
			cout<<"\n";
		}
		time_counter=0;
	}
	else{
		time_counter++;
	}
}
