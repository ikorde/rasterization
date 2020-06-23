#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color= new pixel[width*height];
    state.image_depth= new float[width*height];
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    for(size_t i=0; i<width*height; ++i) {
    	state.image_color[i] = make_pixel(0,0,0);
    	state.image_depth[i] = 1.0; 
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    std::cout<<"TODO: implement rendering."<<std::endl;
    switch(type) {
    	case render_type::triangle : 
    	{
    		size_t j=0;
    		// std::cout << " triangle info: " << state.num_vertices << ", " << state.floats_per_vertex << std::endl; 
    		for(size_t t=0; t<state.num_vertices; t+=3) {
    			data_geometry** triangle = new data_geometry*[3];
	    		for(size_t i=0; i<3; ++i) {
	    			triangle[i] = new data_geometry; 
	    			triangle[i]->data = new float[MAX_FLOATS_PER_VERTEX];
	    			data_vertex vertex; 
	    			vertex.data = new float[MAX_FLOATS_PER_VERTEX];
					j=0;
                    while(j<state.floats_per_vertex) {
                    	vertex.data[j] = state.vertex_data[j + state.floats_per_vertex*i + state.floats_per_vertex*t];
                        triangle[i]->data[j] = state.vertex_data[j + state.floats_per_vertex*i + state.floats_per_vertex*t];
                        ++j;
                	}
                    
                	state.vertex_shader((const data_vertex)vertex, *triangle[i],state.uniform_data);
                    
            	}

            	// std::cout << "    rasteriszing tri: "  <<triangle[0]->gl_Position[0] <<", "<<triangle[0]->gl_Position[1] <<", "<<triangle[0]->gl_Position[2]<<", "<<triangle[0]->gl_Position[3] <<"1: , "
                // <<triangle[1]->gl_Position[0]<< ", "<<triangle[1]->gl_Position[1] <<", "<<triangle[1]->gl_Position[2]<<", "<<triangle[1]->gl_Position[3] <<"2: , " 
                // <<triangle[2]->gl_Position[0]<< ", "<<triangle[2]->gl_Position[2] <<", "<<triangle[2]->gl_Position[2]<<", "<<triangle[2]->gl_Position[3] <<std::endl; 
                clip_triangle(state,(const data_geometry**)triangle,0);
    		}
    		break;
    	}

    	case render_type::indexed :
    	{
            data_geometry temp[3];
            const data_geometry *triangle[3];
            data_vertex vertex[3];

            for (int i = 0; i < state.num_triangles * 3; i += 3) {
                for(int j = 0; j < 3; j++) {
                    vertex[j].data = &state.vertex_data[state.floats_per_vertex *state.index_data[i+j]];
                    temp[j].data = vertex[j].data;
                    state.vertex_shader(vertex[j],temp[j],state.uniform_data);
                    triangle[j] = &temp[j];
                }
                clip_triangle(state, triangle,0);
            }
    	    break; 
    	}

    	case render_type::fan : 
    	{
            data_geometry temp[3];
            data_vertex vertex[3];
            int fan=0;
            for (int i = 0; i < state.num_vertices; i++ ){
                const data_geometry *triangle[3];
                for(int j = 0; j < 3; j++) {
                    if(j==0) fan=1;
                    else fan=i+j;
                    vertex[j].data = &state.vertex_data[state.floats_per_vertex * fan];
                    temp[j].data = vertex[j].data;
                    state.vertex_shader(vertex[j],temp[j],state.uniform_data);
                    triangle[j] = &temp[j];
                }
                clip_triangle(state, triangle,0);
            }

    	    break;
    	}

    	case render_type::strip : 
    	{
            //data_geometry temp[3];
            data_vertex vertex[3];
            for (int i = 0; i < (state.num_vertices - 2); i++) {
                const data_geometry* triangle[3];
                for (int j = 0; j < 3; j++) {
                    data_geometry* temp = new data_geometry();
                    temp->data = &state.vertex_data[ (i+j) * state.floats_per_vertex];
                    vertex[j].data = temp->data;
                    state.vertex_shader(vertex[j], *temp, state.uniform_data);
                    triangle[j] = temp;
                }
                clip_triangle(state, triangle, 0);

            }
    	    break;
    	}

    	default: {
    		break;
    	}

    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply new_triangle the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    // std::cout << "face" << face << std::endl;
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }

    //find how many vertices intersect 
    int intersect = 0;
    vec3 notIntersect = {0,0,0};


    for (size_t i=0; i<3; ++i) {
        // std::cout << "iin[i]->gl_Position[int(face/2)] " << int(face/2) << ",, " <<  in[i]->gl_Position[0] << in[i]->gl_Position[int(face/2)] << ",   " << -1*in[i]->gl_Position[3] << std::endl; 
        if (face%2==1 && in[i]->gl_Position[int(face/2)]>= -1*in[i]->gl_Position[3] ) {
            intersect++; 
            notIntersect[i]=1;
        }
        else if (face%2==0 && in[i]->gl_Position[int(face/2)]<= in[i]->gl_Position[3]) {
            intersect++;
            notIntersect[i]=1;

       }
    }


    // std::cout<< "intersect count " << intersect << std::endl;
    //case: 0 points inside/outside
    if(intersect==0) return;

    //case: 1 point inside, clip along edge
    else if(intersect==1)  {
        // std::cout << "1" << std::endl; 

        const data_geometry* a = 0; //in[0]->gl_Position;
        const data_geometry* b = 0; //in[1]->gl_Position;
        const data_geometry* c = 0; //in[2]->gl_Position;
        const data_geometry* new_triangle[3];
        data_geometry* intersect_1 = new data_geometry;
        data_geometry* intersect_2 = new data_geometry;
        vec4 p_alpha1, p_alpha2;
        float alpha1=0.0;
        float alpha2 = 0.0;
        int sign= (face%2==0) ? 1 : -1;

        intersect_1->data = new float[state.floats_per_vertex];
        intersect_2->data = new float[state.floats_per_vertex];

        for (size_t i = 0; i < 3; i++) {
            if (notIntersect[i] == 1) {
                a = in[i] ; //->gl_Position;
                b = in[i+1 % 3] ; //->gl_Position;
                c = in[i+2 % 3] ; //->gl_Position;
                // std::cout << "in:  " << i << ", " << i+1%3 <<", " << i+2%3 << std::endl; 
                break;
            }
        }
        
        alpha1 = (sign*b->gl_Position[3]-b->gl_Position[int(face/2)]) / (a->gl_Position[int(face/2)]-sign*a->gl_Position[3]+sign*b->gl_Position[3]-b->gl_Position[int(face/2)]);
        alpha2 = (sign*c->gl_Position[3]-c->gl_Position[int(face/2)]) / (a->gl_Position[int(face/2)]-sign*a->gl_Position[3]+sign*c->gl_Position[3]-c->gl_Position[int(face/2)]);

        for (size_t i=0; i<state.floats_per_vertex; ++i) {
            switch(state.interp_rules[i]) {
            case(interp_type::flat): {
                intersect_1->data[i] = a->data[i];
                intersect_2->data[i] = a->data[i];
                break;
            } 
            case(interp_type::noperspective): {
                float temp = (a->gl_Position[3]*alpha1)/ ((a->gl_Position[3]*alpha1) + (1-alpha1)*b->gl_Position[3]);
                intersect_1->data[i] = (temp*a->data[i]) +((1-temp)*b->data[i]);
                
                temp = (alpha2*a->gl_Position[3])/ ((alpha2*a->gl_Position[3]) + (1 - alpha2)*c->gl_Position[3]);
                intersect_2->data[i] = (temp*a->data[i]) + (1-temp)*c->data[i];
                break;
            }
            case(interp_type::smooth): {
                intersect_1->data[i] = (alpha1*a->data[i]) + (1-alpha1)*b->data[i];
                intersect_2->data[i] = (alpha2*a->data[i]) + (1 - alpha2)*c->data[i];
                break;
            } 
            }
        }
        for (int i = 0; i < 4; i++) {
            intersect_1->gl_Position[i] = (alpha1 * a->gl_Position[i]) + ((1 - alpha2) * b->gl_Position[i]);
            intersect_2->gl_Position[i] = (alpha2 * a->gl_Position[i]) + ((1 - alpha2) * c->gl_Position[i]);
        }

        new_triangle[0] = a;
        new_triangle[1] = intersect_1;
        new_triangle[2] = intersect_2;
        clip_triangle(state,new_triangle,face+1);

        return;
    }

    //case: 2 points inside, call clip on 2 new triangles
    else if(intersect==2) {
        // std::cout << "2" << std::endl; 
        const data_geometry* a = 0;
        const data_geometry* b = 0;
        const data_geometry* c = 0;
        const data_geometry* new_triangle[3];
        data_geometry* intersect_1 = new data_geometry;
        data_geometry* intersect_2 = new data_geometry;
        vec4 p_alpha1, p_alpha2;
        float alpha1=0.0;
        float alpha2 = 0.0;
        int sign= (face%2==0) ? 1 : -1; 

        intersect_1->data = new float[state.floats_per_vertex];
        intersect_2->data = new float[state.floats_per_vertex];

        for (int i = 0; i < 3; i++) {

            if (notIntersect[i] == 0) {
                a = in[i];
                b = in[(i + 1) % 3];
                c = in[(i + 2) % 3];
                // std::cout << "in:  " << c->gl_Position[0] <<", " << c->gl_Position[1] <<", " << c->gl_Position[2] <<std::endl; 

                break;
            }

        }
        

        alpha1 = (sign*b->gl_Position[3]-b->gl_Position[int(face/2)]) / (a->gl_Position[int(face/2)]-sign*a->gl_Position[3]+sign*b->gl_Position[3]-b->gl_Position[int(face/2)]);
        alpha2 = (sign*c->gl_Position[3]-c->gl_Position[int(face/2)]) / (a->gl_Position[int(face/2)]-sign*a->gl_Position[3]+sign*c->gl_Position[3]-c->gl_Position[int(face/2)]);
    
        for (size_t i=0; i<state.floats_per_vertex; i++) {
            switch(state.interp_rules[i]) {
            case(interp_type::flat): {
                intersect_1->data[i] = a->data[i];
                intersect_2->data[i] = a->data[i];
                break;
            } 
            case(interp_type::noperspective): {
                float temp = (a->gl_Position[3]*alpha1)/ ((a->gl_Position[3]*alpha1) + (1-alpha1)*b->gl_Position[3]);
                intersect_1->data[i] = (temp*a->data[i]) +((1-temp)*b->data[i]);
                
                temp = (alpha2*a->gl_Position[3])/ ((alpha2*a->gl_Position[3]) + (1 - alpha2)*c->gl_Position[3]);
                intersect_2->data[i] = (temp*a->data[i]) + (1-temp)*c->data[i];
                break;
            }
            case(interp_type::smooth): {
                intersect_1->data[i] = (alpha1*a->data[i]) + (1-alpha1)*b->data[i];
                intersect_2->data[i] = (alpha2*a->data[i]) + (1 - alpha2)*c->data[i];
                break;
            } 
            }
            
        }
        for (size_t i = 0; i < 4; i++) {
            intersect_1->gl_Position[i] = (alpha1 * a->gl_Position[i]) + ((1 - alpha1) * b->gl_Position[i]);
            intersect_2->gl_Position[i] = (alpha2 * a->gl_Position[i]) + ((1 - alpha2) * c->gl_Position[i]);
        }

        new_triangle[0] = intersect_1;
        new_triangle[1] = b;
        new_triangle[2] = c;
        clip_triangle(state,new_triangle,face+1);

        new_triangle[0] = c;
        new_triangle[1] = intersect_2;
        new_triangle[2] = intersect_1;
        clip_triangle(state,new_triangle,face+1);

        return;
    } 

    else if (intersect == 3) {
        clip_triangle(state,in,face+1);
        return;
    } 

}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // std::cout<<" TODO: implement rasterization"<<std::endl;
    // std::cout<<"   in: "<< in[0]->gl_Position[0] << std::endl;
    
    float ax = state.image_width /2.0 * (in[0]->gl_Position[0]/in[0]->gl_Position[3]) + state.image_width /2.0 - 0.5; 
    float ay = state.image_height/2.0 * (in[0]->gl_Position[1]/in[0]->gl_Position[3]) + state.image_height/2.0 - 0.5; //in[0]->gl_Position[1];
    float bx = state.image_width /2.0 * (in[1]->gl_Position[0]/in[1]->gl_Position[3]) + state.image_width /2.0 - 0.5; //in[1]->gl_Position[0];
    float by = state.image_height/2.0 * (in[1]->gl_Position[1]/in[1]->gl_Position[3]) + state.image_height/2.0 - 0.5; //in[1]->gl_Position[1];
    float cx = state.image_width /2.0 * (in[2]->gl_Position[0]/in[2]->gl_Position[3]) + state.image_width /2.0 - 0.5; //in[2]->gl_Position[0];
    float cy = state.image_height/2.0 * (in[2]->gl_Position[1]/in[2]->gl_Position[3]) + state.image_height/2.0 - 0.5; //in[2]->gl_Position[1];
    float abc_area = 0.5 * ((bx*cy - cx*by) - (ax*cy - cx*ay) + (ax*by-bx*ay));
    int ind=0;
    float alpha=0.0,beta=0.0,gamma=0.0,curr_area=0.0;
    // std::cout << ax << ", " << ay << ", " << bx << ", " << by << ", " << cx << ", " << cy << std::endl;  
    
    for(size_t i=0; i<=state.image_width; ++i) {
    	for(size_t j=0; j<=state.image_height; ++j) {
			curr_area = 0.5 * ((bx*cy - cx*by) - (i*cy - cx*j) + (i*by-bx*j));
			alpha = curr_area/abc_area; 
			curr_area = 0.5 * ((i*cy - cx*j) - (ax*cy - cx*ay) + (ax*j-i*ay));
			beta = curr_area/abc_area;

			curr_area = 0.5 * ((bx*j - i*by) - (ax*j - i*ay) + (ax*by-bx*ay));
			gamma = curr_area/abc_area;
			ind = state.image_width * j + i; 

			float z_buff = alpha*in[0]->gl_Position[2]/in[0]->gl_Position[3] + beta*in[1]->gl_Position[2]/in[1]->gl_Position[3] + gamma*in[2]->gl_Position[2]/in[2]->gl_Position[3];
			if(alpha>=-0.001 && beta>=-0.001 && gamma>=-0.001 && state.image_depth[ind] > z_buff) {
				data_fragment* frags = new data_fragment(); 
    			data_output* out = new data_output();
    			float* interp = new float[MAX_FLOATS_PER_VERTEX];

				for(size_t k=0; k<state.floats_per_vertex; k++) {
					switch(state.interp_rules[k]) {
					case(interp_type::flat): {
						interp[k]= in[0]->data[k];
						break;
					} 
					case(interp_type::noperspective): {
						interp[k]= alpha*in[0]->data[k]+beta*in[1]->data[k]+gamma*in[2]->data[k];
						break;
					}
                    case(interp_type::smooth): {
                        float denom = alpha/in[0]->gl_Position[3] + beta/in[1]->gl_Position[3] + gamma/in[2]->gl_Position[3];
                        interp[k]=((alpha*in[0]->data[k] / in[0]->gl_Position[3]) + (beta*in[1]->data[k] / in[1]->gl_Position[3]) + (gamma*in[2]->data[k]/in[2]->gl_Position[3]) )/ denom;

                    }
					}
				}
				frags->data = interp;
				state.fragment_shader((const data_fragment)*frags, *out, state.uniform_data);
				state.image_color[ind] = make_pixel(out->output_color[0]*255, out->output_color[1]*255, out->output_color[2]*255);
				state.image_depth[ind] = z_buff;
			}	
    	}
    }
    // std::cout << "   rast done " << std::endl;
}

