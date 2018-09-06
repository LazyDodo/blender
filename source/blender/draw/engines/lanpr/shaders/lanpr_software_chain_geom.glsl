layout(lines_adjacency) in;
layout(triangle_strip, max_vertices = 6) out;

in vec2 gOffset[];
in int gType[];
in int gLevel[];

in vec3 gNormal[];
uniform int normal_mode;
uniform int normal_effect_inverse;
uniform vec3 normal_direction; // also used as point position
uniform float normal_ramp_begin;
uniform float normal_ramp_end;
uniform float normal_thickness_begin;
uniform float normal_thickness_end;

uniform float thickness;
uniform float thickness_crease;
uniform float thickness_material;
uniform float thickness_edge_mark;
uniform float thickness_intersection;

uniform int enable_contour;
uniform int enable_crease;
uniform int enable_material;
uniform int enable_edge_mark;
uniform int enable_intersection;

uniform int occlusion_level_begin;
uniform int occlusion_level_end;

// implement these later.
//uniform float depth_width_influence;
//uniform float depth_width_curve;
//uniform float depth_alpha_influence;
//uniform float depth_alpha_curve;
//uniform float zNear;
//uniform float zFar;

uniform vec4 color;
uniform vec4 crease_color;
uniform vec4 material_color;
uniform vec4 edge_mark_color;
uniform vec4 intersection_color;

uniform float taper_l_dist;
uniform float taper_r_dist;
uniform float taper_l_strength;
uniform float taper_r_strength;

// for line width correction
uniform vec4 output_viewport;
uniform vec4 preview_viewport;

out vec4 out_color;

float use_thickness;

#define M_PI 3.1415926535897932384626433832795

vec4 END_POINT = vec4(vec2(3e30f), 0, 1);// end point flag

vec4 MakeLeftTaperLinear(vec4 L, vec4 a, float offset){
	if (offset >= taper_l_dist) return a;
	a = mix(mix(a, L, taper_l_strength), a, offset / taper_l_dist);
	return a;
}

vec4 MakeRightTaperLinear(vec4 R, vec4 c, float offset){
	if (offset >= taper_r_dist) return c;
	c = mix(mix(c, R, taper_r_strength), c, offset / taper_r_dist);
	return c;
}

void draw_line(vec4 LL, vec4 L, vec4 R, vec4 RR){

	float LAngle, RAngle;


	float OffsetL = gOffset[1].x;
	float OffsetR = gOffset[2].x;
	float OffsetL2 = gOffset[1].y;
	float OffsetR2 = gOffset[2].y;



	if (L == R) return;

	vec4 a;
	vec4 b;
	vec4 c;
	vec4 d;
	vec4 Line = R - L;
	vec4 Normal = normalize(vec4(-Line.y, Line.x, 0, 0));

	a = L - use_thickness * Normal * 0.001;
	b = L + use_thickness * Normal * 0.001;
	c = R - use_thickness * Normal * 0.001;
	d = R + use_thickness * Normal * 0.001;

	float lim = use_thickness * 0.002;

	float x_scale = preview_viewport.w / preview_viewport.z;

	if (LL.x < 3e20) {
		vec4 avg = normalize(L - LL) + normalize(R - L);
		if (length(avg) > 0.001) {
			vec4 Tangent = normalize(avg);
			vec4 Minter = normalize(vec4(-Tangent.y, Tangent.x, 0, 0));
			float length = use_thickness / (dot(Minter, Normal)) * 0.001;
			if (length < 4 * lim) {
				Minter.x *= x_scale;
				a = L - length * Minter;
				b = L + length * Minter;
			}
		}
	}

	if (RR.x < 3e20) {
		vec4 avg = normalize(RR - R) + normalize(R - L);
		if (length(avg) > 0.001) {
			vec4 Tangent = normalize(avg);
			vec4 Minter = normalize(vec4(-Tangent.y, Tangent.x, 0, 0));
			float length = use_thickness / (dot(Minter, Normal)) * 0.001;
			if (length < 4 * lim) {
				Minter.x *= x_scale;
				c = R - length * Minter;
				d = R + length * Minter;
			}
		}
	}

	a = MakeLeftTaperLinear(L, a, OffsetL);
	b = MakeLeftTaperLinear(L, b, OffsetL);
	c = MakeLeftTaperLinear(R, c, OffsetR);
	d = MakeLeftTaperLinear(R, d, OffsetR);

	a = MakeRightTaperLinear(L, a, OffsetL2);
	b = MakeRightTaperLinear(L, b, OffsetL2);
	c = MakeRightTaperLinear(R, c, OffsetR2);
	d = MakeRightTaperLinear(R, d, OffsetR2);

	a.w = 1;
	b.w = 1;
	c.w = 1;
	d.w = 1;	

	gl_Position = a;
	EmitVertex();
	gl_Position = b;
	EmitVertex();
	gl_Position = c;
	EmitVertex();
	EndPrimitive();

	gl_Position = c;
	EmitVertex();
	gl_Position = d;
	EmitVertex();
	gl_Position = b;
	EmitVertex();
	EndPrimitive();
}

float factor_to_thickness(float factor){
	float r = (factor - normal_ramp_begin)/(normal_ramp_end - normal_ramp_begin);
	if(r>1) r=1;
	if(r<0) r=0;
	float thickness = normal_effect_inverse==1 ?
					  mix(normal_thickness_begin,normal_thickness_end,r) :
					  mix(normal_thickness_end,normal_thickness_begin,r);
	return thickness;
}

void decide_line_style(int component_id){
	float th=thickness;
	if(normal_mode == 0){
		th=thickness;
	}else if(normal_mode == 1){
		float factor = dot(gNormal[0],normal_direction);
		th = factor_to_thickness(factor);
	}else if(normal_mode == 2){
		float factor = dot(gNormal[0],normal_direction);
		th = factor_to_thickness(factor);
	}

	if (component_id == 0) { out_color = color;              use_thickness = th * enable_contour;                               return; }
	if (component_id == 1) { out_color = crease_color;       use_thickness = th * thickness_crease * enable_crease;             return; }
	if (component_id == 2) { out_color = material_color;     use_thickness = th * thickness_material * enable_material;         return; }
	if (component_id == 3) { out_color = edge_mark_color;    use_thickness = th * thickness_edge_mark * enable_edge_mark;       return; }
	if (component_id == 4) { out_color = intersection_color; use_thickness = th * thickness_intersection * enable_intersection; return; }
}

void main() {
	int level = gLevel[1];

	if (occlusion_level_begin > level || occlusion_level_end < level) return;

	float asp1 = output_viewport.z / output_viewport.w;
	float asp2 = preview_viewport.z / preview_viewport.w;
	float x_scale = asp1 / asp2;
	
	vec4 LL = vec4(gl_in[0].gl_Position.xy, 0, 1),
	     L  = vec4(gl_in[1].gl_Position.xy, 0, 1),
	     R  = vec4(gl_in[2].gl_Position.xy, 0, 1),
	     RR = vec4(gl_in[3].gl_Position.xy, 0, 1);
	
	LL.x *= x_scale;
	L.x  *= x_scale;
	R.x  *= x_scale;
	RR.x *= x_scale;

	int type = gType[1];

	decide_line_style(type);

	draw_line(LL, L, R, RR);
}

