layout(lines) in;
layout(triangle_strip, max_vertices = 6) out;

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

// for line width correction
uniform vec4 output_viewport;
uniform vec4 preview_viewport;

out vec4 out_color;

float use_thickness;

void draw_line(vec4 p1, vec4 p2){

	vec4 Line = p2 - p1;
	vec4 Normal = normalize(vec4(-Line.y, Line.x, 0, 0));

	vec4 a, b, c, d;

	float x_scale = preview_viewport.w / preview_viewport.z;
	Normal.x *= x_scale;
	vec4 offset = Normal * use_thickness * 0.001;

	a = p1 + offset;
	b = p1 - offset;
	c = p2 + offset;
	d = p2 - offset;

	gl_Position = vec4(a.xy, 0, 1); EmitVertex();
	gl_Position = vec4(b.xy, 0, 1); EmitVertex();
	gl_Position = vec4(c.xy, 0, 1); EmitVertex();

	gl_Position = vec4(b.xy, 0, 1); EmitVertex();
	gl_Position = vec4(c.xy, 0, 1); EmitVertex();
	gl_Position = vec4(d.xy, 0, 1); EmitVertex();

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

void decide_color_and_thickness(float component_id){
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

	if (component_id < 1.5) { out_color = color;              use_thickness = th;                          return; }
	if (component_id < 2.5) { out_color = crease_color;       use_thickness = th * thickness_crease;       return; }
	if (component_id < 3.5) { out_color = material_color;     use_thickness = th * thickness_material;     return; }
	if (component_id < 4.5) { out_color = edge_mark_color;    use_thickness = th * thickness_edge_mark;    return; }
	if (component_id < 5.5) { out_color = intersection_color; use_thickness = th * thickness_intersection; return; }
}

void main() {

	float asp1 = output_viewport.z / output_viewport.w;
	float asp2 = preview_viewport.z / preview_viewport.w;
	float x_scale = asp1 / asp2;

	vec4 p1 = vec4(gl_in[0].gl_Position.xy, 0, 1);
	vec4 p2 = vec4(gl_in[1].gl_Position.xy, 0, 1);

	p1.x *= x_scale;
	p2.x *= x_scale;

	decide_color_and_thickness(gl_in[0].gl_Position.z);

	draw_line(p1, p2);
}