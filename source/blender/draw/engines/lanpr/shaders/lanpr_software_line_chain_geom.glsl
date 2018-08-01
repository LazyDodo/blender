layout(lines) in;
layout(triangle_strip, max_vertices = 6) out;


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

out vec4 out_color;

float use_thickness;

void draw_line(vec4 p1, vec4 p2){

	vec4 Line = p2 - p1;
	vec4 Normal = normalize(vec4(-Line.y, Line.x, 0, 0));

	vec4 a, b, c, d;

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

void decide_color_and_thickness(float component_id){
	if (component_id < 1.5) { out_color = color;              use_thickness = thickness;                          return; }
	if (component_id < 2.5) { out_color = crease_color;       use_thickness = thickness * thickness_crease;       return; }
	if (component_id < 3.5) { out_color = material_color;     use_thickness = thickness * thickness_material;     return; }
	if (component_id < 4.5) { out_color = edge_mark_color;    use_thickness = thickness * thickness_edge_mark;    return; }
	if (component_id < 5.5) { out_color = intersection_color; use_thickness = thickness * thickness_intersection; return; }
}

void main() {
	vec4 p1 = vec4(gl_in[0].gl_Position.xy, 0, 1);
	vec4 p2 = vec4(gl_in[1].gl_Position.xy, 0, 1);

	decide_color_and_thickness(gl_in[0].gl_Position.z);

	draw_line(p1, p2);
}