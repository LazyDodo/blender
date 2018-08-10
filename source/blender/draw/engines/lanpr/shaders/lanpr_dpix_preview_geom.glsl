layout(points) in;
layout(triangle_strip, max_vertices = 6) out;

uniform sampler2D vert0_tex;//L
uniform sampler2D vert1_tex;//R
uniform sampler2D edge_mask_tex;

//uniform float uValue0; // buffer_w
uniform vec4 viewport; // viewport
uniform float depth_offset;

// these are for depth related thickness control;
uniform float line_thickness;
uniform float depth_width_influence;
uniform float depth_width_curve;
uniform float depth_alpha_influence;
uniform float depth_alpha_curve;
uniform float z_near;
uniform float z_far;

uniform vec4 color;
uniform vec4 background_color;

uniform vec4 crease_color;
uniform vec4 material_color;
uniform vec4 edge_mark_color;
uniform vec4 intersection_color;

uniform float line_thickness_crease;
uniform float line_thickness_material;
uniform float line_thickness_edge_mark;
uniform float line_thickness_intersection;

float use_thickness = line_thickness;

out vec4 out_color;

vec4 use_color = color;

float get_linear_depth(float z){
	float ze = 2.0 * z_near * z_far / (z_far + z_near - z * (z_far - z_near));
	return (ze - z_near) / (z_far - z_near);
}

float curve_01(float z, float factor){
	return pow(z, 1 - factor); // factor is -inf~1-eps
}

vec4 apply_scale(vec4 center, vec4 a){
	float lz = get_linear_depth(center.z);
	float depth_factor = mix(0, curve_01(lz, depth_width_curve), depth_width_influence);

	return mix(a, center, depth_factor);
}

void emit_alpha_pre_mul(vec4 a, int is_crease, float crease_fading){
	float lz = get_linear_depth(a.z);
	float alpha_factor = mix(0, curve_01(lz, depth_alpha_curve), depth_alpha_influence);
	float alpha_crease_fading = alpha_factor;
	if (is_crease > 0) alpha_crease_fading = mix(alpha_factor, 1, crease_fading * 2); // fading=0.5 -> fade all
	out_color = mix(use_color, background_color, alpha_crease_fading);
}

void draw_line(vec4 p1, vec4 p2, int is_crease){

	vec4 Line = p2 - p1;
	vec4 Normal = normalize(vec4(-Line.y, Line.x, 0, 0));

	vec4 a, b, c, d;

	vec4 offset = Normal * use_thickness * 0.001;
	a = p1 + offset;
	b = p1 - offset;
	c = p2 + offset;
	d = p2 - offset;

	a = apply_scale(p1, a);
	b = apply_scale(p1, b);
	c = apply_scale(p2, c);
	d = apply_scale(p2, d);

	gl_Position = vec4(a.xy, a.z - depth_offset, 1); emit_alpha_pre_mul(a, is_crease, p2.w); EmitVertex();
	gl_Position = vec4(b.xy, b.z - depth_offset, 1); emit_alpha_pre_mul(b, is_crease, p2.w); EmitVertex();
	gl_Position = vec4(c.xy, c.z - depth_offset, 1); emit_alpha_pre_mul(c, is_crease, p2.w); EmitVertex();

	gl_Position = vec4(b.xy, b.z - depth_offset, 1); emit_alpha_pre_mul(b, is_crease, p2.w); EmitVertex();
	gl_Position = vec4(c.xy, c.z - depth_offset, 1); emit_alpha_pre_mul(c, is_crease, p2.w); EmitVertex();
	gl_Position = vec4(d.xy, d.z - depth_offset, 1); emit_alpha_pre_mul(d, is_crease, p2.w); EmitVertex();

	EndPrimitive();
}

void main() {
	vec4 p1 = texelFetch(vert0_tex, ivec2(gl_in[0].gl_Position.xy), 0);
	vec4 p2 = texelFetch(vert1_tex, ivec2(gl_in[0].gl_Position.xy), 0);

	if (p1.w == 0 && p2.w == 0) return;

	vec4 edge_mask = texelFetch(edge_mask_tex, ivec2(gl_in[0].gl_Position.xy), 0);

	int is_crease = 0;

	if (edge_mask.g > 0)      { use_color = edge_mark_color;     use_thickness = line_thickness * line_thickness_edge_mark; }
	else if (edge_mask.r > 0) { use_color = material_color;      use_thickness = line_thickness * line_thickness_material;  }
	else if (edge_mask.b > 0) { use_color = intersection_color;  use_thickness = line_thickness * line_thickness_intersection;  }
	else if (p1.w != p2.w)    { use_color = crease_color;        use_thickness = line_thickness * line_thickness_crease; is_crease = 1; }

	draw_line(p1, p2, is_crease);
}