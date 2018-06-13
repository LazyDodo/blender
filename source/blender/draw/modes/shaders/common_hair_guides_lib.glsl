/* NOTE: Keep this code in sync with the C version in BKE_hair! */

#ifdef HAIR_SHADER_FIBERS

#define M_PI 3.1415926535897932384626433832795

mat4 translate(vec3 co)
{
	return mat4(1.0, 0.0, 0.0, 0.0,
	            0.0, 1.0, 0.0, 0.0,
	            0.0, 0.0, 1.0, 0.0,
	            co.x, co.y, co.z, 1.0);
}

mat4 rotateX(float angle)
{
	float ca = cos(angle);
	float sa = sin(angle);
	return mat4(1.0, 0.0, 0.0, 0.0,
	            0.0, ca,   sa, 0.0,
	            0.0, -sa,  ca, 0.0,
	            0.0, 0.0, 0.0, 1.0);
}

mat4 rotateY(float angle)
{
	float ca = cos(angle);
	float sa = sin(angle);
	return mat4(ca,  0.0,  sa, 0.0,
	            0.0, 1.0, 0.0, 0.0,
	            -sa, 0.0,  ca, 0.0,
	            0.0, 0.0, 0.0, 1.0);
}

mat4 rotateZ(float angle)
{
	float ca = cos(angle);
	float sa = sin(angle);
	return mat4(ca,  sa,  0.0, 0.0,
	            -sa, ca,  0.0, 0.0,
	            0.0, 0.0, 1.0, 0.0,
	            0.0, 0.0, 0.0, 1.0);
}

/* Hair Displacement */

/* Note: The deformer functions below calculate a new location vector
 * as well as a new direction (aka "normal"), using the partial derivatives of the transformation.
 * 
 * Each transformation function can depend on the location L as well as the curve parameter t:
 *
 *         Lnew = f(L, t)
 *  => dLnew/dt = del f/del L * dL/dt + del f/del t
 *
 * The first term is the Jacobian of the function f, dL/dt is the original direction vector.
 * Some more information can be found here:
 * https://developer.nvidia.com/gpugems/GPUGems/gpugems_ch42.html
 */

struct DeformParams
{
	/* Length where strand reaches final thickness */
	float taper_length;
	/* Relative strand thickness at the tip.
	 * (0.0, 1.0]
	 * 0.0 : Strand clumps into a single line
	 * 1.0 : Strand does not clump at all
	 * (> 1.0 is possible but not recommended)
	 */
	float taper_thickness;

	/* Radius of the curls.
	 * >= 0.0
	 */
	float curl_radius;
	/* Steepness of curls
	 * < 0.0 : Clockwise curls
	 * > 0.0 : Anti-clockwise curls
	 */
	float curl_angle;
};

void deform_taper(DeformParams params, float t, out float taper, out float dtaper)
{
	/* Uses the right half of the smoothstep function */
	float x = (t + params.taper_length) / params.taper_length;
	float dx = 1.0 / params.taper_length;
	if (x > 2.0)
	{
		x = 2.0;
		dx = 0.0;
	}
	taper = 0.5 * x * x * (3 - x) - 1.0;
	dtaper = 1.5 * x * (2.0 - x) * dx;
}

/* Hairs tend to stick together and run in parallel.
 * The effect increases with distance from the root,
 * as the stresses pulling fibers apart decrease.
 */
void deform_clump(DeformParams params, float t, mat4 target_matrix,
                  inout vec3 co, inout vec3 tang)
{
	float taper, dtaper;
	deform_taper(params, t, taper, dtaper);
	float factor = (1.0 - params.taper_thickness) * taper;
	float dfactor = (1.0 - params.taper_thickness) * dtaper;
	
	vec3 target_co = target_matrix[3].xyz;
	vec3 target_tang = target_matrix[0].xyz;
	vec3 nco = co + (target_co - co) * factor;
	vec3 ntang = normalize(tang + (target_tang - tang) * factor + (target_co - co) * dfactor);

	co = nco;
	tang = ntang;
}

/* Hairs often don't have a circular cross section, but are somewhat flattened.
 * This creates the local bending which results in the typical curly hair geometry.
 */
void deform_curl(DeformParams params, float t,
                 inout mat4 target_matrix)
{
	float pitch = 2.0*M_PI * params.curl_radius * tan(params.curl_angle);
	float angle = t / (params.curl_radius * tan(params.curl_angle));
	mat4 local_mat = rotateX(angle) * translate(vec3(0.0, params.curl_radius, 0.0)) * rotateY(params.curl_angle);
	target_matrix = target_matrix * local_mat;
}

void deform_fiber(DeformParams params,
                  float t, mat4 target_matrix,
                  inout vec3 loc, inout vec3 tang)
{
	//deform_curl(params, t, target_matrix);
	deform_clump(params, t, target_matrix, loc, tang);
}

/*===================================*/
/* Hair Interpolation */

uniform sampler2D fiber_data;

uniform int fiber_start;
uniform int strand_map_start;
uniform int strand_vertex_start;

#define INDEX_INVALID -1

vec2 read_texdata(int offset)
{
	ivec2 offset2 = ivec2(offset % HAIR_SHADER_TEX_WIDTH, offset / HAIR_SHADER_TEX_WIDTH);
	return texelFetch(fiber_data, offset2, 0).rg;
}

mat4 mat4_from_vectors(vec3 nor, vec3 tang, vec3 co)
{
	tang = normalize(tang);
	vec3 xnor = normalize(cross(nor, tang));
	return mat4(vec4(tang, 0.0), vec4(xnor, 0.0), vec4(cross(tang, xnor), 0.0), vec4(co, 1.0));
}

void get_strand_data(int index, out int start, out int count, out DeformParams deform_params)
{
	int offset = strand_map_start + index * 2;
	vec2 a = read_texdata(offset);
	vec2 b = read_texdata(offset + 1);

	start = floatBitsToInt(a.r);
	count = floatBitsToInt(a.g);

	deform_params.taper_length = b.r;
	deform_params.taper_thickness = b.g;
	deform_params.curl_radius = 0.1;
	deform_params.curl_angle = 0.2;
}

void get_strand_vertex(int index, out vec3 co, out vec3 nor, out vec3 tang, out float len)
{
	int offset = strand_vertex_start + index * 5;
	vec2 a = read_texdata(offset);
	vec2 b = read_texdata(offset + 1);
	vec2 c = read_texdata(offset + 2);
	vec2 d = read_texdata(offset + 3);
	vec2 e = read_texdata(offset + 4);

	co = vec3(a.rg, b.r);
	nor = vec3(b.g, c.rg);
	tang = vec3(d.rg, e.r);
	len = e.g;
}

void get_strand_root(int index, out vec3 co)
{
	int offset = strand_vertex_start + index * 5;
	vec2 a = read_texdata(offset);
	vec2 b = read_texdata(offset + 1);

	co = vec3(a.rg, b.r);
}

void get_fiber_data(int fiber_index, out ivec4 parent_index, out vec4 parent_weight, out vec3 pos)
{
	int offset = fiber_start + fiber_index * 6;
	vec2 a = read_texdata(offset);
	vec2 b = read_texdata(offset + 1);
	vec2 c = read_texdata(offset + 2);
	vec2 d = read_texdata(offset + 3);
	vec2 e = read_texdata(offset + 4);
	vec2 f = read_texdata(offset + 5);

	parent_index = ivec4(floatBitsToInt(a.rg), floatBitsToInt(b.rg));
	parent_weight = vec4(c.rg, d.rg);
	pos = vec3(e.rg, f.r);
}

struct ParentCurveResult
{
	vec3 co;
	vec3 nor;
	vec3 tang;
	float len;

	// Only needed for the primary parent curve
	vec3 rootco;
	DeformParams deform_params;
};

bool interpolate_parent_curve(int index, float curve_param, out ParentCurveResult result)
{
	if (index == INDEX_INVALID)
	{
		return false;
	}

	int start, count;
	get_strand_data(index, start, count, result.deform_params);

	get_strand_root(start, result.rootco);
	
#if 0 // Don't have to worry about out-of-bounds segment here, as long as lerpfac becomes 0.0 when curve_param==1.0
	float maxlen = float(count - 1);
	float arclength = curve_param * maxlen;
	int segment = min(int(arclength), count - 2);
	float lerpfac = arclength - min(floor(arclength), maxlen - 1.0);
#else
	float maxlen = float(count - 1);
	float arclength = curve_param * maxlen;
	int segment = int(arclength);
	float lerpfac = arclength - floor(arclength);
#endif
	
	vec3 co0, co1;
	vec3 nor0, nor1;
	vec3 tang0, tang1;
	float len0, len1;
	get_strand_vertex(start + segment, co0, nor0, tang0, len0);
	get_strand_vertex(start + segment + 1, co1, nor1, tang1, len1);

	result.co = mix(co0, co1, lerpfac) - result.rootco;
	result.nor = mix(nor0, nor1, lerpfac);
	result.tang = mix(tang0, tang1, lerpfac);
	result.len = mix(len0, len1, lerpfac);

	return true;
}

void interpolate_vertex(int fiber_index, float curve_param,
	                    out vec3 co, out vec3 tang, out float len,
	                    out mat4 target_matrix, out DeformParams deform_params)
{
	ivec4 parent_index;
	vec4 parent_weight;
	vec3 rootco;
	get_fiber_data(fiber_index, parent_index, parent_weight, rootco);

	co = vec3(0.0);
	tang = vec3(0.0);
	len = 0.0;
	target_matrix = mat4(1.0);

	ParentCurveResult p1, p2, p3, p4;
	if (interpolate_parent_curve(parent_index.x, curve_param, p1))
	{
		co += parent_weight.x * p1.co;
		tang += parent_weight.x * normalize(p1.tang);
		len += parent_weight.x * p1.len;
	}
	if (interpolate_parent_curve(parent_index.y, curve_param, p2))
	{
		co += parent_weight.y * p2.co;
		tang += parent_weight.y * normalize(p2.tang);
		len += parent_weight.y * p2.len;
	}
	if (interpolate_parent_curve(parent_index.z, curve_param, p3))
	{
		co += parent_weight.z * p3.co;
		tang += parent_weight.z * normalize(p3.tang);
		len += parent_weight.z * p3.len;
	}
	if (interpolate_parent_curve(parent_index.w, curve_param, p4))
	{
		co += parent_weight.w * p4.co;
		tang += parent_weight.w * normalize(p4.tang);
		len += parent_weight.w * p4.len;
	}
	
	co += rootco;
	tang = normalize(tang);
	target_matrix = mat4_from_vectors(p1.nor, p1.tang, p1.co + p1.rootco);
	deform_params = p1.deform_params;
}

void hair_fiber_get_vertex(
        int fiber_index, float curve_param,
        bool is_persp, vec3 camera_pos, vec3 camera_z,
        out vec3 pos, out vec3 tang, out vec3 binor,
        out float time, out float thickness, out float thick_time)
{
	float len;
	mat4 target_matrix;
	DeformParams deform_params;
	interpolate_vertex(fiber_index, curve_param, pos, tang, len, target_matrix, deform_params);

	vec3 camera_vec = (is_persp) ? pos - camera_pos : -camera_z;
	binor = normalize(cross(camera_vec, tang));

	deform_fiber(deform_params, len, target_matrix, pos, tang);

	time = curve_param;
	thickness = hair_shaperadius(hairRadShape, hairRadRoot, hairRadTip, time);

	// TODO use the uniform for hairThicknessRes
	int hairThicknessRes = 2;
	if (hairThicknessRes > 1) {
		thick_time = float(gl_VertexID % hairThicknessRes) / float(hairThicknessRes - 1);
		thick_time = thickness * (thick_time * 2.0 - 1.0);

		pos += binor * thick_time;
	}
}

#endif /*HAIR_SHADER_FIBERS*/
