/* NOTE: Keep this code in sync with the C version in BKE_hair! */

#ifdef HAIR_SHADER_FIBERS

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

/* Float with single derivative */
struct DualFloat
{
	float v;
	float dv;
};

/* Vector with single derivative */
struct DualVec3
{
	vec3 v;
	vec3 dv;
};

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

void calc_taper_factor(DeformParams params, float t, out DualFloat taper)
{
	/* Uses the right half of the smoothstep function */
	float x = (t + params.taper_length) / params.taper_length;
	float dx = 1.0 / params.taper_length;
	if (x > 2.0)
	{
		x = 2.0;
		dx = 0.0;
	}

	taper.v = 0.5 * x * x * (3 - x) - 1.0;
	taper.dv = 1.5 * x * (2.0 - x) * dx;
}

/* Hairs tend to stick together and run in parallel.
 * The effect increases with distance from the root,
 * as the stresses pulling fibers apart decrease.
 */
void deform_clump(DualFloat taper, DualVec3 target, float thickness, inout vec3 co, inout vec3 tang)
{
	DualFloat factor;
	factor.v = taper.v * thickness;
	factor.dv = taper.dv * thickness;

	vec3 nco;
	vec3 ntang;
	nco = co + (target.v - co) * factor.v;
	ntang = normalize(tang
	                  + (target.dv - tang) * factor.v
	                  + (target.v - co) * factor.dv);

	co = nco;
	tang = ntang;
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
	result.nor = normalize(mix(nor0, nor1, lerpfac));
	result.tang = normalize(mix(tang0, tang1, lerpfac));
	result.len = mix(len0, len1, lerpfac);

	return true;
}

void interpolate_vertex(int fiber_index, float curve_param,
	                    out vec3 co, out vec3 tang)
{
	ivec4 parent_index;
	vec4 parent_weight;
	vec3 rootco;
	get_fiber_data(fiber_index, parent_index, parent_weight, rootco);

	co = vec3(0.0);
	tang = vec3(0.0);

	ParentCurveResult p1, p2, p3, p4;
	if (interpolate_parent_curve(parent_index.x, curve_param, p1))
	{
		vec3 defco = p1.co + rootco;
		vec3 deftang = p1.tang;

		DualFloat taper;
		calc_taper_factor(p1.deform_params, p1.len, taper);
		if (taper.v > 0.0)
		{
			DualVec3 target;
			target.v = p1.co + p1.rootco;
			target.dv = p1.tang;

			deform_clump(taper, target, p1.deform_params.taper_thickness, defco, deftang);

			/* Modulate weights by taper factor,
			 * so influence of the parent increases with taper
			 */
			parent_weight.x   = parent_weight.x   * (1.0 - taper.v) + taper.v;
			parent_weight.yzw = parent_weight.yzw * (1.0 - taper.v);
		}

		co += parent_weight.x * (defco - rootco);
		tang += parent_weight.x * normalize(deftang);
	}
	if (interpolate_parent_curve(parent_index.y, curve_param, p2))
	{
		co += parent_weight.y * p2.co;
		tang += parent_weight.y * p2.tang;
	}
	if (interpolate_parent_curve(parent_index.z, curve_param, p3))
	{
		co += parent_weight.z * p3.co;
		tang += parent_weight.z * p3.tang;
	}
	if (interpolate_parent_curve(parent_index.w, curve_param, p4))
	{
		co += parent_weight.w * p4.co;
		tang += parent_weight.w * p4.tang;
	}
	
	co += rootco;
	tang = normalize(tang);
}

void hair_fiber_get_vertex(
        int fiber_index, float curve_param,
        bool is_persp, vec3 camera_pos, vec3 camera_z,
        out vec3 pos, out vec3 tang, out vec3 binor,
        out float time, out float thickness, out float thick_time)
{
	interpolate_vertex(fiber_index, curve_param, pos, tang);

	vec3 camera_vec = (is_persp) ? pos - camera_pos : -camera_z;
	binor = normalize(cross(camera_vec, tang));

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
