uniform mat4 ModelViewProjectionMatrix;
uniform mat4 ProjectionMatrix;
uniform mat4 ViewProjectionMatrix;
uniform mat4 ViewMatrixInverse;
uniform mat4 ModelViewMatrixInverse;
uniform mat3 NormalMatrix;

#ifndef HAIR_SHADER
in vec3 pos;
in vec3 nor;
in vec2 uv;
#else /* HAIR_SHADER */
#  ifdef OB_TEXTURE
uniform samplerBuffer u; /* active texture layer */
#  endif
flat out float hair_rand;
#  ifdef HAIR_SHADER_FIBERS
in int fiber_index;
in float curve_param;
#  endif
#endif /* HAIR_SHADER */

#ifdef NORMAL_VIEWPORT_PASS_ENABLED
out vec3 normal_viewport;
#endif

#ifdef OB_TEXTURE
out vec2 uv_interp;
#endif

/* From http://libnoise.sourceforge.net/noisegen/index.html */
float integer_noise(int n)
{
	n = (n >> 13) ^ n;
	int nn = (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff;
	return (float(nn) / 1073741824.0);
}

void main()
{
#ifdef HAIR_SHADER
	bool is_persp = (ProjectionMatrix[3][3] == 0.0);

#  ifdef HAIR_SHADER_FIBERS
	vec2 uv = vec2(0.0); /* TODO */
	float time, thick_time, thickness;
	vec3 pos, tang, binor;
	hair_fiber_get_vertex(
	        fiber_index, curve_param,
	        is_persp, ModelViewMatrixInverse[3].xyz, ModelViewMatrixInverse[2].xyz,
	        pos, tang, binor,
	        time, thickness, thick_time);
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	hair_rand = integer_noise(fiber_index);
#  else
#    ifdef OB_TEXTURE
	vec2 uv = hair_get_customdata_vec2(u);
#    endif
	float time, thick_time, thickness;
	vec3 pos, tang, binor;
	hair_get_pos_tan_binor_time(
	        is_persp, ViewMatrixInverse[3].xyz, ViewMatrixInverse[2].xyz,
	        pos, tang, binor, time, thickness, thick_time);
	gl_Position = ViewProjectionMatrix * vec4(pos, 1.0);
	hair_rand = integer_noise(hair_get_strand_id());
#  endif
	/* To "simulate" anisotropic shading, randomize hair normal per strand. */
	tang = normalize(tang);
	vec3 nor = normalize(cross(binor, tang));
	nor = normalize(mix(nor, -tang, hair_rand * 0.10));
	float cos_theta = (hair_rand*2.0 - 1.0) * 0.20;
	float sin_theta = sqrt(max(0.0, 1.0f - cos_theta*cos_theta));
	nor = nor * sin_theta + binor * cos_theta;
#else /* HAIR_SHADER */
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
#endif /* HAIR_SHADER */

#ifdef OB_TEXTURE
	uv_interp = uv;
#endif

#ifdef NORMAL_VIEWPORT_PASS_ENABLED
	normal_viewport = NormalMatrix * nor;
#  ifndef HAIR_SHADER
	normal_viewport = normalize(normal_viewport);
#  endif
#endif
}
