
uniform mat4 ModelViewProjectionMatrix;
uniform mat4 ModelMatrix;
uniform mat4 ModelViewMatrix;
uniform mat4 ModelViewMatrixInverse;
uniform mat4 ProjectionMatrix;

/* keep in sync with DRWManager.view_data */
layout(std140) uniform clip_block {
	vec4 ClipPlanes[1];
};

#ifdef HAIR_SHADER

#ifdef HAIR_SHADER_FIBERS
in int fiber_index;
in float curve_param;
#endif

#else
in vec3 pos;
#endif

void main()
{
#ifdef HAIR_SHADER
	bool is_persp = (ProjectionMatrix[3][3] == 0.0);

#ifdef HAIR_SHADER_FIBERS
	float time, thick_time, thickness;
	vec3 pos, tang, binor;
	hair_fiber_get_vertex(
	        fiber_index, curve_param,
	        is_persp, ModelViewMatrixInverse[3].xyz, ModelViewMatrixInverse[2].xyz,
	        pos, tang, binor,
	        time, thickness, thick_time);

	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	vec4 worldPosition = ModelMatrix * vec4(pos, 1.0);
#else
	float time, thick_time, thickness;
	vec3 pos, tan, binor;
	hair_get_pos_tan_binor_time(
	        is_persp,
	        ViewMatrixInverse[3].xyz, ViewMatrixInverse[2].xyz,
	        pos, tan, binor, time, thickness, thick_time);

	gl_Position = ViewProjectionMatrix * vec4(pos, 1.0);
	vec4 worldPosition = vec4(pos, 1.0);
#endif

#else
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	vec4 worldPosition = (ModelMatrix * vec4(pos, 1.0));
#endif

#ifdef CLIP_PLANES
	gl_ClipDistance[0] = dot(vec4(worldPosition.xyz, 1.0), ClipPlanes[0]);
#endif
	/* TODO motion vectors */
}
