
uniform mat4 ModelViewProjectionMatrix;
uniform mat4 ModelMatrix;
uniform mat4 ModelViewMatrix;

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

#ifdef HAIR_SHADER_FIBERS
	vec3 pos;
	vec3 nor;
	vec2 view_offset;
	hair_fiber_get_vertex(fiber_index, curve_param, ModelViewMatrix, pos, nor, view_offset);
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	gl_Position.xy += view_offset * gl_Position.w;
#else
	float time, thick_time, thickness;
	vec3 pos, tan, binor;
	hair_get_pos_tan_binor_time(
	        (ProjectionMatrix[3][3] == 0.0),
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
