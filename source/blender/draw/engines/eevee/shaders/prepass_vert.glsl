
uniform mat4 ModelViewProjectionMatrix;
uniform mat4 ModelMatrix;
uniform mat4 ModelViewMatrix;

/* keep in sync with DRWManager.view_data */
layout(std140) uniform clip_block {
	vec4 ClipPlanes[1];
};

#ifndef HAIR_SHADER_FIBERS
in vec3 pos;
#else
in int fiber_index;
in float curve_param;
#endif

void main()
{
#ifndef HAIR_SHADER_FIBERS
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
#else
	vec3 pos;
	vec3 nor;
	vec2 view_offset;
	hair_fiber_get_vertex(fiber_index, curve_param, ModelViewMatrix, pos, nor, view_offset);
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	gl_Position.xy += view_offset * gl_Position.w;
#endif

#ifdef CLIP_PLANES
	vec4 worldPosition = (ModelMatrix * vec4(pos, 1.0));
	gl_ClipDistance[0] = dot(vec4(worldPosition.xyz, 1.0), ClipPlanes[0]);
#endif
	/* TODO motion vectors */
}
