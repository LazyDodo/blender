
out vec4 FragColor;

uniform sampler2D strokeColor;
uniform sampler2D strokeDepth;

uniform int blur[2];

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	/* apply blurring, using a 9-tap filter with predefined gaussian weights */
	/* depth */
	float outdepth = 0;
    outdepth += texelFetch(strokeDepth, ivec2(uv.x - 1.0 * blur[0], uv.y + 1.0 * blur[1]), 0).r * 0.0947416;
    outdepth += texelFetch(strokeDepth, ivec2(uv.x - 0.0 * blur[0], uv.y + 1.0 * blur[1]), 0).r * 0.118318;
    outdepth += texelFetch(strokeDepth, ivec2(uv.x + 1.0 * blur[0], uv.y + 1.0 * blur[1]), 0).r * 0.0947416;
    outdepth += texelFetch(strokeDepth, ivec2(uv.x - 1.0 * blur[0], uv.y + 0.0 * blur[1]), 0).r * 0.118318;

    outdepth += texelFetch(strokeDepth, ivec2(uv.x, uv.y), 0).r * 0.147761;

    outdepth += texelFetch(strokeDepth, ivec2(uv.x + 1.0 * blur[0], uv.y + 0.0 * blur[1]), 0).r * 0.118318;
    outdepth += texelFetch(strokeDepth, ivec2(uv.x - 1.0 * blur[0], uv.y - 1.0 * blur[1]), 0).r * 0.0947416;
    outdepth += texelFetch(strokeDepth, ivec2(uv.x + 0.0 * blur[0], uv.y - 1.0 * blur[1]), 0).r * 0.118318;
    outdepth += texelFetch(strokeDepth, ivec2(uv.x + 1.0 * blur[0], uv.y - 1.0 * blur[1]), 0).r * 0.0947416;

	gl_FragDepth = outdepth;

	/* color */	
	vec4 outcolor = vec4(0.0);
    outcolor += texelFetch(strokeColor, ivec2(uv.x - 1.0 * blur[0], uv.y + 1.0 * blur[1]), 0) * 0.0947416;
    outcolor += texelFetch(strokeColor, ivec2(uv.x - 0.0 * blur[0], uv.y + 1.0 * blur[1]), 0) * 0.118318;
    outcolor += texelFetch(strokeColor, ivec2(uv.x + 1.0 * blur[0], uv.y + 1.0 * blur[1]), 0) * 0.0947416;
    outcolor += texelFetch(strokeColor, ivec2(uv.x - 1.0 * blur[0], uv.y + 0.0 * blur[1]), 0) * 0.118318;

    outcolor += texelFetch(strokeColor, ivec2(uv.x, uv.y), 0) * 0.147761;

    outcolor += texelFetch(strokeColor, ivec2(uv.x + 1.0 * blur[0], uv.y + 0.0 * blur[1]), 0) * 0.118318;
    outcolor += texelFetch(strokeColor, ivec2(uv.x - 1.0 * blur[0], uv.y - 1.0 * blur[1]), 0) * 0.0947416;
    outcolor += texelFetch(strokeColor, ivec2(uv.x + 0.0 * blur[0], uv.y - 1.0 * blur[1]), 0) * 0.118318;
    outcolor += texelFetch(strokeColor, ivec2(uv.x + 1.0 * blur[0], uv.y - 1.0 * blur[1]), 0) * 0.0947416;

	FragColor = clamp(outcolor, 0, 1.0);
}
