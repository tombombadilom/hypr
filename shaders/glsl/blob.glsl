/*
 * Original shader from: https://www.shadertoy.com/view/tslyWf
 */

#ifdef GL_ES
precision mediump float;
#endif

// glslsandbox uniforms
uniform float time;
uniform vec2 resolution;

// shadertoy emulation
#define iTime time
#define iResolution resolution

// --------[ Original ShaderToy begins here ]---------- //
struct ray
{
	vec3 o;
	vec3 d;
	float l;
};
	
vec4 geo(vec3 p)
{
	float sphere = length(p - vec3(sin(iTime * 1.)*0.2,0.,1.5)) - 0.2 + sin(iTime * 3. + p.y * 9.) * 0.0999; //---'Sphere blob'
	float plane = ( p.y + 0.35);
	if(sphere < plane) return (vec4(0.2,0.6,1.0, sphere));
	else return (vec4(0.2,0.4,0.5, plane));
}

vec4 march(ray r)
{
	vec3 col = vec3(1.0);
	for(int i = 0; i < 16; i++)
	{
		vec3 p = r.o + r.d * r.l;
		vec4 g = geo(p);
		r.l += g.w;
		col = g.rgb;
		if(r.l > 8.)
			break;
	}
	return vec4(col, r.l);
}

vec3 normal(vec3 p)
{
	vec2 of = vec2(0.001,0.0);
	vec4 copy = geo(p);
	return normalize(copy.w - vec3(geo(p - of.xyy).w, geo(p - of.yxy).w, geo(p - of.yyx).w));
}

float shadows(vec3 ro, vec3 rd, float b)
{
	float t = 0.01;
	float res = 1.0;
	for(int i = 0; i < 32; i++)
	{
		float g = geo(ro+rd*t).w;
		res = min(res, b * g / t); //---Traditional soft shadows
		t += g;
		if(g < 0.0005|| g > 3.0)
			break;
	}
	return clamp(res,0.1,1.0);
}

float lighting(vec3 p, vec3 lp)
{
	vec3 lPos = normalize(lp - p);
	vec3 n = normal(p);
	
	float light = clamp(dot(n, lPos), 0.6, 1.); //---Lighter ambient light
	return light * shadows(p, normalize(lp), 1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.xy - 0.5;
	uv.x *= iResolution.x/iResolution.y;

    ray r;
	r.o = vec3(0.);
	r.d = vec3(uv, 1.0);
	r.l = 0.;
	vec4 mm = march(r);

	vec3 m = r.o + r.d * mm.w;
	vec3 lPos = vec3(sin(iTime * 1.0)*1.5,1.5,cos(iTime * 1.0)*1.5); //---Light position
	vec3 col = mm.rgb * lighting(m, lPos);
	col *= exp(-0.1 * mm.w * mm.w* mm.w* mm.w)*  3.; //---Fog

    fragColor = vec4(col,1.0);
}
// --------[ Original ShaderToy ends here ]---------- //

void main(void)
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}