#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 resolution;
uniform sampler2D backbuffer;

#define PI 3.14159

void main(){
	vec2 p = (gl_FragCoord.xy - 0.5 * resolution) / min(resolution.x, resolution.y);
	vec2 t = vec2(gl_FragCoord.xy / resolution);
	
	vec3 c = vec3(0);
	
	for(int i = 0; i < 30; i++) {
		float t = 2.0 * PI * float(i) / 30.0 * time * 0.5;
		float x = cos(3.0*t);
		float y = sin(4.0*t);
		vec2 o = 0.40 * vec2(x, y);
		float r = fract(x);
		float g = 1.0 - r;
		c += 0.01 / (length(p-o)) * vec3(r, g, 0.9);
	}
	
	gl_FragColor = vec4(c, 1);
}