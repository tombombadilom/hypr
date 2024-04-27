#version 320 es

uniform float iTime;
uniform vec2 iResolution;

#define PI 3.14159265359
#define SUN_COLOR vec3(1., .956, .839)
#define SKY_COLOR vec3(0.3, 0.5, 0.85)
#define HORIZON_COLOR vec3(0.4, 0.65, 1.0)
#define ITERATIONS 10
#define EPSILON 0.0001
#define FREQUENCY 0.8
#define LACUNARITY 2.0

float hash(in vec2 uv) {
    return fract(sin(dot(uv, vec2(12.9898, 78.233))) * 43758.5453123);
}

float noise(in vec2 uv) {
    vec2 i = floor(uv);
    vec2 f = fract(uv);
    float a = hash(i);
    float b = hash(i + vec2(1.0, 0.0));
    float c = hash(i + vec2(0.0, 1.0));
    float d = hash(i + vec2(1.0, 1.0));
    vec2 u = f * f * (3.0 - 2.0 * f);
    return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float fbm(in vec2 uv) {
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = FREQUENCY;
    for (int i = 0; i < ITERATIONS; ++i) {
        value += amplitude * noise(frequency * uv);
        amplitude *= 0.5;
        frequency *= LACUNARITY;
    }
    return value;
}

vec3 getNormal(in vec3 p) {
    vec2 e = vec2(EPSILON, 0.0);
    float centralHeight = fbm(p.xz);
    return normalize(vec3(
        centralHeight - fbm(p.xz - e.xy),
        2.0 * EPSILON,
        centralHeight - fbm(p.xz - e.yx)
    ));
}

float rayMarching(in vec3 ro, in vec3 rd) {
    const float tMax = 20.0;
    float t = 0.1;
    for (int i = 0; i < ITERATIONS; ++i) { // Reduced iterations for performance
        vec3 pos = ro + t * rd;
        float h = pos.y - fbm(pos.xz);
        if (abs(h) < EPSILON || t > tMax) break;
        t += h * 0.4;
    }
    return t;
}

vec3 calculateLighting(in vec3 p, in vec3 normal, in vec3 lightDir, in vec3 viewDir) {
    vec3 diff = max(dot(normal, lightDir), 0.) * SUN_COLOR;
    vec3 refl = reflect(-lightDir, normal);
    float spec = pow(max(dot(refl, viewDir), 0.), 16.);
    return diff + spec * SUN_COLOR;
}

mat3 lookAt(in vec3 origin, in vec3 target, in vec3 up) {
    vec3 w = normalize(target - origin);
    vec3 u = normalize(cross(up, w));
    vec3 v = cross(w, u);
    return mat3(u, v, w);
}

vec3 camerapath(in float t) {
    return vec3(-13.0 + 3.5 * cos(t), 3.3, -1.1 + 2.4 * cos(2.4 * t + 2.0));
}

void main() {
    vec2 uv = (gl_FragCoord.xy - iResolution.xy * 0.5) / iResolution.y;
    const vec3 lightDir = normalize(vec3(-0.8, 0.15, -0.3));

    vec3 camPos = camerapath(iTime);
    vec3 camTarget = vec3(1., 1., 4.);
    vec3 camStep = vec3(lightDir.x, 0., lightDir.z) * iTime;
    camPos += camStep;
    camTarget += camStep;

    mat3 mat = lookAt(camPos, camTarget, vec3(0.0, 1.0, 0.0));
    vec3 rd = normalize(mat * vec3(uv, 1.0));

    float t = rayMarching(camPos, rd);

    vec3 col;
    if (t > 20.0) {
        float sundot = clamp(dot(rd, lightDir), 0.0, 1.0);
        col = mix(SKY_COLOR, SUN_COLOR, pow(sundot, 5.0) * 0.1
                                    + pow(sundot, 64.0) * 0.05
                                    + pow(sundot, 512.0) * 0.025);
    } else {
        vec3 p = camPos + rd * t;
        vec3 normal = getNormal(p);
        col = calculateLighting(p, normal, lightDir, -rd);
    }

    col = pow(col, vec3(2.2)); // sRGB gamma correction
    gl_FragColor = vec4(col, 1.0);
}
