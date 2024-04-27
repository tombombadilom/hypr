#version 320 es
precision highp float;

// Converts a vec2 to a vec3 with a z-component of 0.0
vec3 vec_tovec3(vec2 v) {
    return vec3(v.x, v.y, 0.0);
}

// Hash function converting a vec2 to a single float value
float hash1_2(vec2 v) {
    vec3 v3 = vec3(v, 0.0);
    v3 = fract(v3 * 0.1031);
    v3 += dot(v3, v3 + 33.33);
    return fract((v3.x + v3.y) * v3.z);
}

// Another hash function that returns a vec2 from input vec2
vec2 hash2_1(vec2 v) {
    vec3 v3 = vec3(v, 0.0);
    v3 *= vec3(0.1031, 0.103, 0.0973);
    v3 += dot(v3, v3 + 33.33);
    return fract((v3.xy + v3.yz) * v3.xz);
}

// Smooth interpolation between four values within unit square
float noise_noisemix2(float a, float b, float c, float d, vec2 f) {
    vec2 u = f * f * (3.0 - 2.0 * f);
    return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

// White noise based on hashing input coordinate
float noise_white_1(vec2 p) {
    return hash1_2(p);
}

// Value noise function
float noise_value_1(vec2 p) {
    // Grid cell coordinates and fractional part
    vec2 i = floor(p);
    vec2 f = fract(p);

    vec2 I = i + 1.0;

    // Hashed noise values at grid cell corners
    float a = hash1_2(i);
    float b = hash1_2(vec2(I.x, i.y));
    float c = hash1_2(vec2(i.x, I.y));
    float d = hash1_2(I);

    // Noise interpolation
    return noise_noisemix2(a, b, c, d, f);
}

// Gradient noise function
float noise_gradient_1(vec2 p) {
    // Grid cell coordinates and fractional part
    vec2 i = floor(p);
    vec2 f = fract(p);

    // Pseudorandom gradients at cell corners
    vec4 rnd = vec4(hash1_2(i), hash1_2(i + vec2(1.0, 0.0)),
                    hash1_2(i + vec2(0.0, 1.0)), hash1_2(i + vec2(1.0, 1.0)));

    // Interpolation weights
    vec2 u = f * f * (3.0 - 2.0 * f);

    // Interpolated noise result
    return mix(mix(rnd.x, rnd.y, u.x), mix(rnd.z, rnd.w, u.x), u.y);
}

// Fractal Brownian Motion function using value noise
float fbm_value_1(vec2 p, int oct) {
    float f = 0.0;
    float amp = 0.5;
    for (int i = 0; i < oct; i++) {
        f += amp * noise_value_1(p);
        amp *= 0.5;
        p *= 2.0;
    }
    return f;
}

// Fractal Brownian Motion function using gradient noise
float fbm_gradient_1(vec2 p, int oct) {
    float f = 0.0;
    float amp = 0.5;
    for (int i = 0; i < oct; i++) {
        f += amp * noise_gradient_1(p);
        amp *= 0.5;
        p *= 2.0;
    }
    return f;
}

// Fractal Brownian Motion using white noise
float fbm_white_1(vec2 p, int oct) {
    float f = 0.0;
    float amp = 0.5;
    for (int i = 0; i < oct; i++) {
        f += amp * noise_white_1(p);
        amp *= 0.5;
        p *= 2.0;
    }
    return f;
}

// Uniform variable for screen resolution
uniform vec2 iResolution;

// Output color variable
out vec4 fragColor;

void main() {
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = gl_FragCoord.xy / iResolution.xy;

    // Compute Fractal Brownian Motion for different noise types
    const int octaves = 6;
    float fbm_value = fbm_value_1(uv * 10.0, octaves);
    float fbm_gradient = fbm_gradient_1(uv * 10.0, octaves);
    float fbm_white = fbm_white_1(uv * 10.0, octaves);

    // Map noise to height and create a corresponding color
    float height = fbm_value * 5.0;
    vec3 color = vec3(fbm_gradient, fbm_gradient, fbm_gradient);
    color = mix(color, vec3(height * 0.01, height * 0.1, height * 0.5), 0.5);

    // Output color to screen
    fragColor = vec4(color, 1.0);
}
