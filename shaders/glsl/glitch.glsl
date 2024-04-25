#ifdef GL_ES
precision lowp float;
#endif

uniform float time;
uniform vec2 resolution;
uniform vec2 mouse;

uniform sampler2D backbuffer;

//vec2 u_resolution = resolution;
//float time = time;
//vec2 mouse = u_mouse;//vec2(0.6,0.45);//u_mouse/u_resolution;
vec3 position = vec3(7.,-10.,10.);//vec3(-7.,-8.,2);

#define PI_C 3.14159265358979323846

const int MAX_ITERS = 300;
const int MAX_DIST = 500;
const int MAX_REFLECTIONS = 5;
const float ZOOM = 0.60;
const float EPSI = 0.005;
const float SUN_BRIGHTNESS = 10.;
const float SUN_SIZE = 0.05;
const float RAYLIGH_SCATTER_COEF = 0.0;
const bool ENABLE_EXACT_CLOSEST_APPROACH_SEARCH = false;
const bool ENABLE_AMBIENT_OCCLUSION = false;
const bool ENABLE_DAMN_PYRAMID = true;
const bool ENABLE_DUST = true;
float VCJDECLIPSESIZE = 2.0001;
vec3 LIGHT_DIR = normalize(vec3(0.,5.09,-0.01)); // sin(time*0.25)+(1.-.34873)
vec3 AMBIENT = vec3(0.05,0.08,0.1);//vec3(0.392,0.632,1.000);

float random(vec2 p){
	return(fract(sin(fract(p.x*18.879465+p.y*19.05454))*154.));
}

float noise(vec2 uv, vec2 seed){
	uv += seed + time/200000.;
	vec2 id = floor(uv*10.);
	vec2 lc = smoothstep(0.,1.,fract(uv*10.));

	float a = random(id);
	float b = random(id + vec2(1.,0.));
	float c = random(id + vec2(0.,1.));
	float d = random(id + vec2(1.,1.));

	float ud = mix(a,b,lc.x);
	float lr = mix(c,d,lc.x);
	float fin = mix(ud,lr,lc.y);
	return fin;
}

float octaves(vec2 uv, vec2 seed, int octaves){
	float amp = 0.5;
	float f = 0.;
	for(int i = 1; i < 15;i++)
	{
		f+=noise(uv, seed)*amp;
		uv *= 2.;
		amp *= 0.5;
		if(octaves == i)
			return f;
	}
	return f;
}

float sdRoundBox( vec3 p, vec3 b, float r ) {
  vec3 d = abs(p) - b;
  return length(max(d,0.0)) - r + min(max(d.x,max(d.y,d.z)),0.0);
}

float sdBox(vec3 p, vec3 b) {
  vec3 d = abs(p) - b;
  return length(max(d,0.0)) + min(max(d.x,max(d.y,d.z)),0.0);}

float sdPlane(vec3 p, vec4 n) {
	return dot(p, n.xyz) + n.w;
}

float sdSphere(vec3 p, float r) {
	return length(p)-r;
}

float sdCylinderY(vec3 p, vec3 c) {
  return length(p.xz-c.xy)-c.z;
}

float sdCylinderZ(vec3 p, vec3 c) {
  return length(p.xy-c.xy)-c.z;
}

float sdCylinderX(vec3 p, vec3 c) {
  return length(p.zy-c.xy)-c.z;
}

float smoothmin(float a, float b, float k) {
  float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
  return mix(b, a, h) - k * h * (1.0 - h);
}

float smoothmax(float a, float b, float k){
  return abs(a-b)+smoothmin(a,b,k);
}

float inf_norma(vec3 a){
	return max(a.x,max(a.y,a.z));
}

struct map_types {
	int collision_map;
	int reflection_map;
	int r_map;
	int g_map;
	int b_map;
};
const map_types maps = map_types(1,2,3,4,5);

float trueAngle(float x, float y){
	if(x>0.)
		return atan(y/x);
	if(x<0.)
		return PI_C + atan(y/x);
	if(y>0.)
		return PI_C/2.;
	return 3.*PI_C/2.;
}

vec3 rotateZ(vec3 v, float angle){
	float r = length(v.xz);
	float ang = trueAngle(v.x,v.z);
	ang -= angle;
	return vec3(r*cos(ang), v.y, r*sin(ang));
}

vec3 rotateY(vec3 v, float angle){
	float r = length(v.yz);
	float ang = trueAngle(v.y,v.z);
	ang -= angle;
	return vec3(v.x, r*cos(ang), r*sin(ang));
}

vec3 rotateX(vec3 v, float angle){
	float r = length(v.xy);
	float ang = trueAngle(v.x,v.y);
	ang -= angle;
	return vec3(r*cos(ang), r*sin(ang), v.z);
}

const float wtf = 25.;
const float inverse_wtf = 1./wtf;
const float half_of_wtf= wtf*0.5;

vec3 fractalizer(vec3 p){
	p = fract((p + half_of_wtf)*inverse_wtf)*wtf - half_of_wtf; 
	
	return p;
}

vec2 fractalizer(vec2 p){
	p = fract((p + half_of_wtf)*inverse_wtf)*wtf - half_of_wtf; 
	
	return p;
}

float collision_map(vec3 p)
{
	/*const float time_loop = 100.;
	float fract_time = fract(time / time_loop) * time_loop;
	float inverse_fract_time = time_loop - fract_time;
	float relative_time = fract_time / time_loop;
	
	float combined_z = 0.;
	combined_z -= length(sin(p.xy*0.02 + fract_time*0.25)) * relative_time + (1. - relative_time) * length(sin(p.xy*0.02 + inverse_fract_time*0.25));
	combined_z -= length(0.5*cos(p.yx*0.04 + fract_time*0.25)) * relative_time + (1. - relative_time) * length(0.5*cos(p.yx*0.04 + inverse_fract_time*0.25));
	combined_z -= length(cos(p.z*0.07 + fract_time*0.25)) * relative_time + (1. - relative_time) * length(cos(p.z*0.07 + inverse_fract_time*0.25));
	*/
	float combined_z = length(sin(p.xy*0.02 + time*0.25));
	combined_z -= length(0.5*cos(p.yx*0.04 + time*0.25));
	combined_z -= length(cos(p.z*0.07 + time*0.25));
	
	p.z -= combined_z;
	float plane = sdPlane(p, vec4(0,0,1,0));
	
	return plane;
}

float reflection_map(vec3 p){
	//float plane = sdPlane(p, vec4(0,0,1,0));
	//p.xy=fractalizer(p.xy);
	
	//float pool = sdSphere(p - vec3(half_of_wtf, 0., 0.), half_of_wtf - 5. - 10.*EPSI);
	//pool = min(sdSphere(p - vec3(-half_of_wtf, 0., 0.), half_of_wtf - 5. - 10.*EPSI), pool);
	//if(pool < -2.*EPSI)
	if(true || p.z < -9.45)
		return 0.71;
	return .0;
}

vec3 color_map(vec3 p){
	return vec3(1);
	/*p.xy=fractalizer(p.xy);
	float ceil1 = sdBox(rotateY(p - vec3(half_of_wtf,half_of_wtf,ceiling_height), -0.2*PI_C), vec3(wtf ,half_of_wtf/tan(0.2*PI_C) ,0.75));
	float ceil2 = sdBox(rotateY(p - vec3(half_of_wtf,-half_of_wtf,ceiling_height), 0.2*PI_C), vec3(wtf ,half_of_wtf/tan(0.2*PI_C),0.75));
	
	if(min(ceil1, ceil2) < 2.*EPSI)
		return vec3(1);
	
	if(p.z>1.)
		return vec3(0.5, 0.5, 0.5);
	return vec3(1);*/
}



float general_map(int map_type, vec3 p) {
	if(map_type == maps.collision_map)
		return collision_map(p);
	if(map_type == maps.reflection_map)
		return reflection_map(p);
	if(map_type == maps.r_map)
		return color_map(p).x;
	if(map_type == maps.g_map)
		return color_map(p).x;
	if(map_type == maps.b_map)
		return color_map(p).x;
	
	return -1.;
}

vec3 normal(int map_type, vec3 p) {
	return normalize(vec3(
		general_map(map_type, vec3(p.x + EPSI, p.y, p.z)) - general_map(map_type, vec3(p.x - EPSI, p.y, p.z)),
		general_map(map_type, vec3(p.x, p.y + EPSI, p.z)) - general_map(map_type, vec3(p.x, p.y - EPSI, p.z)),
		general_map(map_type, vec3(p.x, p.y, p.z  + EPSI)) - general_map(map_type, vec3(p.x, p.y, p.z - EPSI))
	));
}

float that_formula(float value){
	const float ssqt = 7.*sqrt(2.);
	return (-ssqt*value + 2.*sqrt(value*value*25. + pow(ssqt+10.,2.)*0.25*(0.25)))/(ssqt + 10.);
}

struct max_approach_point
{
	vec3 position;
	float distance_to_map;
	float distance_from_origin;
	float _debug;
};

max_approach_point rayMarcher(int map_type, 
							  in vec3 start,
							  in vec3 direction,
							  out float dist,
							  out int iters,
							  bool _fill_map_struct){
	max_approach_point _map;
	_map.distance_to_map = 1./0.;
	_map.distance_from_origin = 0.;
	_map._debug = 0.;
	
	dist = 0.;
	vec3 p = start;
	float prev_diff = 0.;
	
	vec3 outway_path = normal(map_type, p);
	float colinearity_coef = dot(outway_path, direction);
	if(colinearity_coef > EPSI)
		p += direction * EPSI/colinearity_coef;
	_map.position = p;
	
	for(int a=0;a<MAX_ITERS;a++){
		float diff = general_map(map_type, p);
		p += direction*diff;
		dist += diff;
		iters = a;
		
		if(prev_diff >= diff && _fill_map_struct)
		{
			if(_map.distance_to_map > diff)
			{
				_map.distance_to_map = diff;
				_map.distance_from_origin = dist;
				_map.position = p;
			}
		}
		
		if(diff<=EPSI && prev_diff >= diff){
			break;
		}
		if(dist>=float(MAX_DIST)){
			break;
		}
		
		prev_diff = diff;
	}
	
	if(_fill_map_struct && ENABLE_EXACT_CLOSEST_APPROACH_SEARCH && _map.distance_from_origin != 0.)
	{
		float dfo = _map.distance_from_origin;
		vec3 _p = _map.position;
		float last_distance = general_map(map_type, _p);
		for(int a = 0; a < 10; a++){
			if(last_distance < EPSI)
				break;
			
			vec3 n = normal(map_type, _p);
			vec3 approx_solid_edge = _p - n * last_distance;
			
			float t = (dot(direction, approx_solid_edge - _p));
			_p += direction * t;
			last_distance = general_map(map_type, _p);
			
			if(abs(t) < EPSI * 0.1)
				break;
			_map._debug = t;
		}
		_map.distance_from_origin = length(_p - start);
		_map.position = _p;
		_map.distance_to_map = last_distance;
	}
	
	return _map;
}

vec4 normalize_color(vec4 color){
	float coef = inf_norma(color.xyz);
	if(coef < 1.) coef = 1.;
	return color / coef;
}

vec3 normalize_color(vec3 color){
	float coef = inf_norma(color);
	if(coef < 1.) coef = 1.;
	return color / coef;
}

vec3 get_color_with_shadows(vec3 p);
vec3 sky_color_in_direction(vec3 p);
vec3 get_diffuse_lighting(vec3 p);
vec4 reflected_raymarching(vec3 start_position, vec3 direction);
	
float skewed_abs_pow(float value, float neg_pow, float pos_pow)
{
	if(value < 0.)
		return pow(-value, neg_pow);
	return pow(value, pos_pow);
}

vec3 smooth_normalize(vec3 p)
{
	float _length = length(p); 
	vec3 normalized = normalize(p);
	float combination_coef = _length / (_length + 1.);
	return (1. - combination_coef) * p + combination_coef * normalized;
}

vec3 sky_color_in_direction(vec3 p){
	p = normalize(p);
	float height = PI_C/2. - acos(p.z/length(p));
	float light_height = PI_C/2. - acos(LIGHT_DIR.z);
	float dawn_color_range = PI_C/10.;
	
	float aligment = dot(p, LIGHT_DIR);
	const vec3 dawn_color = vec3(0.008,0.050,0.290);
	
	vec3 sun_color = vec3(0.7,0.9,1);
	float sun = 0.;
	
	float dawn_color_mul = 0.;
	if(abs(light_height) < dawn_color_range)
		dawn_color_mul = (dawn_color_range - abs(light_height))/dawn_color_range;
	
	sun_color = sun_color + dawn_color_mul*dawn_color;
	sun_color = normalize_color(sun_color);
	
	float aligment_angle = acos(aligment);
	float residual_sun_aligment = (SUN_SIZE - aligment_angle);
	bool is_sun = aligment_angle < SUN_SIZE;
	sun = exp(22.*residual_sun_aligment); 
	
	vec3 dirdiff = LIGHT_DIR - p;
	
	if(VCJDECLIPSESIZE >= 0.)
	{
		float residual_sun_aligment = (SUN_SIZE - aligment_angle) * VCJDECLIPSESIZE;
		
		float angle = trueAngle(dirdiff.x, -dirdiff.z);
		float puff = octaves(vec2(residual_sun_aligment*6., angle*3.), vec2(time, 0.)*0.2, 4);
		puff *= octaves(vec2(residual_sun_aligment*9. + angle*1., angle/3.), vec2(12.*time, -time * 2.)*0.02, 7);
		puff *= octaves(vec2(-residual_sun_aligment*9. - angle*2., angle/4.), vec2(11.*time, -time * 1.7)*0.015, 4);
		puff *= octaves(vec2(-residual_sun_aligment*7. + angle*1., angle/2.), vec2(7.*time, -time * 1.9)*0.032, 4);
		puff = pow(puff, 2.) * 5.25;
		
		if(VCJDECLIPSESIZE <= 1.)
			sun = sun + sun * (1. - VCJDECLIPSESIZE) * aligment_angle/ (abs(height - light_height)) * (0.05 + puff);
		else
		{			
			if(residual_sun_aligment < 0.)
				sun = exp(22.*residual_sun_aligment) * (1. + puff);
			else sun = 0.2*exp(-25.*residual_sun_aligment) * (0.5 + puff);
		}
	}
	
	vec3 pyramid;
	
	if(ENABLE_DAMN_PYRAMID && VCJDECLIPSESIZE > 2.)
	{ // pyramid
		
		const float pyrasize = 0.11;
		const float p_z_cutoff = 0.10269;
		const float thickness = 7.;
		const float fadeout_flicker_time = 1e10;
		float fadeout_time = fadeout_flicker_time / (time + fadeout_flicker_time);
		
		float outtershell = 
			clamp(-abs(dirdiff.z - dirdiff.y + 0.7 * dirdiff.x + pyrasize)*thickness + 0.005, 0., 0.03) +
			0.9*clamp(-abs(dirdiff.z - dirdiff.y - 0.7 * dirdiff.x + pyrasize)*thickness + 0.004, 0., 0.03) + 
			0.15*clamp(-abs(dirdiff.z - dirdiff.y + 2.4 * dirdiff.x + pyrasize)*thickness + 0.010, 0., 0.03);
		
		float flare = 0.;
		if(SUN_SIZE < 0.)
		{
			float flicker = 
				abs(sin(time * 7.88) + cos(time * 11.44) - sin(time * 5.31)) * 
				fadeout_time + time / (time + fadeout_flicker_time);
			float flare_v = clamp(-abs(dirdiff.x)*0.8 + 0.005, 0., 0.09) * 200. * flicker;
			float flare_h = 5.0 * flicker / abs((height - p_z_cutoff)*19.3 + EPSI) * clamp(-pow(dirdiff.x, 2.)*(3. - flicker)*97. + 0.01, 0., 1.);
			
			const float special_dist_p = 0.5;
			float distance_from_pyramid_tip = pow(abs(dirdiff.x), special_dist_p) + pow(abs(height - p_z_cutoff), special_dist_p);
						
			float central_flare = exp(0.25 * flicker / distance_from_pyramid_tip) * 0.01;
			flare += 2. * pow(float(p.z >= p_z_cutoff) * exp(-(p.z - p_z_cutoff)*1.) * (flare_v), 4.);
			//flare += 0.25 * flare_h;
			flare += 1. * central_flare;
			flare = pow(abs(flare), 0.95) * 0.9;
			
			outtershell += 
				0.21*clamp(-abs(dirdiff.z - dirdiff.y - 2.4 * dirdiff.x + pyrasize)*thickness + 0.0024, 0., 0.03);
		}
		
		outtershell *= 
			pow((p.z < p_z_cutoff && p.z > 0.) ? 3.* p.z * (p_z_cutoff - p.z + 0.01) : 0., 1.0) * fadeout_time;
		flare *= fadeout_time * 2.;
		
		pyramid += vec3(3.2,2.2,1.9) * 96000. * outtershell + flare * vec3(1.5,1.9,1.9);
		sun += 96000. * outtershell + flare;
	}
	
	vec3 norhtern_lights;
	
	if(true || ENABLE_DAMN_PYRAMID && VCJDECLIPSESIZE > 2.)
	{ //damn northern lights
		float accurateProjectedHight = tan(p.z * PI_C / 2. - EPSI) / (PI_C / 2.);
		float angle = trueAngle(p.x, p.y);
		float height1 = accurateProjectedHight / 0.1 - 1.5 + 
			sin(27. * angle + time * 0.25) * 0.1 + 
			cos(angle * 21. - time * 0.34) * 0.07;
		float puff = octaves(vec2(height1/9., angle*3.), vec2(time, 0.)*0.2, 4) /
			max(abs(height1) + 0.5, 0.) * 0.3;
		
		float height2 = accurateProjectedHight / 0.07 - 1.7;
		float height4 =  accurateProjectedHight / 0.2 - 5.5 + 
			sin(19. * angle + time * 0.55) * 0.31 + 
			cos(angle * 34. - time * 0.74) * 0.17;		
		float height3 = (height1 + height2 + height4) * 0.5 + 0.18; 
				
		puff += octaves(vec2(height2/2. + angle/9., angle*3.), vec2(12.*time, -time * 2.)*0.02, 7) / max(abs(height2) + 0.5, 0.) * 0.5;
		puff += octaves(vec2(-height3/3. - angle/8., angle*4.), vec2(11.*time, -time * 1.7)*0.015, 4) / max(abs(height2) + 0.5, 0.) * 0.5 ;
		puff *= octaves(vec2(-height4/8. + angle*4., angle*2.), vec2(7.*time, -time * 1.9)*0.032, 4) / max(abs(height4)*0.5 + 0.9, 0.) * 2.5;
		
		puff += octaves(vec2(height2 + angle * 2.1, angle * 1.3)*0.2, vec2(0.02*time,0.), 4) * 1.7;

		puff = pow(abs(puff), 1.95) * 0.25 * (clamp(-pow(height3, 2.) + 4.5, 0., 1.) + clamp(abs(height4 - 2.) + 0.5, 0., 1.));
		puff *= max(accurateProjectedHight / (0.1 + accurateProjectedHight) * (0.1) / (0.1 + accurateProjectedHight) / 0.25, 0.);
		
		norhtern_lights = vec3(0.01, puff, 0.5) * puff;
	}
	
	float luminocity = 0.;
	if(light_height > -PI_C/8.){
		luminocity = (light_height + dawn_color_range);
		if(luminocity > 1.)
			luminocity = 1.;
	}
	
	vec3 sky_dawn_color = aligment * dawn_color * dawn_color_mul;
	vec3 color = sqrt(luminocity) * AMBIENT + sky_dawn_color + sun*sun_color + pyramid + norhtern_lights;
	if(sun>EPSI)
		return color;
	return normalize_color(color);
}

vec3 get_current_ambient() {
	return sky_color_in_direction(vec3(-LIGHT_DIR.z, LIGHT_DIR.x, LIGHT_DIR.y));
}

vec3 get_color_with_shadows(vec3 p){
	float dist;
	int iter_count;
	vec3 normal = normal(maps.collision_map, p);
	
	float normal_light_dot = dot(normal, LIGHT_DIR);
	if(normal_light_dot < EPSI)
		return color_map(p) * AMBIENT;
	
	max_approach_point _map = rayMarcher(maps.collision_map, p, LIGHT_DIR, dist, iter_count, true);
	float shade_marched = 0.;
	
	if(_map.distance_to_map > EPSI) {
		if(ENABLE_EXACT_CLOSEST_APPROACH_SEARCH) {
			float visible_sun_size = sqrt(SUN_SIZE); 
			float part_of_visible_sun = (_map.distance_to_map / _map.distance_from_origin) / tan(visible_sun_size);
			float clamped_visibility = clamp(part_of_visible_sun, 0., 1.);

			shade_marched = sqrt(clamped_visibility) * SUN_BRIGHTNESS;
		}
		else{
			shade_marched =  SUN_BRIGHTNESS;
		}
	}
	
	shade_marched *= pow(abs(dot(normal, LIGHT_DIR)),3.);
		
	if(false && _map._debug < -50.)
		return vec3(1.,0.,0.);
	
	vec3 end_pos = p + LIGHT_DIR*dist;
	
	vec3 color = color_map(p);
	
	if(collision_map(end_pos) > EPSI)
		return color * ((shade_marched) + AMBIENT);
	
	return color * AMBIENT;
}

vec4 reflected_raymarching(vec3 start_position, vec3 direction){
	vec4 color = vec4(0.);
	int iter_count;
	
	float accumulating_reflection_coef = 1.;
	for(int i = 0; i < MAX_REFLECTIONS; i++){
		float dist;
		iter_count=0;
		rayMarcher(maps.collision_map, start_position, direction, dist, iter_count, false);
		float shading = float(iter_count)/float(MAX_ITERS);
		vec3 endpoint = start_position + dist*direction;
		vec3 norm = normal(maps.collision_map, endpoint);
		float reflection_coef = general_map(maps.reflection_map, endpoint);
		vec4 cur_color;
		bool is_sky = length(start_position - endpoint) > float(MAX_DIST);
		
		vec4 visibility = pow(vec4(0.9, 0.8, 0.6, 1.), vec4(dist*RAYLIGH_SCATTER_COEF));
		
		if(!is_sky){
			cur_color = vec4(get_color_with_shadows(endpoint),1.0) * (1. - 0.0*shading);
				
			color += cur_color*accumulating_reflection_coef*(1. - reflection_coef)*visibility;
			accumulating_reflection_coef *= reflection_coef;
		}
		else{
			cur_color = vec4(sky_color_in_direction(endpoint - start_position),1.);
			color += cur_color*accumulating_reflection_coef*visibility;
		}
		
		if(reflection_coef < float(EPSI) || is_sky)
			break;
		
		direction = reflect(direction, norm);
		start_position = endpoint + direction*2.*EPSI;
		
	}
	return color;
}

void _mainImage(out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = (fragCoord-0.5*resolution.xy)/resolution.x * ZOOM;
	vec3 start_position = position;
	vec3 direction = normalize(vec3(uv.x,1.,uv.y)); 
	direction = rotateY(direction, -2.*PI_C*(mouse.y - 0.5));
	direction = rotateX(direction, 2.*PI_C*(mouse.x - 0.5));   
	
	fragColor = reflected_raymarching(position, direction);
}

vec4 postprocess(vec4 localColor, vec2 fragCoord)
{
    vec2 position = ( fragCoord / resolution.xy );
    // Replace u_resolution with resolution in pixel calculation
    vec2 pixel = 1./resolution;
    
    if(ENABLE_DUST)
    {
        float dust1 = pow(octaves(vec2(-time*0.009) + (position)*1.5, vec2((-time*0.15)), 4), 0.3)*0.1;
        float dust2 = pow(octaves(
            vec2(-time*0.07, time*0.02) + 
                vec2(position.x * 0.6, position.y * 0.7),
            vec2(-time*0.19, time*0.1),
            4), 0.3)*0.2;
        float dust3 = pow(octaves(
            vec2(time*0.09, -time*0.03) + 
                vec2(position.x * 0.4, -position.y * 0.5),
            vec2(-time*0.17, time*0.12),
            12), 0.3)*0.1;
        
        float dust = dust1 + dust2;
        localColor = max(localColor - 
                 vec4(dust1,dust1,dust1,0.) +
                 vec4(dust2,dust2,dust2,0.) -
                 vec4(dust3,dust3,dust3,0.),
            0.);
    }
    
    return localColor;
}

void main(void) {
    // Directly use the `time` and `resolution` uniforms in your logic.
    AMBIENT=get_current_ambient();
    vec4 outputVec;
    vec2 fragCoord = gl_FragCoord.xy;
    _mainImage(outputVec, fragCoord);
    outputVec = postprocess(outputVec, fragCoord);
    gl_FragColor = outputVec;
}
