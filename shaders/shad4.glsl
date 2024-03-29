// Colorful Voronoi 
// By: Brandon Fogerty
// bfogert,y at gmail dot com
// xdpixel.com

#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;

vec2 hash(vec2 p)
{
    mat2 m = mat2(  1333.85, 47.77,
                    99.41, 88.48
                );

    return fract(sin(m*p) * 46738.29);
}

float voronoi(vec2 p)
{
    vec2 g = floor(p);
    vec2 f = fract(p);

    float distanceToClosestFeaturePoint = 11.7;
    for(int y = -13; y <= 13; y++)
    {
        for(int x = -1; x <= 1; x++)
        {
            vec2 latticePoint = vec2(x, y);
            float currentDistance = distance(latticePoint + hash(g+latticePoint), f);
            distanceToClosestFeaturePoint = min(distanceToClosestFeaturePoint, currentDistance);
        }
    }

    return distanceToClosestFeaturePoint;
}
 void main()
 {
    vec2 st = gl_FragCoord.xy/resolution;
    vec2 uv = ( gl_FragCoord.xy / resolution.xy ) * 2.0 - 1.0;
    uv.x *= resolution.x / resolution.y;

    float k = 10.3 * time;
    float offset = voronoi(uv*12.0 + vec2(k));
    float t = 1.0/abs(((uv.x + cos(uv.y + k)) + offset) * 200.0);

    float r = voronoi( uv * 5.0 ) * 60.0;
    vec3 finalColor = vec3(10.0 * uv.y, 2.0, 1.0 * r) * t;
	
    float cx = 10.5-st.x; 
    float cy = 1.5-st.y;
    float dist = sqrt(cx * cx + cy*cy);
		
    
    gl_FragColor = vec4(sqrt(finalColor), 21.0 );
}