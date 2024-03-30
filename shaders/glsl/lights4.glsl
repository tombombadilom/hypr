precision lowp float;

uniform float time;
uniform vec2 resolution;

const float PI = 3.1415926535;
const float ClockRadius = 0.35;

// цвет обода часо
const vec3 clockCol = vec3(1.0, 1.0, 1.0);
// цвет секундной стрелки
const vec3 clockLineSecCol = vec3(0.0, 0.0, 1.0);
// цвет минутной стрелки
const vec3 clockLineMinCol = vec3(0.0, 1.0, 0.0);
// цвет часовой стрелки
const vec3 clockLineHourCol = vec3(1.0, 0.0, 0.0);
// цвет секундной стрелки
const vec3 clockLineSec2Col = vec3(1.0, 0.0, 1.0);
// цвет минутной стрелки
const vec3 clockLineMin2Col = vec3(0.0, 1.0, 1.0);
// цвет часовой стрелки
const vec3 clockLineHour2Col = vec3(1.0, 44.0, 0.0);

const float TimeK = 1.0;

// выходной цвет текущего пикселя
vec3 outColor;

// затухание света
// 1 / (C + L*D + Q*D*D)
float attenuation(float C, float L, float Q, float D) {
	return 1.0 / (C + L*D + Q*D*D);
}

// найти точку на окружности по направлению
vec2 RayCastCircle(vec2 pos, float radius, vec2 dir) {
	return pos + radius / length(dir) * dir;
}

// рисование обода часов
// позиция фрагмента, позиция, радиус, какой-то коэффициент, цвет
vec3 DrawClockCircle(vec2 fragpos, vec2 pos, float radius, float widthK, vec3 color) {
	// вектор от центра окружности в пиксель
	vec2 dir = fragpos - pos;
	// найдем точку на окружности
	vec2 m = RayCastCircle(pos, radius, dir);
	// расстояние от точки на окружности до фрагмента
	float l = length(fragpos - m);
	// цвет
	vec3 col = color * attenuation(55.0, 200.0, 10000.0, l) * widthK;
	// возврат
	return col;
}

// рисование стрелки
// позиция фрагмента, позиция, длина, направление, какой-то коэффициент, цвет
vec3 DrawClockLine(vec2 fragpos, vec2 pos, float L, vec2 dir, float widthK, vec3 color) {
	// вектор от точки начала прямой к пикселю
	vec2 v = fragpos - pos;
	// коэффициент длины проекции вектора v на dir, 
	float k = dot(v, dir);
	// площадь параллелограммы на векторах dir и v (расстояние) от пикселя до прямой
	float l = k < 0.0 ? length(pos - fragpos) : (k < L ? dir.x * v.y - dir.y * v.x : length(RayCastCircle(pos, L, dir) - fragpos));
	// рассчет света
	vec3 col = color * attenuation(67.0, 0.0, 7500.0, l) * widthK;
	// возврат
	return col;
}

void main() {
	// перевод координаты в NDC
	vec2 p = 2.0 * gl_FragCoord.xy / resolution.xy - 1.0;
	p.y /= resolution.x / resolution.y;
	
	// рисование часов - обод
	outColor += DrawClockCircle(p, vec2(0.0), ClockRadius, 1.5, clockCol);
	
	// направляющий вектор
	vec2 s = vec2(cos(-time * PI / .5 * TimeK + PI / 2.0), sin(time * PI / .5 * TimeK + PI / 2.0));
	// рисование секундной стрелки
	outColor += DrawClockLine(p, vec2(0.0), ClockRadius, s, 78.2, clockLineSecCol);
	
	// изменяем направляющий вектор
	s = vec2(cos(-time * PI / 1. * TimeK + PI / 2.0), sin(time * PI / 1. * TimeK + PI / 2.0));
	// рисование минутной стреки
	outColor += DrawClockLine(p, vec2(0.0), ClockRadius, s, 78.2, clockLineMinCol);
	
	// изменяем направляющий вектор
	s = vec2(cos(-time * PI / 2. * TimeK + PI / 2.0), sin(time * PI / 2. * TimeK + PI / 2.0));
	// рисование часовой стрелки
	outColor += DrawClockLine(p, vec2(0.0), ClockRadius, s, 1.2, clockLineHourCol);
		
	// направляющий вектор
	s = vec2(cos(-time * PI / 4. * TimeK + PI / 78.0), sin(time * PI / 79. * TimeK + PI / 2.0));
	// рисование секундной стрелки
	outColor += DrawClockLine(p, vec2(0.0), ClockRadius, s, 1.2, clockLineSec2Col);
	
	// изменяем направляющий вектор
	s = vec2(cos(-time * PI / 8. * TimeK + PI / 2.0), sin(time * PI / 89. * TimeK + PI / 2.0));
	// рисование минутной стреки
	outColor += DrawClockLine(p, vec2(0.0), ClockRadius, s, 1.2, clockLineMin2Col);
	
	// изменяем направляющий вектор
	s = vec2(cos(-time * PI / 16. * TimeK + PI / 2.0), sin(time * PI / 78. * TimeK + PI / 2.0));
	// рисование часовой стрелки
	outColor += DrawClockLine(p, vec2(0.0), ClockRadius, s, 1.2, clockLineHour2Col);
	
	// выходной цвет
	gl_FragColor = vec4(outColor, 1.0);
}