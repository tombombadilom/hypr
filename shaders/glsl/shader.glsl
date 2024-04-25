// Vertex Shader Code (as a string)
const GLchar* vertexShaderSource = "#version 300 es\n"
		"precision mediump float;\n"
		"in vec3 a_position;\n"
		"in vec2 a_texCoord;\n"
		"out vec2 v_texCoord;\n"
		"void main() {\n"
		"   gl_Position = vec4(a_position, 1.0);\n"
		"   v_texCoord = a_texCoord;\n"
		"}\n";

// Fragment Shader Code (as a string)
const GLchar* fragmentShaderSource = "#version 300 es\n"
		"precision mediump float;\n"
		"in vec2 v_texCoord;\n"
		"uniform sampler2D u_texture;\n"
		"uniform float u_time;\n"
		"out vec4 fragColor;\n"
		"void main() {\n"
		"   vec4 texColor = texture(u_texture, v_texCoord);\n"
		"   float red = texColor.r + sin(u_time);\n"
		"   float green = texColor.g + cos(u_time);\n"
		"   float blue = texColor.b;\n"
		"   fragColor = vec4(red, green, blue, 1.0);\n"
		"}\n";

// Create and compile vertex shader
GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
glCompileShader(vertexShader);

// Check for vertex shader compile errors...
GLint vertexCompiled;
glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &vertexCompiled);
if (vertexCompiled != GL_TRUE) {
		GLsizei log_length = 0;
		GLchar message[1024];
		glGetShaderInfoLog(vertexShader, 1024, &log_length, message);
		// Write the error to a log or std::cerr
}

// Create and compile fragment shader
GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
glCompileShader(fragmentShader);

// Check for fragment shader compile errors...
GLint fragmentCompiled;
glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &fragmentCompiled);
if (fragmentCompiled != GL_TRUE) {
		GLsizei log_length = 0;
		GLchar message[1024];
		glGetShaderInfoLog(fragmentShader, 1024, &log_length, message);
		// Write the error to a log or std::cerr
}

// Create program, attach shaders, link and use...
GLuint shaderProgram = glCreateProgram();
glAttachShader(shaderProgram, vertexShader);
glAttachShader(shaderProgram, fragmentShader);
glLinkProgram(shaderProgram);

// Check for linking errors...
GLint programLinked;
glGetProgramiv(shaderProgram, GL_LINK_STATUS, &programLinked);
if (programLinked != GL_TRUE) {
		GLsizei log_length = 0;
		GLchar message[1024];
		glGetProgramInfoLog(shaderProgram, 1024, &log_length, message);
		// Write the error to a log or std::cerr
}

glUseProgram(shaderProgram);

// Clean up shaders since they're linked into our program now and no longer necessary
glDeleteShader(vertexShader);
glDeleteShader(fragmentShader);
