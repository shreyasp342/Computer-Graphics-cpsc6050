// main.frag ..............................
#define BLEND_REPLACE 1
#define BLEND_NONE 2

#define M_PI 3.141592653589793238462643383279

varying vec3 ec_vnormal, ec_vposition;

uniform vec4 diffuse, specular;
uniform float shininess;
uniform int numLights;

void main() {

  vec3 P, N, V, L, H;
  vec4 diffuse_m = diffuse;
  vec4 specular_m = specular;
  vec4 diffuse_color, specular_color;

  P = ec_vposition;
  N = normalize(ec_vnormal);
  V = normalize(-P);

  L = normalize(gl_LightSource[0].position.xyz - P);
  H = normalize(L+V);
    
  diffuse_color = gl_LightSource[0].diffuse*max(dot(N,L),0.0)*diffuse_m; 
  specular_color = gl_LightSource[0].specular*(((shininess + 2.0) / 8.0 * M_PI)
	  *pow(max(dot(H,N),0.0),shininess))*specular_m; 
	 
	float d = distance(gl_LightSource[0].position.xyz, P);
	float attenuation = 1.0 / (gl_LightSource[0].constantAttenuation + gl_LightSource[0].linearAttenuation * d + gl_LightSource[0].quadraticAttenuation * d * d);  

  gl_FragColor = (diffuse_color + specular_color) * attenuation; 
  gl_FragColor = vec4(gl_FragColor.xyz, 1.0);
}
