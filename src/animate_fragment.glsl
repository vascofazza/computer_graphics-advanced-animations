#version 120

varying vec3 pos;                   // [from vertex shader] position in world space
varying vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
varying vec2 texcoord;              // [from vertex shader] texture coordinate

uniform vec3 camera_pos;            // camera position (center of the camera frame)

uniform vec3 ambient;               // scene ambient

uniform int lights_num;             // number of lights
uniform vec3 light_pos[16];         // light positions
uniform vec3 light_intensity[16];   // light intensities

uniform vec3 material_kd;           // material kd
uniform vec3 material_ks;           // material ks
uniform float material_n;           // material n
uniform bool material_is_lines;     // whether the material is lines or meshes
uniform bool material_double_sided; // whether the material is double sided

uniform bool material_kd_txt_on;    // material kd texture enabled
uniform sampler2D material_kd_txt;  // material kd texture
uniform bool material_ks_txt_on;    // material ks texture enabled
uniform sampler2D material_ks_txt;  // material ks texture
uniform bool material_norm_txt_on;    // material norm texture enabled
uniform sampler2D material_norm_txt;  // material norm texture

// main
void main() {
    // re-normalize normals
    vec3 n = normalize(norm);
    // lookup normal map if needed
    if(material_norm_txt_on) n = normalize(2*texture2D(material_norm_txt,texcoord).xyz-vec3(1));
    // check for double sided
    //FIXED
    if(material_double_sided) n = faceforward(n,normalize(pos-camera_pos),n);

    
    // compute material values by looking up textures is necessary
    vec3 kd = material_kd * ( (material_kd_txt_on)?texture2D(material_kd_txt,texcoord).xyz:vec3(1) );
    vec3 ks = material_ks * ( (material_ks_txt_on)?texture2D(material_ks_txt,texcoord).xyz:vec3(1) );
    // accumulate ambient
    vec3 c = ambient * kd;
    // foreach light
    for(int i = 0; i < lights_num; i ++) {
        // compute point light color at pos
        vec3 cl = light_intensity[i] / pow(length(light_pos[i]-pos),2);
        // compute light direction at pos
        vec3 l = normalize(light_pos[i]-pos);
        // compute view direction using camera_pos and pos
        vec3 v = normalize(camera_pos-pos);
        // compute h
        vec3 h = normalize(v+l);
        // accumulate blinn-phong model
        if(material_is_lines) {
            c += cl * kd * sqrt(1-dot(l,n)*dot(l,n));
        } else {
            c += cl * max(0,dot(l,n)) * (kd + ks * pow(max(0,dot(h,n)),material_n));
        }
        //c = vec3(c[0]*(1-texcoord.x), c[1]*(1-texcoord.y), c[2]);
    }
    // output final color by setting gl_FragColor
    gl_FragColor = vec4(c,1);
}
