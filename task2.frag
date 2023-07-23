//============================================================================
// PROJECT ID: SWS3005_12
//
// GROUP NUMBER: 12
//
// STUDENT NAME: HU ZHANPENG
// NUS User ID.: t0930152
//
// STUDENT NAME: NING JUNTING
// NUS User ID.: t0930085
//
// STUDENT NAME: LUO TANGWEN
// NUS User ID.: t0930258
//
// COMMENTS TO GRADER:
//
//============================================================================
const int SAMPLES_PER_PIXEL = 10;
const float PI = 3.1415926536;
const int BOUNCES = 3;
const float P_TERMINATE = 0.75;
const float EPSILON = 0.01;


#define saturate(x) clamp(x, 0.0, 1.0)
// sRGB, linear space conversions
#define stol1(x) (x <= 0.04045 ? x / 12.92 : pow((x + 0.055) / 1.055, 2.4))
#define stol3(x, y, z) vec3(stol1(x), stol1(y), stol1(z))
#define ltos1(x) (x <= 0.0031308 ? x * 12.92 : 1.055 * pow(x, 0.4166667) - 0.055)
#define ltos3(x, y, z) vec3(ltos1(x), ltos1(y), ltos1(z))


const vec3 BACKGROUND_COLOR = vec3( 0.1, 0.2, 0.6 );

// Vertical field-of-view angle of camera. In radians.
const float FOVY = 50.0 * PI / 180.0;

// Use this for avoiding the "epsilon problem" or the shadow acne problem.
const float DEFAULT_TMIN = 10.0e-4;

// Use this for tmax for non-shadow ray intersection test.
const float DEFAULT_TMAX = 10.0e6;


// Constants for the scene objects.
const int NUM_LIGHTS = 3;
const int NUM_MATERIALS = 17;
const int NUM_SPHERES = 22;
const int NUM_CUBES = 17;
const int NUM_PLANES = 5;


//============================================================================
// Define new struct types.
//============================================================================
struct Ray_t {
    vec3 o;  // Ray Origin.
    vec3 d;  // Ray Direction. A unit vector.
    float t; // Length of Ray.
};

struct Camera {
    vec3 pos;
    vec3 lookat;
    vec3 axis[3]; // x,y,z
};

struct Material_t {
    vec3 albedo;
    float roughness;
    float metalness;
    vec3 emission;
};

struct Hit {
    bool hit;
    vec3 position;
    vec3 normal;
    int MaterialID;
};

struct Plane_t {
    // The plane equation is Ax + By + Cz + D = 0.
    float A, B, C, D;
    int materialID;
};

struct Sphere_t {
    vec3 center;
    float radius;
    int materialID;
};


struct Cube_t {
    vec3 center;
    int materialID;
    vec3 size;
};


struct Light_t {
    vec3 position;
    vec3 E;
};
//============================================================================
// Global scene data.
//============================================================================

Plane_t Plane[NUM_PLANES];
Sphere_t Sphere[NUM_SPHERES];
Cube_t Cube[NUM_CUBES];
Light_t Light[NUM_LIGHTS];
Material_t Material[NUM_MATERIALS];
Camera camera;

void InitPlane() {
    
    // Bottom Plane
    Plane[0].A = 0.0;
    Plane[0].B = 1.0;
    Plane[0].C = 0.0;
    Plane[0].D = 9.65;
    Plane[0].materialID = 12;
    
    // Left Plane
    Plane[1].A = -1.0;
    Plane[1].B = 0.0;
    Plane[1].C = 0.0;
    Plane[1].D = 30.0;
    Plane[1].materialID = 10;
    
    // Right Plane
    Plane[2].A = 1.0;
    Plane[2].B = 0.0;
    Plane[2].C = 0.0;
    Plane[2].D = 30.0;
    Plane[2].materialID = 10;
    
    // Front Plane
    Plane[3].A = 0.0;
    Plane[3].B = 0.0;
    Plane[3].C = -1.0;
    Plane[3].D = 25.0;
    Plane[3].materialID = 10;
    
    // Back Plane
    Plane[4].A = 0.0;
    Plane[4].B = 0.0;
    Plane[4].C = 1.0;
    Plane[4].D = 25.0;
    Plane[4].materialID = 10;
}

void InitBall() {
    float timeBias = 0.06;
    
    // Black Ball
    Sphere[0].center = vec3(-15.0, 0.5 + 2.0 * abs(cos(5.0 * iTime)), 0.0);
    Sphere[0].radius = 0.5;
    Sphere[0].materialID = 0;
    

    // Red Ball * 15
    for (int i = 1; i <= 15; i++) {
        Sphere[i].radius = 0.5;
        Sphere[i].materialID = 1;
    }
    float time = iTime;
    Sphere[1].center = vec3(-12.0, 0.5 + 2.0 * abs(cos(5.0 * (iTime + timeBias))), 2.0);
    
    Sphere[2].center = vec3(-12.0, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 2.0 * timeBias))), 1.0);
    
    Sphere[3].center = vec3(-12.0, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 3.0 * timeBias))), 0.0);
    
    Sphere[4].center = vec3(-12.0, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 4.0 * timeBias))), -1.0);
    
    Sphere[5].center = vec3(-12.0, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 5.0 * timeBias))), -2.0);
    
    Sphere[6].center = vec3(-11.13, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 6.0 * timeBias))), 1.5);
    
    Sphere[7].center = vec3(-11.14, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 7.0 * timeBias))), 0.5);
    
    Sphere[8].center = vec3(-11.14, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 8.0 * timeBias))), -0.5);
    
    Sphere[9].center = vec3(-11.14, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 9.0 * timeBias))), -1.5);
    
    Sphere[10].center = vec3(-10.26, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 10.0 * timeBias))), 1.0);
    
    Sphere[11].center = vec3(-10.28, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 11.0 * timeBias))), 0.0);
    
    Sphere[12].center = vec3(-10.28, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 12.0 * timeBias))), -1.0);
    
    Sphere[13].center = vec3(-9.39, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 13.0 * timeBias))), 0.5);
    
    Sphere[14].center = vec3(-9.39, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 14.0 * timeBias))), -0.5);
    
    Sphere[15].center = vec3(-8.52, 0.5 + 2.0 * abs(cos(5.0 * (iTime + 15.0 * timeBias))), 0.0);
    

    // Purple Ball
    Sphere[16].center = vec3(-7.52, 0.5 + 2.0 * abs(sin(5.0 * iTime)), 0.0);
    Sphere[16].radius = 0.5;
    Sphere[16].materialID = 2;
    
    
    // Blue Ball
    Sphere[17].center = vec3(0.0, 0.5, 0.0);
    Sphere[17].radius = 0.5;
    Sphere[17].materialID = 3;
    

    // Green Ball
    Sphere[18].center = vec3(11.0, 0.5, 2.5);
    Sphere[18].radius = 0.5;
    Sphere[18].materialID = 4;
    
    
    // Brown Ball
    Sphere[19].center = vec3(11.0, 0.5, 0.0);
    Sphere[19].radius = 0.5;
    Sphere[19].materialID = 5;
    

    // Yellow Ball
    Sphere[20].center = vec3(11.0, 0.5, -2.5);
    Sphere[20].radius = 0.5;
    Sphere[20].materialID = 6;
    

    // White Ball
    Sphere[21].center = vec3(-4.0, 0.5, 6.0);
    Sphere[21].radius = 0.5;
    Sphere[21].materialID = 7;

}

void InitCube() {
    // Back Left Baffle
    Cube[0].center = vec3(-8.75, 0.38, -9.6);
    Cube[0].size = vec3(16.3, 0.76, 1.2);
    Cube[0].materialID = 9;
    
    // Back Right Baffle
    Cube[1].center = vec3(8.75, 0.38, -9.6);
    Cube[1].size = vec3(16.3, 0.76, 1.2);
    Cube[1].materialID = 9;
    
    // Front Left Baffle
    Cube[2].center = vec3(-8.75, 0.38, 9.6);
    Cube[2].size = vec3(16.3, 0.76, 1.2);
    Cube[2].materialID = 9;
    
    // Front Right Baffle
    Cube[3].center = vec3(8.75, 0.38, 9.6);
    Cube[3].size = vec3(16.3, 0.76, 1.2);
    Cube[3].materialID = 9;
    
    // Left Baffle
    Cube[4].center = vec3(-18.6, 0.38, 0.0);
    Cube[4].size = vec3(1.2, 0.76, 15.8);
    Cube[4].materialID = 9;
    
    // Right Baffle
    Cube[5].center = vec3(18.6, 0.38, 0.0);
    Cube[5].size = vec3(1.2, 0.76, 15.8);
    Cube[5].materialID = 9;
    
    // Table Plane
    Cube[6].center = vec3(0.0, -0.05, 0.0);
    Cube[6].size = vec3(36.0, 0.1, 18.0);
    Cube[6].materialID = 8;
    
    // Hole Plane
    Cube[7].center = vec3(0.0, -0.15, 0.0);
    Cube[7].size = vec3(38.4, 0.1, 20.4);
    Cube[7].materialID = 13;
    
    // Second Plane
    Cube[8].center = vec3(0.0, -1.7, 0.0);
    Cube[8].size = vec3(39.2, 3.0, 21.2);
    Cube[8].materialID = 14;
    
    // Front Fence
    Cube[9].center = vec3(0.0, 0.28, 10.4);
    Cube[9].size = vec3(39.2, 0.96, 0.4);
    Cube[9].materialID = 9;
    
    // Back Fence
    Cube[10].center = vec3(0.0, 0.28, -10.4);
    Cube[10].size = vec3(39.2, 0.96, 0.4);
    Cube[10].materialID = 9;
    
    // Left Fence
    Cube[11].center = vec3(-19.4, 0.28, 0.0);
    Cube[11].size = vec3(0.4, 0.96, 20.4);
    Cube[11].materialID = 9;
    
    // Right Fence
    Cube[12].center = vec3(19.4, 0.28, 0.0);
    Cube[12].size = vec3(0.4, 0.96, 20.4);
    Cube[12].materialID = 9;
    
    // First Leg
    Cube[13].center = vec3(-18.15, -6.4, 9.15);
    Cube[13].size = vec3(2.5, 6.4, 2.5);
    Cube[13].materialID = 15;
    
    // Second Leg
    Cube[14].center = vec3(18.15, -6.4, 9.15);
    Cube[14].size = vec3(2.5, 6.4, 2.5);
    Cube[14].materialID = 15;
    
    // Third Leg
    Cube[15].center = vec3(-18.15, -6.4, -9.15);
    Cube[15].size = vec3(2.5, 6.4, 2.5);
    Cube[15].materialID = 15;
    
    // Fourth Leg
    Cube[16].center = vec3(18.15, -6.4, -9.15);
    Cube[16].size = vec3(2.5, 6.4, 2.5);
    Cube[16].materialID = 15;
}

void InitLight() {
    // Light 0
    Light[0].position = vec3(9.0 * cos(iTime), 15.0, 9.0 * sin(iTime));
    Light[0].E = vec3(1, 1, 1) * vec3(1500);
    
    // Light 1
    Light[1].position = vec3(-29.0, 12.0, 19.0);
    Light[1].E = vec3(1, 1, 1) * vec3(2000);
    // Light 2
    Light[2].position = vec3(-29.0, 12.0, 19.0);
    Light[2].E = vec3(1, 1, 1) * vec3(2000);


}
const vec2 s = vec2(1, 1.7320508);

float hex(in vec2 p){
    p = abs(p);
    return max(dot(p, s*.5), p.x);
}

vec4 getHex(vec2 p){
    vec4 hC = floor(vec4(p, p - vec2(.5, 1))/s.xyxy) + .5;
    vec4 h = vec4(p - hC.xy*s, p - (hC.zw + .5)*s);
    return dot(h.xy, h.xy)<dot(h.zw, h.zw) ? vec4(h.xy, hC.xy) : vec4(h.zw, hC.zw + vec2(.5, 1)); 
}

float aafract(float x) {
    float v = fract(x),
          w = fwidth(x);
    return v < 1.-w ? v/(1.-w) : (1.-v)/w;
}
void InitMaterial() {
    
    // Black Plastic Material
    Material[0].albedo = vec3(0, 0, 0);
    Material[0].roughness = 0.5;
    Material[0].metalness = 0.6;
    Material[0].emission = vec3(0.0);
    // Red Plastic Material
    Material[1].albedo = vec3(0.95, 0.09, 0.05);
    Material[1].roughness = 0.5;
    Material[1].metalness = 0.6;
    Material[1].emission = vec3(0.0);
    // Purple Plastic Material
    Material[2].albedo = vec3(0.86, 0.44, 0.84);
    Material[2].roughness = 0.5;
    Material[2].metalness = 0.6;
    Material[2].emission = vec3(0.0);
    // Blue Plastic Material
    Material[3].albedo = vec3(0.0, 0.0, 1.0);
    Material[3].roughness = 0.5;
    Material[3].metalness = 0.6;
    Material[3].emission = vec3(0.0);
    // Green Plastic Material
    Material[4].albedo = vec3(0.13, 0.55, 0.13);
    Material[4].roughness = 0.5;
    Material[4].metalness = 0.6;
    Material[4].emission = vec3(0.0);
    // Brown Plastic Material
    Material[5].albedo = vec3(0.78, 0.38, 0.08);
    Material[5].roughness = 0.5;
    Material[5].metalness = 0.6;
    Material[5].emission = vec3(0.0);
    // Yellow Plastic Material
    Material[6].albedo = vec3(1.0, 1.0, 0.0);
    Material[6].roughness = 0.5;
    Material[6].metalness = 0.6;
    Material[6].emission = vec3(0.0);
    // White Plastic Material
    Material[7].albedo = vec3(1.0, 1.0, 1.0);
    Material[7].roughness = 0.5;
    Material[7].metalness = 0.6;
    Material[7].emission = vec3(0.0);
    // Table Material
    Material[8].albedo = vec3(0.02, 0.36, 0.04);
    Material[8].roughness = 0.7;
    Material[8].metalness = 0.2;
    Material[8].emission = vec3(0.0);
    // Baffle Material
    Material[9].albedo = vec3(0.2, 0.05, 0.04);
    Material[9].roughness = 0.8;
    Material[9].metalness = 0.2;
    Material[9].emission = vec3(0.0);
    // Room Material
    Material[10].albedo = vec3(1, 1, 1);
    Material[10].roughness = 1.0;
    Material[10].metalness = 0.2;
    Material[10].emission = vec3(0.0);
    // NUS Material
    Material[11].albedo = vec3(0.0, 0.0, 0.6);
    Material[11].roughness = 0.8;
    Material[11].metalness = 0.2;
    Material[11].emission = vec3(0.4);
    // Hole Material
    Material[13].albedo = vec3(0.02, 0.02, 0.02);
    Material[13].roughness = 1.0;
    Material[13].metalness = 0.0;
    Material[13].emission = vec3(0.0);
    // Second Plane Material
    Material[14].albedo = vec3(0.95, 0.90, 0.2);
    Material[14].roughness = 0.5;
    Material[14].metalness = 0.6;
    Material[14].emission = vec3(0.0);
    // Leg Material
    Material[15].albedo = vec3(0.46, 0.63, 0.69);
    Material[15].roughness = 0.5;
    Material[15].metalness = 0.5;
    Material[15].emission = vec3(0.0);

}

void getFloorMaterial(vec2 xz) {
    // Floor Material
    xz = vec2((xz.x+30.0)/60.0, (xz.y+25.0)/50.0);
    vec2 uv = 20.* xz / s;
    vec4 h1 = getHex(uv);
    vec4 h2 = getHex(uv - 1./s);
    vec4 h3 = getHex(uv + 1./s);
    
    float v1 = aafract(hex(h1.xy)/.2);
    float v2 = aafract(hex(1.5*h2.xy)/0.45);
    float v3 = aafract(hex(2.*h3.xy)/0.3);
    
    Material[12].albedo = vec3(vec3(v1+v2+v3).r/8., 0, 0);
    Material[12].roughness = 0.8;
    Material[12].metalness = 0.8;
    Material[12].emission = vec3(0.0);
}

float smin(float a,float b,float k){
    float h = clamp(0.5+0.5*(a-b)/k,0.0,1.0);
    return mix(a,b,h)-k*h*(1.0-h);
}
float smax(float a,float b,float k){
    return -smin(-a,-b,k);
}
//Signed Distance Function of segment shape
float udSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 ba = b-a;
    vec2 pa = p-a;
    float h = clamp(dot(pa,ba)/dot(ba,ba), 0.0, 1.0);
    return length(pa-h*ba);
}

//Signed Distance Function of horsehoe shape
float sdHorseshoe(in vec2 p, in vec2 c, in float r, in vec2 w) {
    p.x = abs(p.x);
    float l = length(p);
    p = mat2(-c.x, c.y, 
              c.y, c.x)*p;
    p = vec2((p.y>0.0 || p.x>0.0)?p.x:l*sign(-c.x),
             (p.x>0.0)?p.y:l);
    p = vec2(p.x,abs(p.y-r))-w;
    return length(max(p, 0.0)) + min(0.0, max(p.x, p.y));
}

//Signed Distance Function of Letter 'N'
float N_sdf(in vec2 p) {
    vec2 offset = vec2(-2.0, 0.0);
    float d = udSegment(p-offset, vec2(-0.6, 2.0), vec2(-0.6, 4.0)) - 0.2;
    float d1 = udSegment(p-offset, vec2(-0.6, 4.0), vec2(0.6, 2.0)) - 0.2;
    float d2 = udSegment(p-offset, vec2(0.6, 2.0), vec2(0.6, 4.0)) - 0.2;
    return smin(d2,smin(d, d1, 0.04), 0.04); 
}

//Signed Distance Function of Letter 'U'
float U_sdf(in vec2 p) {
    float d = udSegment(p, vec2(-0.6, 2.6), vec2(-0.6, 4.0)) - 0.2;
    float d1 = sdHorseshoe(p-vec2(0.0, 2.6),vec2(cos(1.6), sin(1.6)), 0.6, vec2(0.2,0.2));
    float d2 = udSegment(p, vec2(0.6, 2.6), vec2(0.6, 4.0)) - 0.2;
    return smin(d2,smin(d, d1, 0.04), 0.04);
}

////Signed Distance Function of Letter 'S'
float S_sdf(in vec2 p) {
    vec2 offset = vec2(2.0, 0.0);
    float d1 = udSegment(p-offset, vec2(-0.6, 2.0), vec2(0.2, 2.0)) - 0.2;
    float d2 = udSegment(p-offset, vec2(-0.2, 4.0), vec2(0.6, 4.0)) - 0.2;
    float d3 = udSegment(p-offset, vec2(-0.2, 3.0), vec2(0.2, 3.0)) - 0.2;
    float d4 = sdHorseshoe(p.yx-offset.yx-vec2(3.5, -0.2),vec2(cos(1.6),sin(1.6)), 0.5, vec2(0.2,0.2));
    float d5 = sdHorseshoe(-p.yx+offset.yx-vec2(-2.5, -0.2),vec2(cos(1.6),sin(1.6)), 0.5, vec2(0.2,0.2));
    
    return smin(d1,smin(d2,smin(d3,smin(d4,d5, 0.04), 0.04), 0.04), 0.04);
}

bool IntersectNUS(in vec2 p) {
    return N_sdf(p)<0.0||U_sdf(p)<0.0||S_sdf(p)<0.0;
}
void InitScene() {
 
    InitPlane();
    InitBall();
    InitCube();
    InitLight();
    InitMaterial();
}
/////////////////////////////////////////////////////////////////////////////
// Initializes the scene.
/////////////////////////////////////////////////////////////////////////////


// Generate basis matrix for given normal
mat3 formBasis(vec3 n) {
    // Make vector q that is non-parallel to n
    vec3 q = n;
    vec3 aq = abs(q);
    if (aq.x <= aq.y && aq.x <= aq.z) {
        q.x = 1.f;
    } else if (aq.y <= aq.x && aq.y <= aq.z) {
        q.y = 1.f;
    } else {
        q.z = 1.f;
    }

    // Generate two vectors perpendicular to n
    vec3 t = normalize(cross(q, n));
    vec3 b = normalize(cross(n, t));

    // Construct the rotation matrix
    mat3 m;
    m[0] = t;
    m[1] = b;
    m[2] = n;
    return m;
}


// --SHADING-------------------------------------------------------------------
// Lambert diffuse term
vec3 lambertBRFD(vec3 albedo) {
    return albedo / PI;
}

// GGX distribution function
float ggx(float NoH, float roughness) {
    float a2 = roughness * roughness;
    a2 *= a2;
    float denom = NoH * NoH * (a2 - 1.0) + 1.0;
    return a2 / (PI * denom * denom);
}

// Schlick fresnel function
vec3 schlickFresnel(float VoH, vec3 f0) {
    return f0 + (1.0 - f0) * pow(1.0 - VoH, 5.0);
}

// Schlick-GGX geometry function
float schlick_ggx(float NoL, float NoV, float roughness) {
    float k = roughness + 1.0;
    k *= k * 0.125;
    float gl = NoL / (NoL * (1.0 - k) + k);
    float gv = NoV / (NoV * (1.0 - k) + k);
    return gl * gv;
}

// Evaluate the Cook-Torrance specular BRDF
vec3 cookTorranceBRDF(float NoL, float NoV, float NoH, float VoH, vec3 F, float roughness) {
    vec3 DFG = ggx(NoH, roughness) * F * schlick_ggx(NoL, NoV, roughness);
    float denom = 4.0 * NoL * NoV + 0.0001;
    return DFG / denom;
}

// Evaluate combined diffuse and specular BRDF
vec3 evalBRDF(vec3 n, vec3 v, vec3 l, Material_t m) {
    // Common dot products
    float NoV = saturate(dot(n, v));
    float NoL = saturate(dot(n, l));
    vec3 h = normalize(v + l);
    float NoH = saturate(dot(n, h));
    float VoH = saturate(dot(v, h));

    // Use standard approximation of default fresnel
    vec3 f0 = mix(vec3(0.04), m.albedo, m.metalness);
    vec3 F = schlickFresnel(VoH, f0);

    // Diffuse amount
    vec3 Kd = (1.0 - F) * (1.0 - m.metalness);

    return (Kd * lambertBRFD(m.albedo) + cookTorranceBRDF(NoL, NoV, NoH, VoH, F, m.roughness)) * NoL;
}
//-------------------------------------------------------------------------

//---Intersection-----------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is such an intersection, outputs the value of t, the position
// of the intersection (hitPos) and the normal vector at the intersection
// (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane(in Plane_t pln, in Ray_t ray, in float tmin, in float tmax,
                     out float t, out vec3 hitPos, out vec3 hitNormal) {
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;

    // We have a hit -- output results.
    t = t0;
    hitPos = ray.o + t0 * ray.d;
    hitNormal = normalize( N );
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane(in Plane_t pln, in Ray_t ray, in float tmin, in float tmax) {
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is one or two such intersections, outputs the value of the
// smaller t, the position of the intersection (hitPos) and the normal
// vector at the intersection (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere(in Sphere_t sph, in Ray_t ray, in float tmin, in float tmax,
                      out float t, out vec3 hitPos, out vec3 hitNormal) {
    /////////////////////////////////
    // TASK: WRITE YOUR CODE HERE. //
    /////////////////////////////////
    float a = dot(ray.d, ray.d);
    float b = 2.0 * dot(ray.o - sph.center, ray.d);
    float c = dot(ray.o - sph.center, ray.o - sph.center) - sph.radius * sph.radius;
    float det_square = b * b - 4.0 * a * c;
    if (det_square < 0.0) return false;
    float det = sqrt(det_square);
    float t1 = (-b - det) / (2.0 * a);
    float t2 = (-b + det) / (2.0 * a);
    
    if (tmin <= t1 && t1 <= tmax) t = t1;
    else if (tmin <= t2 && t2 <= tmax) t = t2;
    else return false;
    hitPos = ray.o + t * ray.d;
    hitNormal = normalize(hitPos - sph.center);
    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere(in Sphere_t sph, in Ray_t ray, in float tmin, in float tmax) {
    float a = dot(ray.d, ray.d);
    float b = 2.0 * dot(ray.o - sph.center, ray.d);
    float c = dot(ray.o - sph.center, ray.o - sph.center) - sph.radius * sph.radius;
    float det_square = b * b - 4.0 * a * c;
    if (det_square < 0.0) return false;
    float det = sqrt(det_square);
    float t1 = (-b - det) / (2.0 * a);
    float t2 = (-b + det) / (2.0 * a);
    if (tmin <= t1 && t1 <= tmax) return true;
    else if (tmin <= t2 && t2 <= tmax) return true;
    else return false;

}

bool IntersectCube(in Cube_t cube, in Ray_t ray, in float tmin, in float tmax, out float t, out vec3 hitPos, out vec3 hitNormal) {
    vec3 minBound = cube.center - (cube.size / 2.0); // Compute the minimum boundaries of the cube.
    vec3 maxBound = cube.center + (cube.size / 2.0); // Compute the maximum boundaries of the cube.
    // Compute the intersection of ray and planes of cube
    vec3 t1 = (minBound - ray.o) / ray.d;
    vec3 t2 = (maxBound - ray.o) / ray.d;
    // Compute the smallest and largest t that intersects the cube
    vec3 t3 = min(t1, t2);
    vec3 t4 = max(t1, t2);
    // Compute entry and exit points.
    float tstart = max(max(t3.x, t3.y), max(t3.y, t3.z));
    float tend = min(min(t4.x, t4.y), min(t4.y, t4.z));
    if(tend < tstart || tstart > tmax || tend < tmin)
        return false;
        
    t = tstart < tmin ? tend : tstart;
    // Compute the hit position
    hitPos = ray.o + ray.d * t;
    // Compute the length between the cube center and the hit position.
    vec3 d = abs(hitPos - cube.center) / (cube.size / 2.0);
    float maxD = max(max(d.x, d.y), d.z); // Maximum dimension

    if(maxD == d.x)
        hitNormal = vec3(sign(hitPos.x - cube.center.x), 0.0, 0.0);
    else if(maxD == d.y)
        hitNormal = vec3(0.0, sign(hitPos.y - cube.center.y), 0.0);
    else
        hitNormal = vec3(0.0, 0.0, sign(hitPos.z - cube.center.z));
    hitNormal = normalize(hitNormal);

    return true;
}


bool IntersectCube( in Cube_t cube, in Ray_t ray, in float tmin, in float tmax) {
    vec3 minBound = cube.center - (cube.size / 2.0);
    vec3 maxBound = cube.center + (cube.size / 2.0);

    vec3 t1 = (minBound - ray.o) / ray.d;
    vec3 t2 = (maxBound - ray.o) / ray.d;
    
    vec3 t3 = min(t1, t2);
    vec3 t4 = max(t1, t2);
    
    float tstart = max(max(t3.x, t3.y), max(t3.y, t3.z));
    float tend = min(min(t4.x, t4.y), min(t4.y, t4.z));
    
    if(tend < tstart || tstart > tmax || tend < tmin)
        return false;

    return true;
}
Hit IntersectScene(Ray_t ray) { 
    bool hasHitSomething = false;
    float nearest_t = ray.t;   // The ray parameter t at the nearest hit point.
    vec3 nearest_hitPos;              // 3D position of the nearest hit point.
    vec3 nearest_hitNormal;           // Normal vector at the nearest hit point.
    int nearest_hitMatID;             // MaterialID of the object at the nearest hit point.

    float temp_t;
    vec3 temp_hitPos;
    vec3 temp_hitNormal;
    bool temp_hasHit;
    
    for (int i = 0; i < NUM_PLANES; i++) {
        if (IntersectPlane(Plane[i], ray, DEFAULT_TMIN, nearest_t, temp_t, temp_hitPos, temp_hitNormal) == true) {
            temp_hasHit = true;
            if (temp_t < nearest_t) {
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Plane[i].materialID;
                hasHitSomething = temp_hasHit;
            }
        }
    }
    for (int i = 0; i < NUM_SPHERES; i++) {
        if (IntersectSphere(Sphere[i], ray, DEFAULT_TMIN, nearest_t, temp_t, temp_hitPos, temp_hitNormal) == true) {
            temp_hasHit = true;
            if (temp_t < nearest_t) {
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Sphere[i].materialID;
                hasHitSomething = temp_hasHit;
            }
        }
    }
    for (int i = 0; i < NUM_CUBES; i++) {
        if(IntersectCube(Cube[i], ray, DEFAULT_TMIN, nearest_t,
                      temp_t, temp_hitPos, temp_hitNormal)) {
                         hasHitSomething = true;
                         nearest_t = temp_t;
                         nearest_hitPos = temp_hitPos;
                         nearest_hitNormal = temp_hitNormal;
                         nearest_hitMatID = Cube[i].materialID;
                      }
    }
    if(nearest_hitMatID==10 && IntersectNUS(nearest_hitPos.zy))nearest_hitMatID = 11;
    if(nearest_hitMatID==12) getFloorMaterial(nearest_hitPos.xz);
    return Hit(hasHitSomething, nearest_hitPos, nearest_hitNormal, nearest_hitMatID);
}
//---Random------------------------------------------------------
float random (vec2 st) {
    return fract(sin(dot(st.xy,vec2(12.9898,78.233)))*43758.5453123);
}

float randomInt(float seed) {
    vec2 temp = vec2(seed, 0.0);
    float rand = random(temp);
    rand = rand * 6.0; 
    rand = ceil(rand);
    return rand;
}

vec3 mapHolePos(in float num) {
    if (num == 1.0) {
        return vec3(18.0, 0.5, -9.0);
    }else if(num == 2.0) {
        return vec3(18.0, 0.5, 9.0);
    }else if(num == 3.0) {
        return vec3(0.0, 0.5, -9.0);
    }else if(num == 4.0) {
        return vec3(0.0, 0.5, 9.0);
    }else if(num == 5.0) {
        return vec3(-18.0, 0.5, -9.0);
    }else if(num == 6.0) {
        return vec3(-18.0, 0.5, 9.0);
    }
}

vec3 mapWhiteBall(in float num) {
    if (num == 1.0) {
        return vec3(-5.0, 0.5, 4.0);
    }else if(num == 2.0) {
        return vec3(-5.0, 0.5, -4.0);
    }else if(num == 3.0) {
        return vec3(15.0, 0.5, 6.0);
    }else if(num == 4.0) {
        return vec3(15.0, 0.5, -6.0);
    }else if(num == 5.0) {
        return vec3(15.0, 0.5, 2.0);
    }else if(num == 6.0) {
        return vec3(15.0, 0.5, -2.0);
    }
}

//----------------------------------------------------------------------

//---Movement-----------------------------------------------------------
void CalcMove(in int sph, in float rand) {
    vec3 holePos = mapHolePos(rand);
    Sphere[21].center = mapWhiteBall(rand);
    vec3 centerPos = Sphere[sph].center;
    vec3 TargetMove = holePos - centerPos;
    float TargetMovex = TargetMove.x;
    float TargetMovez = TargetMove.z;
    float TargetMoveLength = sqrt(TargetMovex * TargetMovex + TargetMovez * TargetMovez);
    float S = TargetMovez / TargetMoveLength;
    float C = TargetMovex / TargetMoveLength;
    vec3 arrivePos = vec3(centerPos.x - 1.0 * C, 0.5, centerPos.z - 1.0 * S);
    vec3 HitMove = arrivePos - Sphere[21].center;
    float HitMovex = abs(HitMove.x);
    float HitMovez = abs(HitMove.z);
    float HitMoveLength = sqrt(HitMovex * HitMovex + HitMovez * HitMovez);
    float s = HitMovez / HitMoveLength;
    float c = HitMovex / HitMoveLength;
    float T = 3.0;
    float t = mod(iTime, T);
    float Speed = (HitMoveLength + 2.0 * TargetMoveLength) / T;
    float DownSpeed = 1.0 / 0.15;
    float t1 = HitMoveLength / Speed;
    float xDir = 1.0;
    if (Sphere[sph].center.x < Sphere[21].center.x) {
        xDir = -1.0;
    }
    float zDir = 1.0;
    if (Sphere[sph].center.z < Sphere[21].center.z) {
        zDir = -1.0;
    }
    if (t <= t1) {
        
        Sphere[21].center.x += xDir * Speed * t * c;
        Sphere[21].center.z += zDir * Speed * t * s;
    }
    else {
        Sphere[21].center = arrivePos;
        float delta = t - t1;
        Sphere[sph].center.x += Speed * 0.5 * delta * C;
        Sphere[sph].center.z += Speed * 0.5 * delta * S;
        if (t >= 2.8) {
            float downDelta = t - 2.85;
            Sphere[sph].center.y -= DownSpeed * downDelta;
        }
        vec3 dir = HitMove - 0.5 * TargetMove;
        float dirLength = sqrt(dir.x * dir.x + dir.z * dir.z);
        float ss = dir.z / dirLength;
        float cc = dir.x / dirLength;
        Sphere[21].center.x += Speed * 0.13 * delta * cc;
        Sphere[21].center.z += Speed * 0.13 * delta * ss;
    }
}

int Movement() {
    float T = 3.0;
    int turn = int(mod(iTime / T, 4.0));
    float t = floor((iTime + 3.0) / T);
    float rand = randomInt(t);
    
    CalcMove(17 + turn, rand);
    return 17 + turn;
}
//----------------------------------------------------------------------

vec3 Whitted_Raytracing() {
    // Scale pixel 2D position such that its y coordinate is in [-1.0, 1.0].
    vec3 result = vec3(0);
    vec3 pixel_pos = vec3((2.0 * gl_FragCoord.xy - iResolution.xy) / iResolution.y, -1.0 / tan(FOVY / 2.0));
    Ray_t pRay;
    pRay.o = camera.pos;
    pRay.d = normalize(pixel_pos.x * camera.axis[0]  +  pixel_pos.y * camera.axis[1]  +  pixel_pos.z * camera.axis[2]);
    pRay.t = DEFAULT_TMAX;
    
    int bounce = 1;
    vec3 throughput = vec3(1);

    while (true) {
        Hit hit = IntersectScene(pRay);
        if(!hit.hit)break;
        Material_t m = Material[hit.MaterialID];
        vec3 n = hit.normal;
        vec3 p = hit.position + hit.normal * EPSILON;

        // Add hacky emission on first hit to draw lights
        if (bounce == 1)
            result += throughput * m.emission;

        // Sample lights
        for (int i = 0; i < NUM_LIGHTS; i++) {
            vec3 lightPos = Light[i].position;

            // Generate shadow ray
            Ray_t sray;
            sray.o = p;
            sray.t = length(lightPos - p);
            sray.d = (lightPos - p) / sray.t;

            // Test visibility
            Hit sHit = IntersectScene(sray);
            if (!sHit.hit) {
                // Add light contribution when visible
                float rSquare = sray.t * sray.t;
                vec3 E = Light[i].E;
                result += throughput * evalBRDF(hit.normal, -pRay.d, sray.d, m) * E / rSquare;
                
            }
        }

        // Russian roulette for termination
        if (bounce >= BOUNCES)
            break;
        
        throughput *= evalBRDF(hit.normal, -pRay.d, reflect(pRay.d, n), m);
        pRay.d = reflect(pRay.d, n);
        pRay.o = p;
        bounce++;
    }
    return result;
}


void SetCamera(in int num) {
    
    // Position the camera.
    vec3 OriginPos = vec3(-4.0, 0.5, 6.0);
    camera.pos = vec3( OriginPos.x + (Sphere[21].center.x - OriginPos.x) * 0.6 - 23.0, 12.5, OriginPos.z + (Sphere[21].center.z - OriginPos.z) * 0.6 + 15.0 );
    camera.lookat = Sphere[num].center;
    vec3 cam_up_vec = vec3(0.0, 1.0, 0.0);

    // Set up camera coordinate frame in world space.
    camera.axis[2] = normalize(camera.pos - camera.lookat);
    camera.axis[0] = normalize(cross(cam_up_vec, camera.axis[2]));
    camera.axis[1] = normalize(cross(camera.axis[2], camera.axis[0]));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord)
{
    InitScene();
    int ballNo = Movement();
    SetCamera(ballNo);
    vec3 result = Whitted_Raytracing();
    fragColor = vec4(ltos3(result.x, result.y, result.z), 1.0); // Gamma correct
}
