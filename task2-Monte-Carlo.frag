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
const int SAMPLES_PER_PIXEL = 200;
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

float seed = 1.0; //seed initialized in main
float rnd() { return fract(sin(seed++)*43758.5453123); }

const vec3 BACKGROUND_COLOR = vec3( 0.1, 0.2, 0.6 );

// Vertical field-of-view angle of camera. In radians.
const float FOVY = 50.0 * PI / 180.0;

// Use this for avoiding the "epsilon problem" or the shadow acne problem.
const float DEFAULT_TMIN = 10.0e-4;

// Use this for tmax for non-shadow ray intersection test.
const float DEFAULT_TMAX = 10.0e6;


// Constants for the scene objects.
const int NUM_LIGHTS = 3;
const int NUM_MATERIALS = 11;
const int NUM_SPHERES = 22;
const int NUM_CUBES = 7;
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
    mat4 worldMat;
    vec2 size;
    vec3 E;
    vec3 normal;
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
    Plane[0].D = 3.0;
    Plane[0].materialID = 10;
    
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
    Plane[3].D = 20.0;
    Plane[3].materialID = 10;
    
    // Back Plane
    Plane[4].A = 0.0;
    Plane[4].B = 0.0;
    Plane[4].C = 1.0;
    Plane[4].D = 20.0;
    Plane[4].materialID = 10;
}

void InitBall() {
    // Black Ball
    Sphere[0].center = vec3(-15.0, 0.5, 0.0);
    Sphere[0].radius = 0.5;
    Sphere[0].materialID = 0;
    

    // Red Ball * 15
    for (int i = 1; i <= 15; i++) {
        Sphere[i].radius = 0.5;
        Sphere[i].materialID = 1;
    }
    Sphere[1].center = vec3(-12.0, 0.5, 2.0);
    
    Sphere[2].center = vec3(-12.0, 0.5, 1.0);
    
    Sphere[3].center = vec3(-12.0, 0.5, 0.0);
    
    Sphere[4].center = vec3(-12.0, 0.5, -1.0);
    
    Sphere[5].center = vec3(-12.0, 0.5, -2.0);
    
    Sphere[6].center = vec3(-11.13, 0.5, 1.5);
    
    Sphere[7].center = vec3(-11.14, 0.5, 0.5);
    
    Sphere[8].center = vec3(-11.14, 0.5, -0.5);
    
    Sphere[9].center = vec3(-11.14, 0.5, -1.5);
    
    Sphere[10].center = vec3(-10.26, 0.5, 1.0);
    
    Sphere[11].center = vec3(-10.28, 0.5, 0.0);
    
    Sphere[12].center = vec3(-10.28, 0.5, -1.0);
    
    Sphere[13].center = vec3(-9.39, 0.5, 0.5);
    
    Sphere[14].center = vec3(-9.39, 0.5, -0.5);
    
    Sphere[15].center = vec3(-8.52, 0.5, 0.0);
    

    // Purple Ball
    Sphere[16].center = vec3(-7.52, 0.5, 0.0);
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
    Sphere[21].center = vec3(3.0, 0.5, 4.5);
    Sphere[21].radius = 0.5;
    Sphere[21].materialID = 7;

}

void InitCube() {
    // Back Left Baffle
    Cube[0].center = vec3(-8.75, 0.5, -9.5);
    Cube[0].size = vec3(16.3, 1.0, 1.0);
    Cube[0].materialID = 9;
    
    // Back Right Baffle
    Cube[1].center = vec3(8.75, 0.5, -9.5);
    Cube[1].size = vec3(16.3, 1.0, 1.0);
    Cube[1].materialID = 9;
    
    // Front Left Baffle
    Cube[2].center = vec3(-8.75, 0.5, 9.5);
    Cube[2].size = vec3(16.3, 1.0, 1.0);
    Cube[2].materialID = 9;
    
    // Front Right Baffle
    Cube[3].center = vec3(8.75, 0.5, 9.5);
    Cube[3].size = vec3(16.3, 1.0, 1.0);
    Cube[3].materialID = 9;
    
    // Left Baffle
    Cube[4].center = vec3(-18.5, 0.5, 0.0);
    Cube[4].size = vec3(1.0, 1.0, 15.8);
    Cube[4].materialID = 9;
    
    // Right Baffle
    Cube[5].center = vec3(18.5, 0.5, 0.0);
    Cube[5].size = vec3(1.0, 1.0, 15.8);
    Cube[5].materialID = 9;
    
    // Table Plane
    Cube[6].center = vec3(0.0, -0.05, 0.0);
    Cube[6].size = vec3(38.0, 0.1, 20.0);
    Cube[6].materialID = 8;
}

void InitLight() {
    // Light 0
    Light[0].worldMat = mat4(1, 0, 0, 0,
                             0, 1, 0, 0,
                             0, 0, 1, 0,
                             0, 5, 0, 1);
    Light[0].size = vec2(20);
    Light[0].E = vec3(1, 1, 1) * vec3(3);
    Light[0].normal = vec3(0, -1 ,0);
    
    // Light 1
    Light[1].worldMat = mat4(1, 0, 0, 0,
                             0, 1, 0, 0,
                             0, 0, 1, 0,
                             0, 100, 0, 1);
    Light[1].size = vec2(20);
    Light[1].E = vec3(1, 1, 1) * vec3(5);
    Light[1].normal = normalize(vec3(0, -1, -1));
    // Light 2
    Light[2].worldMat = mat4(1, 0, 0, 0,
                             0, 1, 0, 0,
                             0, 0, 1, 0,
                             0, 20, -20, 1);
    Light[2].size = vec2(2);
    Light[2].E = vec3(1, 1, 1) * vec3(100);
    Light[2].normal = normalize(vec3(0, -1, 1));


}
void InitMaterial() {
    
    // Black Plastic Material
    Material[0].albedo = vec3(0, 0, 0);
    Material[0].roughness = 0.8;
    Material[0].metalness = 0.2;
    Material[0].emission = vec3(0.0);
    // Red Plastic Material
    Material[1].albedo = vec3(0.89, 0.09, 0.05);
    Material[1].roughness = 0.8;
    Material[1].metalness = 0.2;
    Material[1].emission = vec3(0.0);
    // Purple Plastic Material
    Material[2].albedo = vec3(0.86, 0.44, 0.84);
    Material[2].roughness = 0.8;
    Material[2].metalness = 0.2;
    Material[2].emission = vec3(0.0);
    // Blue Plastic Material
    Material[3].albedo = vec3(0.0, 0.0, 1.0);
    Material[3].roughness = 0.8;
    Material[3].metalness = 0.2;
    Material[3].emission = vec3(0.0);
    // Green Plastic Material
    Material[4].albedo = vec3(0.13, 0.55, 0.13);
    Material[4].roughness = 0.8;
    Material[4].metalness = 0.2;
    Material[4].emission = vec3(0.0);
    // Brown Plastic Material
    Material[5].albedo = vec3(0.78, 0.38, 0.08);
    Material[5].roughness = 0.8;
    Material[5].metalness = 0.2;
    Material[5].emission = vec3(0.0);
    // Yellow Plastic Material
    Material[6].albedo = vec3(1.0, 1.0, 0.0);
    Material[6].roughness = 0.8;
    Material[6].metalness = 0.2;
    Material[6].emission = vec3(0.0);
    // White Plastic Material
    Material[7].albedo = vec3(1.0, 1.0, 1.0);
    Material[7].roughness = 0.2;
    Material[7].metalness = 0.8;
    Material[7].emission = vec3(0.0);
    // Table Material
    Material[8].albedo = vec3(0.42, 0.56, 0.14);
    Material[8].roughness = 0.8;
    Material[8].metalness = 0.5;
    Material[8].emission = vec3(0.0);
    // Baffle Material
    Material[9].albedo = vec3(0.18, 0.17, 0.23);
    Material[9].roughness = 1.0;
    Material[9].metalness = 0.2;
    Material[9].emission = vec3(0.0);
    // Room Material
    Material[10].albedo = vec3(1, 1, 1);
    Material[10].roughness = 0.3;
    Material[10].metalness = 1.0;
    Material[10].emission = vec3(0.4);
}
/////////////////////////////////////////////////////////////////////////////
// Initializes the scene.
/////////////////////////////////////////////////////////////////////////////

void InitScene() {
    seed = iTime + gl_FragCoord.y * gl_FragCoord.x / iResolution.x + gl_FragCoord.y / iResolution.y;
    InitPlane();
    InitBall();
    InitCube();
    InitLight();
    InitMaterial();
}

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

// --SAMPLING------------------------------------------------------------------
vec4 sampleLight(int i) {
    Light_t light = Light[i];
    float pdf = 1.0 / (4.0 * light.size.x * light.size.y);
    mat4 S = mat4(light.size.x, 0, 0, 0,
                  0, light.size.y, 0, 0,
                  0,            0, 1, 0,
                  0,            0, 0, 1);
    mat4 M = light.worldMat * S;
    return vec4((M * vec4(vec2(rnd(), rnd()) * 2.0 - 1.0, 0, 1)).xyz, pdf);
}

// From http://www.rorydriscoll.com/2009/01/07/better-sampling/
vec3 cosineSampleHemisphere() {
    vec2 u = vec2(rnd(), rnd());
    float r = sqrt(u.x);
    float theta = 2.0 * PI * u.y;
    return vec3(r * cos(theta), r * sin(theta), sqrt(saturate(1.0 - u.x)));
}

float cosineHemispherePDF(float NoL) {
    return NoL / PI;
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
    return Hit(hasHitSomething, nearest_hitPos, nearest_hitNormal, nearest_hitMatID);
}

//----------------------------------------------------------------------

//---Movement-----------------------------------------------------------
void CalcMove(in int sph) {
    vec3 holePos = vec3(18.0, 0.5, -9.0);
    vec3 centerPos = Sphere[sph].center;
    vec3 TargetMove = holePos - centerPos;
    float TargetMovex = abs(TargetMove.x);
    float TargetMovez = abs(TargetMove.z);
    float TargetMoveLength = sqrt(TargetMovex * TargetMovex + TargetMovez * TargetMovez);
    float S = TargetMovez / TargetMoveLength;
    float C = TargetMovex / TargetMoveLength;
    vec3 arrivePos = vec3(centerPos.x - 1.0 * C, 0.5, centerPos.z + 1.0 * S);
    vec3 HitMove = arrivePos - Sphere[21].center;
    float HitMovex = abs(HitMove.x);
    float HitMovez = abs(HitMove.z);
    float HitMoveLength = sqrt(HitMovex * HitMovex + HitMovez * HitMovez);
    float s = HitMovez / HitMoveLength;
    float c = HitMovex / HitMoveLength;
    float T = 3.0;
    float t = mod(iTime, T);
    float Speed = (HitMoveLength + 2.0 * TargetMoveLength) / T;
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
        Sphere[sph].center.z -= Speed * 0.5 * delta * S;
        vec3 dir = HitMove - 0.5 * TargetMove;
        float dirLength = sqrt(dir.x * dir.x + dir.z * dir.z);
        float ss = dir.z / dirLength;
        float cc = dir.x / dirLength;
        Sphere[21].center.x += Speed * 0.15 * delta * cc;
        Sphere[21].center.z += Speed * 0.15 * delta * ss;
    }
}

void Movement() {
    float T = 3.0;
    int turn = int(mod(iTime / T, 4.0));
    CalcMove(17 + turn);
}
//----------------------------------------------------------------------

bool rayOutScene(vec3 pos){
    //return pos.z < -5.0; // Change it!
    return false;
}

vec3 Monte_Carlo_Raytracing() {
    // Scale pixel 2D position such that its y coordinate is in [-1.0, 1.0].
    vec3 result = vec3(0);
    for (int j = 0; j < SAMPLES_PER_PIXEL; j++) {
        // Generate ray
        vec2 sample_px = gl_FragCoord.xy + vec2(rnd(), rnd());
        vec3 pixel_pos = vec3((2.0 * sample_px.xy - iResolution.xy) / iResolution.y, -1.0 / tan(FOVY / 2.0));
        Ray_t pRay;
        pRay.o = camera.pos;
        pRay.d = normalize(pixel_pos.x * camera.axis[0]  +  pixel_pos.y * camera.axis[1]  +  pixel_pos.z * camera.axis[2]);
        pRay.t = DEFAULT_TMAX;
        
        int bounce = 1;
        vec3 throughput = vec3(1);

        while (true) {
            Hit hit = IntersectScene(pRay);
            if(!hit.hit || rayOutScene(hit.position))break;
            Material_t m = Material[hit.MaterialID];
            vec3 n = hit.normal;
            vec3 p = hit.position + hit.normal * EPSILON;

            // Add hacky emission on first hit to draw lights
            if (bounce == 1)
                result += throughput * m.emission;

            // Sample lights
            for (int i = 0; i < NUM_LIGHTS; i++) {
                // Generate point on light surface
                vec4 ls = sampleLight(i);
                vec3 lightPos = ls.xyz;
                float pdf = ls.w;

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
                    if (dot(Light[i].normal, -sray.d) > 0.0) {
                    	vec3 E = Light[i].E;
                    	result += throughput * evalBRDF(hit.normal, -pRay.d, sray.d, m) * E / (rSquare * pdf);
                    }
                }
            }

            // Russian roulette for termination
            if (bounce >= BOUNCES && rnd() < P_TERMINATE)
                break;

            // Get random direction for reflection ray
            vec3 randomDir = cosineSampleHemisphere();
            // Rotate by normal frame
            randomDir = normalize(formBasis(n) * randomDir);
            float pdf = cosineHemispherePDF(dot(n, randomDir));
            // TODO: Multiple importance sampling on diffuse and specular?
            throughput *= evalBRDF(hit.normal, -pRay.d, randomDir, m) / pdf;
            pRay.d = randomDir;
            pRay.o = p;
            bounce++;
        }
    }
    return result / float(SAMPLES_PER_PIXEL);
}


void SetCamera() {
    // Position the camera.
    camera.pos = vec3( Sphere[21].center.x - 20.0, 10.0, Sphere[21].center.z + 15.0 );
    camera.lookat = Sphere[21].center;
    vec3 cam_up_vec = vec3(0.0, 1.0, 0.0);

    // Set up camera coordinate frame in world space.
    camera.axis[2] = normalize(camera.pos - camera.lookat);
    camera.axis[0] = normalize(cross(cam_up_vec, camera.axis[2]));
    camera.axis[1] = normalize(cross(camera.axis[2], camera.axis[0]));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord)
{
    InitScene();
    Movement();
    SetCamera();
    vec3 result = Monte_Carlo_Raytracing();
    fragColor = vec4(ltos3(result.x, result.y, result.z), 1.0); // Gamma correct
}
