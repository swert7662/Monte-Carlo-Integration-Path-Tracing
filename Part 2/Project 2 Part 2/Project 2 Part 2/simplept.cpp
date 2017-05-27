#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#include <iostream>
#define PI 3.1415926535897932384626433832795

/*
* Thread-safe random number generator
*/

struct RNG {
	RNG() : distrb(0.0, 1.0), engines() {}

	void init(int nworkers) {
		std::random_device rd;
		engines.resize(nworkers);
		for (int i = 0; i < nworkers; ++i)
			engines[i].seed(rd());
	}

	double operator()() {
		int id = omp_get_thread_num();
		return distrb(engines[id]);
	}

	std::uniform_real_distribution<double> distrb;
	std::vector<std::mt19937> engines;
} rng;


/*
* Basic data types
*/

struct Vec {
	double x, y, z;

	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }

	Vec operator+ (const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator- (const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator* (double b) const { return Vec(x*b, y*b, z*b); }

	Vec mult(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec& normalize() { return *this = *this * (1.0 / std::sqrt(x*x + y*y + z*z)); }
	double dot(const Vec &b) const { return x*b.x + y*b.y + z*b.z; }
	Vec cross(const Vec&b) const { return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
};

struct Ray {
	Vec o, d;
	Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

struct BRDF {
	virtual Vec eval(const Vec &n, const Vec &o, const Vec &i) const = 0;
	virtual void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const = 0;
	virtual bool isSpecular() const = 0;
};


/*
* Utility functions
*/

inline double clamp(double x) {
	return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x) {
	return static_cast<int>(std::pow(clamp(x), 1.0 / 2.2) * 255 + .5);
}


/*
* Shapes
*/

struct Sphere {
	Vec p, e;           // position, emitted radiance
	double rad;         // radius
	const BRDF &brdf;   // BRDF

	Sphere(double rad_, Vec p_, Vec e_, const BRDF &brdf_) :
		rad(rad_), p(p_), e(e_), brdf(brdf_) {}

	double intersect(const Ray &r) const { // returns distance, 0 if nohit
		Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, eps = 1e-4, b = op.dot(r.d), det = b*b - op.dot(op) + rad*rad;
		if (det<0) return 0; else det = sqrt(det);
		return (t = b - det)>eps ? t : ((t = b + det)>eps ? t : 0);
	}
};


/*
* Sampling functions
*/

inline void createLocalCoord(const Vec &n, Vec &u, Vec &v, Vec &w) {
	w = n;
	u = ((std::abs(w.x)>.1 ? Vec(0, 1) : Vec(1)).cross(w)).normalize();
	v = w.cross(u);
}

Vec uniformRandPSA(const Vec &n) {
	double z = sqrt(rng());
	double r = sqrt(1.0 - z * z);
	double phi = 2.0 * PI * rng();
	double x = r * cos(phi);
	double y = r * sin(phi);
	Vec u, v, w;
	createLocalCoord(n, u, v, w);
	return u * x + v * y + w * z;
}
/*
* BRDFs
*/

// Ideal diffuse BRDF
struct DiffuseBRDF : public BRDF {
	DiffuseBRDF(Vec kd_) : kd(kd_) {}

	Vec eval(const Vec &n, const Vec &o, const Vec &i) const {
		return kd * (1.0 / PI);
	}

	void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {
		i = uniformRandPSA(n);
		pdf = i.dot(n) * (1 / PI);
	}
	bool isSpecular() const{
		return false;
	}
	Vec kd;
};

// Ideal Specular BRDF
struct SpecularBRDF : public BRDF {
	SpecularBRDF(Vec ks_) : ks(ks_) {}

	Vec mirroredDir(const Vec &n, const Vec &o) const{
		return ((n * 2.0) * (n.dot(o))) - o;
	}

	Vec eval(const Vec &n, const Vec &o, const Vec &i) const {
		Vec m = mirroredDir(n, i);
		if (abs(o.x - m.x) < 1e-5 && abs(o.y - m.y) < 1e-5 && abs(o.z - m.z) < 1e-5) //o == m
			return (ks * (1 / n.dot(i)));
		else
			return 0.0;
	}

	void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {
		i = mirroredDir(n, o).normalize();
		pdf = 1.0;
	}
	bool isSpecular() const{
		return true;
	}
	Vec ks;
};


/*
* Scene configuration
*/

// Pre-defined BRDFs
const DiffuseBRDF leftWall(Vec(.75, .25, .25)),
rightWall(Vec(.25, .25, .75)),
otherWall(Vec(.75, .75, .75)),
blackSurf(Vec(0.0, 0.0, 0.0)),
brightSurf(Vec(0.9, 0.9, 0.9));
const SpecularBRDF shiny(Vec(0.999, 0.999, 0.999));

// Scene: list of spheres
const Sphere spheres[] = {
	Sphere(1e5,  Vec(1e5 + 1,40.8,81.6),   Vec(),         leftWall),   // Left
	Sphere(1e5,  Vec(-1e5 + 99,40.8,81.6), Vec(),         rightWall),  // Right
	Sphere(1e5,  Vec(50,40.8, 1e5),      Vec(),         otherWall),  // Back
	Sphere(1e5,  Vec(50, 1e5, 81.6),     Vec(),         otherWall),  // Bottom
	Sphere(1e5,  Vec(50,-1e5 + 81.6,81.6), Vec(),         otherWall),  // Top
	Sphere(16.5, Vec(27,16.5,47),        Vec(),         brightSurf), // Ball 1
	Sphere(16.5, Vec(73,16.5,78),        Vec(),         shiny), // Ball 2
	Sphere(5.0,  Vec(50,70.0,81.6),      Vec(50,50,50), blackSurf)   // Light
};

// Camera position & direction
const Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());


/*
* Global functions
*/

bool intersect(const Ray &r, double &t, int &id) {
	double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
	for (int i = int(n); i--;) if ((d = spheres[i].intersect(r)) && d<t) { t = d; id = i; }
	return t<inf;
}


void luminaireSample(const Vec &X, Vec &Y, Vec &ny, double &pdf, Vec &light)
{
	const Sphere l = spheres[7];
	double z = 2 * rng() - 1.0;
	double xi2 = rng();
	double x = sqrt(1.0 - z * z) * sin(2 * PI * xi2);
	double y = sqrt(1.0 - z * z) * cos(2 * PI * xi2);

	Y = l.p + (Vec(x, y, z)*l.rad);
	pdf = 1 / (4 * PI *l.rad * l.rad);
	ny = Vec(x, y, z);
	light = l.e;

}

double visiblity(const Vec &x, const Vec &y)
{
	Ray r(x, (y - x).normalize());
	double t;
	int id = 0;
	intersect(r, t, id);
	Vec i = r.o + r.d*t;
	if (abs(i.x - y.x) < 1e-5 && abs(i.y - y.y) < 1e-5 && abs(i.z - y.z) < 1e-5) //i == y
		return 1.0;
	return 0.0;
};
/*
* KEY FUNCTION: radiance estimator
*/
Vec receivedRadiance(const Ray &r, int depth, bool flag) {
	double t;                                   // Distance to intersection
	int id = 0;                                 // id of intersected sphere

	if (!intersect(r, t, id)) return Vec();   // if miss, return black
	const Sphere &obj = spheres[id];            // the hit object
	const BRDF &brdf = obj.brdf;                // Surface BRDF at x
	Vec rad = obj.e;							// Emitted radiance

	Vec x = r.o + r.d*t;                        // The intersection point
	Vec o = (Vec() - r.d).normalize();          // The outgoing direction (= -r.d)
	
	Vec n = (x - obj.p).normalize();            // The normal direction
	if (n.dot(o) < 0) n = n*-1.0;

	const int rrd = 5;
	const double sp = 0.9;
	double p = 1.0;

	Vec directRadiance;
	//Direct Radiance Section
	if (!brdf.isSpecular())      
	{
		Vec y1, ny, light;
		double pdf1;
		luminaireSample(x, y1, ny, pdf1, light);
		Vec o1 = (y1 - x).normalize();
		double r2 = (x - y1).dot(x - y1);
		directRadiance = light.mult(brdf.eval(n, o1, o)) * visiblity(x, y1) * n.dot(o1) * (ny.dot(o1*-1) / (r2 * pdf1));
	}

	//Indirect Radiance Section
	if (depth > rrd)
		p = sp;
	Vec indirectRadiance;
	if (rng() < p)
	{
		Vec o2;
		double pdf2;
		brdf.sample(n, o, o2, pdf2);
		Ray y(x, o2.normalize());
		indirectRadiance = receivedRadiance(y, depth + 1, brdf.isSpecular()).mult(brdf.eval(n, o, o2))*(n.dot(o2)/(pdf2*p));
	}
	if (flag)
		return rad + directRadiance + indirectRadiance;
	return directRadiance + indirectRadiance;
}

/*
* Main function (do not modify)
*/

int main(int argc, char *argv[]) {
	int nworkers = omp_get_num_procs();
	omp_set_num_threads(nworkers);
	rng.init(nworkers);

	int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
	Vec cx = Vec(w*.5135 / h), cy = (cx.cross(cam.d)).normalize()*.5135;
	std::vector<Vec> c(w*h);

#pragma omp parallel for schedule(dynamic, 1)
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			const int i = (h - y - 1)*w + x;

			for (int sy = 0; sy < 2; ++sy) {
				for (int sx = 0; sx < 2; ++sx) {
					Vec r;
					for (int s = 0; s<samps; s++) {
						double r1 = 2 * rng(), dx = r1<1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double r2 = 2 * rng(), dy = r2<1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
						Vec d = cx*(((sx + .5 + dx) / 2 + x) / w - .5) +
							cy*(((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
						r = r + receivedRadiance(Ray(cam.o, d.normalize()), 1, true)*(1. / samps);
					}
					c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
				}
			}
		}
#pragma omp critical
		fprintf(stderr, "\rRendering (%d spp) %6.2f%%", samps * 4, 100.*y / (h - 1));
	}
	fprintf(stderr, "\n");

	// Write resulting image to a PPM file
	FILE *f = fopen("image.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i<w*h; i++)
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	fclose(f);

	return 0;
}