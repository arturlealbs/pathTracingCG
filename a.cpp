#include <climits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

#include "atom.h"
#include "scene.h"
#define ccout std::cout <<
#define ccin std::cin >>
using namespace std;
#define N_THREADS 3
#define MAX_RECURSION_DEPTH 5
#define PI 3.14159265358979323846
// Buffer buf;
// Recebe 800x600 como defalt mas � alterada dependendo do arquivo de entrada
// static GLfloat window_width = 800.0;
// static GLfloat window_height = 600.0;

// Loading camera paramenters
point olho;

// Loading view paramenters
Screen janela;

// Loading scene paramenters
Scene cena;

// Loading illumination paramenters
Light luz;

// color** texture;
// int textureX = 300;
// int textureY = 300;

// Loading objects of the scene
std::vector<Object> objetos;

class Intersection {
 public:
  Object objeto;
  point p;
  vec normal;
  bool hit;

 private:
  Face f;
};

//Intersection bInt1[300][300];
point testPoints[500];
int countPoints = 0;

// Color texture[50][50];

///-------------------------OTHER PEOPLE'S
/// STUFF-----------------------------------------///
// http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
// http://www.cs.virginia.edu/~gfx/Courses/2003/ImageSynthesis/papers/Acceleration/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf


class Triangle {
public:
    vec v0, v1, v2;
    vec normal;

    Triangle(vec v0_, vec v1_, vec v2_) : v0(v0_), v1(v1_), v2(v2_) {
        // Calculate normal
        vec e1 = v1 - v0;
        vec e2 = v2 - v0;
        normal = e1.cross(e2).normalized();
    }

    double intersect(const Ray &ray) const {
        double t;
        const double EPSILON = 0.0000001;
        vec edge1, edge2, h, s, q;
        double a, f, u, v;
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        h = ray.direcao.cross(edge2);
        a = edge1.dot(h);
        if (a > -EPSILON && a < EPSILON)
            return INT_MAX; // Ray is parallel to this triangle.
        f = 1.0 / a;
        s = ray.posicao - v0;
        u = f * s.dot(h);
        if (u < 0.0 || u > 1.0)
            return INT_MAX;
        q = s.cross(edge1);
        v = f * ray.direcao.dot(q);
        if (v < 0.0 || u + v > 1.0)
            return INT_MAX;
        // Compute t to find out where the intersection point is on the line.
        t = f * edge2.dot(q);
        if (t > EPSILON) // ray intersection
            return t;
        else // There is a line intersection but not a ray intersection.
            return t;
    }
};

///------------------------- END OF OTHER PEOPLE'S
/// STUFF-----------------------------------------///

bool retorno = false;
double dist = INT_MAX;
double distTemp = -1;

vec calcularRefracao(double n1, double n2, vec i, vec n) {
  float cosI = -i.dot(n);
  float sen2t = std::pow(n1 / n2, 2) * (1 - std::pow(cosI, 2));

  vec t = (((n1 / n2) * i) + (((n1 / n2) * cosI - std::sqrt(1 - sen2t)) * n));
  return t;
}

color difuso(color ip, double kd, vec lightDir, vec normal, color corObjeto) {
  color retorno;
  double prodEscalar = lightDir.dot(normal);
  double aux = kd * prodEscalar;
  retorno.x = ip.x * aux * corObjeto.x;
  retorno.y = ip.y * aux * corObjeto.y;
  retorno.z = ip.z * aux * corObjeto.z;

  return retorno;
};

double rand01() {
  random_device rd;
  return ((double)rd() / (double)rd.max());
}


Intersection closestObject(Ray ray, Scene scene) {
  // TODO
  Intersection out;
  Object current;
  vector<Face> currentFaces;
  Face face;
  double closestDist = INT_MAX;

  for (int i = 0; i < objetos.size(); i++) {
    current = objetos.at(i);
    currentFaces = current.faces;
    for (int j = 0; j < currentFaces.size(); j++) {
      Face f = currentFaces.at(j);
      point* x = f.v1;
      double v0[3] = {x->x, x->y, x->z};
      point* y = f.v2;
      double v1[3] = {y->x, y->y, y->z};
      point* z = f.v3;
      double v2[3] = {z->x, z->y, z->z};

      double p[3] = {ray.posicao.x, ray.posicao.y, ray.posicao.z};
      double d[3] = {ray.direcao.x, ray.direcao.y, ray.direcao.z};
      Triangle tt({v0[0], v0[1], v0[2]},{v1[0], v1[1], v1[2]} , {v2[0], v2[1], v2[2]});
      double t = tt.intersect({{p[0], p[1], p[2]}, {d[0], d[1], d[2]}});
      //double t = rayIntersectsTriangle(p, d, v0, v1, v2);
      if (t > 0 && t < closestDist) {
        closestDist = t;
        out.objeto = current;
        out.p.x = ray.direcao.x * t + ray.posicao.x;
        out.p.y = ray.direcao.y * t + ray.posicao.y;
        out.p.z = ray.direcao.z * t + ray.posicao.z;
        out.normal = current.point_normal(out.p, f);
      }
    }
  }

  if (closestDist != INT_MAX) {
    out.hit = true;
  } else {
    out.hit = false;
  }

  return out;
}

Ray cameraRay(int x, int y, Screen jan, point o) {
  /*
   * Trace a ray from the camera to the pixel (x,y) on the window
   */

  Ray output;  // Output variable
  double xw, yw, zw, sizexw,
      sizeyw;            // Pixel location and window size on world coordinates
  vec direcao, posicao;  // Ray's direction and position

  sizexw = jan.x1 - jan.x0;
  sizeyw = jan.y0 - jan.y1;
  xw = ((double)x / jan.sizeX) * sizexw + jan.x0;
  yw = ((double)y / jan.sizeY) * sizeyw + jan.y1;
  zw = 0;

  posicao.x = o.x;
  posicao.y = o.y;
  posicao.z = o.z;

  direcao.x = xw - o.x;
  direcao.y = yw - o.y;
  direcao.z = zw - o.z;

  output.posicao = posicao;
  output.direcao = direcao;
  output.direcao = output.direcao.normalized();
  return output;
}

bool shadowRay(Ray ray, Scene scene) {
  bool retorno = false;
  Intersection intersection = closestObject(ray, scene);
  Object closest = intersection.objeto;
  if (closest.isLight) return false;

  if (intersection.hit) retorno = true;
  return retorno;
}

vec rotacionar(double angle, vec v, vec eixo) {
  double matrix[3][3];
  double co = cos(angle);
  double si = sin(angle);
  double a = eixo.x;
  double b = eixo.y;
  double c = eixo.z;
  matrix[0][0] = co + (1 - co) * pow(a, 2);
  matrix[0][1] = (1 - co) * a * b + si * c;
  matrix[0][2] = (1 - co) * a * c - si * b;
  matrix[1][0] = (1 - co) * b * a - si * c;
  matrix[1][1] = co + (1 - co) * pow(b, 2);
  matrix[1][2] = (1 - co) * b * c + si * a;
  matrix[2][0] = (1 - co) * a * c + si * b;
  matrix[2][1] = (1 - co) * b * c - si * a;
  matrix[2][2] = co + (1 - co) * pow(b, 2);
  vec out;
  out.x = matrix[0][0] * v.x + matrix[0][1] * v.y + matrix[0][2] * v.z;
  out.y = matrix[1][0] * v.x + matrix[1][1] * v.y + matrix[1][2] * v.z;
  out.z = matrix[2][0] * v.x + matrix[2][1] * v.y + matrix[2][2] * v.z;
  return out;
}

color trace_path(int depth, Ray ray, Scene scene, Light luz, int i, int j,
                 int nSample) {
  if (depth >= MAX_RECURSION_DEPTH) return {0.0, 0.0, 0.0};
  color output = {0.0,0.0,0.0};

  // --------------------------check intersections--------------------------
  // ray intersects triangle for each triangle. Get the one closest to the eye
  // that returns true.
  ccout "raio " << ray.posicao << " " << ray.direcao << ray.tamanho << '\n';
  Intersection intersection;
  intersection = closestObject(ray, scene);

  // if (depth == 0 && nSample == 0) {
  //   intersection = closestObject(ray, scene);
  //   bInt1[i][j] = intersection;
  // } else if (depth == 0 && nSample != 0) {
  //   intersection = bInt1[i][j];
  // } else {
  // }

  if (intersection.hit == false) {
    ccout "sem colisão " << scene.bg << '\n';
    return scene.bg;
  } else if (intersection.objeto.isLight) {
    return luz.col;
  }
  ccout "com colisão\n";

  Object closest = intersection.objeto;
  point inters = intersection.p;
  vec normal = intersection.normal;
  normal = normal.normalized();

  double bias = 1e-4;  // Use in order to avoid noise

  // Sortear um ponto de luz aleat�rio dentro do obj Luz
  int triangulo = rand() % 2;
  Face triLuz = objetos.at(0).faces.at(triangulo);
  double alpha = rand() % 100;
  double beta = rand() % 100;
  double gama = rand() % 100;
  double sum = alpha + beta + gama;
  alpha = alpha / sum;
  beta = beta / sum;
  gama = gama / sum;

  point* v1 = triLuz.v1;
  point* v2 = triLuz.v2;
  point* v3 = triLuz.v3;

  point lightRand;

  lightRand.x = alpha * v1->x + beta * v2->x + gama * v3->x;
  lightRand.y = v1->y;
  lightRand.z = alpha * v1->z + beta * v2->z + gama * v3->z;

  vec toLight = cena.light.pos - inters;
  toLight = toLight.normalized();


  //---------------------------------------------------------
  float kd = closest.kd, ks = closest.ks, kt = closest.kt;

  // ---------------------------color calculation----------------------------
  // Rambiente = Ia*kar
  float iA = scene.la;
  vec ambiente = iA * closest.ka * closest.col;

  // Rdifuso = Ip*kd(L.N)r
  color difusa = difuso(luz.col, closest.kd, toLight, normal, closest.col);

  // Respecular = Ip*ks*(R.V)^n
  vec rVetor = (2 * toLight.dot(normal) * normal) - toLight;
  rVetor = rVetor.normalized();
  vec vVetor = -1 * ray.direcao;
  vVetor = vVetor.normalized();
  color especular;
  float aux = pow(rVetor.dot(vVetor), closest.n);
  aux = luz.lp * closest.ks * aux;
  especular.x = luz.col.x * aux;
  especular.y = luz.col.y * aux;
  especular.z = luz.col.z * aux;

  // Shadow Ray
  Ray ray2;
  ray2.direcao = lightRand - inters;

  // Walk a little bit in the normal direction in order to avoid self
  // intersection
  vec dist = bias * normal;

  ray2.posicao.x = inters.x + dist.x;
  ray2.posicao.y = inters.y + dist.y;
  ray2.posicao.z = inters.z + dist.z;

  bool sombra = shadowRay(ray2, scene);

  ////Definindo o valor da cor local
  color corLocal;
  if (sombra) {
    corLocal = {0.0, 0.0, 0.0};
  } else {
    corLocal = difusa + color(ambiente) + especular;
  }

  // Ray cast from the recursion
  Ray novoRaio;
  novoRaio.posicao.x = inters.x;
  novoRaio.posicao.y = inters.y;
  novoRaio.posicao.z = inters.z;

  // -------------------------recursion for contribution from other
  // objects---------------------------------
  double ktot = kd + ks + kt;
  double r = rand01() * ktot;
  vec direcao = vec(), posicao = vec();
  if (r < kd) {
    // raio difuso
    float r1 = 2 * PI * rand01();  // random angle around
    float r2 = rand01();           // random distance from center
    float r2s = sqrt(r2);          // square root of distance from center

    vec w = normal;  // set first axis equal to normal
    vec v1;
    v1.x = 0;
    v1.y = 1;
    v1.z = 0;
    vec v2;
    v2.x = 1;
    v2.y = 0;
    v2.z = 0;
    vec u = fabs(w.x) > 0.1 ? v1 : v2;
    u = u.cross(w).normalized();  // second axis
    vec v = w.cross(u);           // final axis

    // random direction
    vec psi;
    psi.x = u.x * cos(r1) * r2s + v.x * sin(r1) * r2s + w.x * sqrt(1 - r2);
    psi.y = u.y * cos(r1) * r2s + v.y * sin(r1) * r2s + w.y * sqrt(1 - r2);
    psi.z = u.z * cos(r1) * r2s + v.z * sin(r1) * r2s + w.z * sqrt(1 - r2);

    psi = psi.normalized();

    direcao = psi;
  } else if (r < kd + ks) {
    // raio especular
    // direcao: R=2N(NL) - L
    direcao = (2 * normal.dot(toLight) * normal) - toLight;
  } else {
    // raio transmitido
    // TODO
    // objeto opaco? nenhuma cor transmitida
    // caso contrario... verificar refracao

    float cos = ray.direcao.dot(normal);
    float n1, n2;
    if (cos > 0) {
      // Inside the object
      n1 = closest.kt;  // coeficienteRefracao;
      n2 = 1;
      direcao = calcularRefracao(n1, n2, ray.direcao, -1 * normal);
      normal = -1 * normal;
    } else {
      n1 = 1;
      n2 = closest.kt;
      direcao = calcularRefracao(n1, n2, ray.direcao, normal);
    }

    vec dist2 = bias * normal;
    novoRaio.posicao.x = inters.x - dist2.x;
    novoRaio.posicao.y = inters.y - dist2.y;
    novoRaio.posicao.z = inters.z - dist2.y;
  }

  // Use the direction and position vector to make a ray

  novoRaio.direcao = direcao;
  novoRaio.direcao = (novoRaio.direcao).normalized();
  float cos_theta = direcao.dot(normal);

  color recursion =
      trace_path(depth + 1, novoRaio, scene, scene.light, i, j, nSample);
  
  float factor = max(max(kd, ks), kt);
  output = (recursion + corLocal) * factor;

  return output;
}

float clamp(float x) { return x < 0 ? 0 : x > 1 ? 1 : x; };

float tColor(float x) { return pow(clamp(x), 1 / 2.2); };

std::vector<std::vector<color>> render(Screen jan, Scene scene, point o,
                                       Light luz) {
  int xsize = jan.sizeX;  // width in pixels
  int ysize = jan.sizeY;  // height in pixels

  int nSamples = scene.npaths;  // number of color samples per pixel

  float count = 0, maxCount = xsize * ysize * nSamples;
  Ray ray;  // Camera to window variable, used for each different pixel
  std::vector<std::vector<color>> img = {};
  for (int i = 0; i < xsize; i++) {
    std::vector<color> temp = {};
    for (int j = 0; j < ysize; j++) {
      temp.push_back({0.0, 0.0, 0.0});
    }
    img.push_back(temp);
  }
   time_t startTime;
  time(&startTime);


  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      color sum, sample;

      sum = {0.0, 0.0, 0.0};
      ray = cameraRay(i, j, jan, o);
      for (int k = 0; k < nSamples; k++) {
        sample = trace_path(0, ray, scene, luz, i, j, k);
        sum = sum + sample;
      }
      //cout << "color  " << sum << '\n'; 
      color out = sum / nSamples;  
      color newOut = out / (out + cena.tonemapping);
      img[i][j] = newOut;
    }
    auto fds = img[i][0];
    //cout << "done " << i << " " << fds.x << " " << fds.y << " " << fds.z
    //     << '\n';
  }

  return img;
}

namespace ppm {
typedef std::vector<std::vector<color>> display;

void save2ppm(std::string path, display screen, std::pair<int, int> res) {
  int w = res.first, h = res.second;
  std::ofstream file(path, std::ios_base::binary);
  if (!file) return;

  std::stringstream content;
  content << "P6\n"
          << "#\n#made by casg\n#\n"
          << res.first << ' ' << res.second << '\n'
          << "255\n";
  for (auto row : screen)
    for (auto pixel : row)
      content << (char)(((int)pixel.x * 255) % 256)
              << (char)(((int)pixel.y * 255) % 256)
              << (char)(((int)pixel.z * 255) % 256);
  file << content.rdbuf();
  file.close();
};
};  // namespace ppm

int main(int argc, char** argv) {
  // Lendo arquivo sdl que descreve a cena utilizada e calculando a normal ap�s
  loadScene("cornellroom.sdl", cena);

  // Loading camera paramenters
  olho = cena.eye;

  // Loading view paramenters
  janela = cena.screen;

  //cena.light.pos = {0, 3.8360, -24};

  // Loading objects of the scene
  objetos = cena.objects;

  // window_height = janela.sizeY;
  // window_width = <janela.sizeX>;

  // Carregando os objetos no vector de objetos
  for (int i = 0; i < objetos.size(); i++) {
    char realPath[100] = "";
    strcat(realPath, objetos.at(i).path.c_str());
    loadObject(realPath, objetos.at(i));
    objetos.at(i).compute_vertices_normals();
  }
  for(auto i: objetos){
    cout << i.col <<" "<< i.ka << " "<< i.kd << " "<< i.ks << " " << i.kt << '\n';
    for(auto v: i.vertices)
      ccout v << " || ";
    ccout '\n';
  }

  auto kk = render(janela, cena, olho, luz);
  ppm::save2ppm("uau.ppm", kk, {janela.sizeY, janela.sizeX});

  return 0;
}