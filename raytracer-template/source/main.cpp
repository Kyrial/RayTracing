//////////////////////////////////////////////////////////////////////////////
//
//  --- main.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
//sp�cifique a W10 !!!!

#include "common.h"
#include "SourcePath.h"

#include <math.h>

using namespace Angel;

typedef vec4  color4;
typedef vec4  point4;


//Scene variables
enum{_SPHERE, _SQUARE, _BOX};
int scene = _SPHERE; //Simple sphere, square or cornell box
std::vector < Object * > sceneObjects;
point4 lightPosition;
color4 lightColor;
point4 cameraPosition;

//Recursion depth for raytracer
int maxDepth = 10;

void initGL();

namespace GLState {
  int window_width, window_height;

  bool render_line;

  std::vector < GLuint > objectVao;
  std::vector < GLuint > objectBuffer;

  GLuint vPosition, vNormal, vTexCoord;

  GLuint program;

  // Model-view and projection matrices uniform location
  GLuint  ModelView, ModelViewLight, NormalMatrix, Projection;

  //==========Trackball Variables==========
  static float curquat[4],lastquat[4];
  /* current transformation matrix */
  static float curmat[4][4];
  mat4 curmat_a;
  /* actual operation  */
  static int scaling;
  static int moving;
  static int panning;
  /* starting "moving" coordinates */
  static int beginx, beginy;
  /* ortho */
  float ortho_x, ortho_y;
  /* current scale factor */
  static float scalefactor;

  mat4  projection;
  mat4 sceneModelView;

  color4 light_ambient;
  color4 light_diffuse;
  color4 light_specular;

};

/* ------------------------------------------------------- */
/* -- PNG receptor class for use with pngdecode library -- */
class rayTraceReceptor : public cmps3120::png_receptor
{
private:
  const unsigned char *buffer;
  unsigned int width;
  unsigned int height;
  int channels;

public:
  rayTraceReceptor(const unsigned char *use_buffer,
                   unsigned int width,
                   unsigned int height,
                   int channels){
    this->buffer = use_buffer;
    this->width = width;
    this->height = height;
    this->channels = channels;
  }
  cmps3120::png_header get_header(){
    cmps3120::png_header header;
    header.width = width;
    header.height = height;
    header.bit_depth = 8;
    switch (channels)
    {
      case 1:
      header.color_type = cmps3120::PNG_GRAYSCALE;break;
      case 2:
      header.color_type = cmps3120::PNG_GRAYSCALE_ALPHA;break;
      case 3:
      header.color_type = cmps3120::PNG_RGB;break;
      default:
      header.color_type = cmps3120::PNG_RGBA;break;
    }
    return header;
  }
  cmps3120::png_pixel get_pixel(unsigned int x, unsigned int y, unsigned int level){
    cmps3120::png_pixel pixel;
    unsigned int idx = y*width+x;
    /* pngdecode wants 16-bit color values */
    pixel.r = buffer[4*idx]*257;
    pixel.g = buffer[4*idx+1]*257;
    pixel.b = buffer[4*idx+2]*257;
    pixel.a = buffer[4*idx+3]*257;
    return pixel;
  }
};

/* -------------------------------------------------------------------------- */
/* ----------------------  Write Image to Disk  ----------------------------- */
bool write_image(const char* filename, const unsigned char *Src,
                 int Width, int Height, int channels){
  cmps3120::png_encoder the_encoder;
  cmps3120::png_error result;
  rayTraceReceptor image(Src,Width,Height,channels);
  the_encoder.set_receptor(&image);
  result = the_encoder.write_file(filename);
  if (result == cmps3120::PNG_DONE)
    std::cerr << "finished writing "<<filename<<"."<<std::endl;
  else
    std::cerr << "write to "<<filename<<" returned error code "<<result<<"."<<std::endl;
  return result==cmps3120::PNG_DONE;
}


/* -------------------------------------------------------------------------- */
/* -------- Given OpenGL matrices find ray in world coordinates of ---------- */
/* -------- window position x,y --------------------------------------------- */
std::vector < vec4 > findRay(GLdouble x, GLdouble y){

  y = GLState::window_height-y;

  int viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  GLdouble modelViewMatrix[16];
  GLdouble projectionMatrix[16];
  for(unsigned int i=0; i < 4; i++){
    for(unsigned int j=0; j < 4; j++){
      modelViewMatrix[j*4+i]  =  GLState::sceneModelView[i][j];
      projectionMatrix[j*4+i] =  GLState::projection[i][j];
    }
  }


  GLdouble nearPlaneLocation[3];
  _gluUnProject(x, y, 0.0, modelViewMatrix, projectionMatrix,
                viewport, &nearPlaneLocation[0], &nearPlaneLocation[1],
                &nearPlaneLocation[2]);

  GLdouble farPlaneLocation[3];
  _gluUnProject(x, y, 1.0, modelViewMatrix, projectionMatrix,
                viewport, &farPlaneLocation[0], &farPlaneLocation[1],
                &farPlaneLocation[2]);


  vec4 ray_origin = vec4(nearPlaneLocation[0], nearPlaneLocation[1], nearPlaneLocation[2], 1.0);
  vec3 temp = vec3(farPlaneLocation[0]-nearPlaneLocation[0],
                   farPlaneLocation[1]-nearPlaneLocation[1],
                   farPlaneLocation[2]-nearPlaneLocation[2]);
  temp = normalize(temp);
  vec4 ray_dir = vec4(temp.x, temp.y, temp.z, 0.0);

  std::vector < vec4 > result(2);
  result[0] = ray_origin;
  result[1] = ray_dir;

  return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
bool intersectionSort(Object::IntersectionValues i, Object::IntersectionValues j){
  return (i.t < j.t);
}

/* -------------------------------------------------------------------------- */
/* ---------  Some debugging code: cast Ray = p0 + t*dir  ------------------- */
/* ---------  and print out what it hits =                ------------------- */
void castRayDebug(vec4 p0, vec4 dir){
    
  std::vector < Object::IntersectionValues > intersections;

  for(unsigned int i=0; i < sceneObjects.size(); i++){
    intersections.push_back(sceneObjects[i]->intersect(p0, dir));
    intersections[intersections.size()-1].ID_ = i;
  }

  for(unsigned int i=0; i < intersections.size(); i++){
    if(intersections[i].t != std::numeric_limits< double >::infinity()){
      std::cout << "Hit " << intersections[i].name << " " << intersections[i].ID_ << "\n";
      std::cout << "P: " <<  intersections[i].P << "\n";
      std::cout << "N: " <<  intersections[i].N << "\n";
      vec4 L = lightPosition-intersections[i].P;
      L  = normalize(L);
      std::cout << "L: " << L << "\n";
    }
  }

}

/* -------------------------------------------------------------------------- */
double shadowFeeler(vec4 p0, Object* object, vec4 E, double dist) {
    int precision = 4;
    double inShadow = 1.0 +  (1.0 / ((precision * 2) * 3 + 1) * precision*2);
    double valMax = inShadow;
    double espacement = 1.4;
    
    //TODO: Shadow code here


    float minTransparance = 1;
    for (unsigned int i = 0; i < sceneObjects.size(); i++) {
        Object::IntersectionValues inter = (sceneObjects[i]->intersect(p0, normalize(lightPosition - p0)));
        if ((/*(sceneObjects[i] == object && inter.t >= EPSILON/2 ) ||*/ sceneObjects[i] != object) && inter.t < dist) {
            minTransparance = minTransparance-(1 -max(0, sceneObjects[i]->shadingValues.Kt - 0.5));
            if (minTransparance <= 0) 
                break;
        }
    }
    
    inShadow = inShadow - (valMax / ((precision * 2) * 3 + 1)) *max(0, 1 - minTransparance);
        
    

     minTransparance = 1;
  /*  if (inShadow ==1 + (1.0 / ((precision * 2) * 3 + 1) * precision * 2)/4)
        return 1.0;*/

    for (int x = -precision; x <= precision; x++)
        if (x != 0) {
            double valX = (double(x) / precision)* espacement;
            for (unsigned int i = 0; i < sceneObjects.size(); i++) {
                Object::IntersectionValues inter = (sceneObjects[i]->intersect(p0, normalize((lightPosition - vec4(valX, 0., 0.)) - p0)));
                if (((sceneObjects[i] == object && inter.t >= EPSILON/2) || sceneObjects[i] != object) && inter.t < dist) {
                    minTransparance = minTransparance -(1-max(0, sceneObjects[i]->shadingValues.Kt - 0.5));
                    if (minTransparance <= 0)
                        break;
                }
            }
        }
    inShadow = inShadow - (valMax / ((precision * 2) * 3 + 1)) * max(0, 1 - minTransparance);


     minTransparance = 1;
    for (int y = -precision; y <= precision; y++)
        if (y != 0) {
            double valY = (double(y) / precision) * espacement;
            for (unsigned int i = 0; i < sceneObjects.size(); i++) {
                Object::IntersectionValues inter = (sceneObjects[i]->intersect(p0, normalize((lightPosition - vec4(0., valY, 0.)) - p0)));
                if ((/*(sceneObjects[i] == object && inter.t >= EPSILON/2) ||*/ sceneObjects[i] != object) && inter.t < dist) {
                    minTransparance = minTransparance -(1 - max(0, sceneObjects[i]->shadingValues.Kt - 0.5));
                    if (minTransparance <= 0)
                        break;
                }
            }
        }
    inShadow = inShadow - (valMax / ((precision * 2) * 3 + 1)) *  max(0, 1 - minTransparance);


     minTransparance = 1;
    for (int z = -precision; z <= precision; z++)
        if (z != 0) {
            double valZ = (double(z) / precision) * espacement;

            for (unsigned int i = 0; i < sceneObjects.size(); i++) {
                Object::IntersectionValues inter = (sceneObjects[i]->intersect(p0, normalize((lightPosition - vec4(0., 0., valZ)) - p0)));
                if ((/*(sceneObjects[i] == object && inter.t >= EPSILON/2) ||*/ sceneObjects[i] != object) && inter.t < dist) {
                    minTransparance = minTransparance - (1 - max(0, sceneObjects[i]->shadingValues.Kt -0.5));
                    if (minTransparance <= 0)
                        break;
                }
            }
        }
    inShadow = inShadow - (valMax / ((precision * 2) * 3 + 1)) * max(0, 1 - minTransparance);

        
    return max( 0 , min( 1, inShadow));
}
/*
bool shadowFeeler(vec4 p0, Object *object, vec4 E, double dist){
  bool inShadow = false;

  //TODO: Shadow code here

  for (unsigned int i = 0; i < sceneObjects.size(); i++) {
      Object::IntersectionValues inter = (sceneObjects[i]->intersect(p0, normalize(E)));
      if (sceneObjects[i] != object  && inter.t < dist)
          return true;
  }

  return inShadow;
}*/


double dotMajorer(vec4 L, vec4 N) {
    double value = dot(L, N);
    if (value > 1)
        return 1;
    else
        return value;
}




/* -------------------------------------------------------------------------- */
/* ----------  cast Ray = p0 + t*dir and intersect with sphere      --------- */
/* ----------  return color, right now shading is approx based      --------- */
/* ----------  depth                                                --------- */
vec4 castRay(vec4 p0, vec4 E, Object *lastHitObject, int depth){
  vec4 color = vec4(0.,0.0,0.0,0.);
 

  if(depth > maxDepth){ return color; }

  std::vector < Object::IntersectionValues > intersections;
  
  double minT = std::numeric_limits< double >::infinity();
  int plusProcheObject=0;

  for(unsigned int i=0; i < sceneObjects.size(); i++){
    intersections.push_back(sceneObjects[i]->intersect(p0, E));
    intersections[intersections.size()-1].ID_ = i;

    if (intersections[intersections.size() - 1].t != std::numeric_limits< double >::infinity() && /*!(sceneObjects[intersections[intersections.size() - 1].ID_] == lastHitObject) &&*/ intersections[intersections.size() - 1].t < minT) {
        minT = intersections[intersections.size() - 1].t;
        plusProcheObject = i;
        }
    }
 



  //for(unsigned int i=0; i < intersections.size(); i++){
      
     // if (intersections[i].t != std::numeric_limits< double >::infinity()  && !(sceneObjects[intersections[i].ID_] == lastHitObject)) {
          //std::cout << "t et minT " << intersections[i].t << " et  " << minT << std::endl;
          //if (intersections[i].t < minT ) {
  if ( minT < std::numeric_limits< double >::infinity()){
              //minT = intersections[i].t;
              //  std::cout << "Hit " << intersections[i].name << " " << intersections[i].ID_ << "\n";
              //  std::cout << "P: " <<  intersections[i].P << "\n";
              //  std::cout << "N: " <<  intersections[i].N << "\n";
              vec4 L = lightPosition - intersections[plusProcheObject].P;
              vec4 V = cameraPosition - intersections[plusProcheObject].P;
              vec4 N = intersections[plusProcheObject].N;
             
             


             /* if (shadowFeeler(intersections[plusProcheObject].P, sceneObjects[plusProcheObject], L, length(L))) {
                  color = vec4(0., 0.0, 0.0, 1.);
                  
              }*/
             // else {
              double ombre = shadowFeeler(intersections[plusProcheObject].P, sceneObjects[plusProcheObject], L, length(L));
            // double ombre = 1;

                  // L.w = 0;
                 //  V.w = 0;
                 //  N.w = 0;
              N.w = 0; L.w = 0; V.w = 0;
                  N = normalize(N),
                      L = normalize(L);
                  V = normalize(V);
                 


                  //   std::cout << "L: " << L << "\n";
                  //   std::cout << "V: " << V << "\n";
                  //   std::cout << "N: " << N << "\n";
                   //color = sceneObjects[intersections[i].ID_]->shadingValues.color;
                   //modele de Phong ->   (ambiante + diffuse + sp�culaire )
                   //color = phong(sceneObjects[intersections[i].ID_]);
                  Object* object = sceneObjects[intersections[plusProcheObject].ID_];
                  //if (dot(V, N) > 0 && dot(V, intersections[i].P) >0 && dot(V, L)>0) {
                  vec4 material_ambient(object->shadingValues.color.x * object->shadingValues.Ka,
                      object->shadingValues.color.y * object->shadingValues.Ka,
                      object->shadingValues.color.z * object->shadingValues.Ka,
                      object->shadingValues.color.w * object->shadingValues.Ka);

                  vec4 ambient_product = GLState::light_ambient * material_ambient;


                  vec4 material_diffuse(object->shadingValues.color.x * object->shadingValues.Kd,
                      object->shadingValues.color.y * object->shadingValues.Kd,
                      object->shadingValues.color.z * object->shadingValues.Kd, 1);
                     // object->shadingValues.color.w * object->shadingValues.Kd);


                  color4 diffuse_product = GLState::light_diffuse * material_diffuse * min(1, max(dot(N, L), 0));
                  //if (min(1, max(dot(N, L), 0))==0)
                      diffuse_product.w = 1;

                  //color4 diffuse_product = GLState::light_diffuse * material_diffuse * sin(angle(L,N));
                 // std::cout << "test  " << diffuse_product << "\n ";

                  vec4 material_specular(object->shadingValues.Ks,
                      object->shadingValues.Ks,
                      object->shadingValues.Ks,
                      object->shadingValues.Ks);
                  vec4 R = normalize(2 * (min(1, max(dot(N, L), 0))) * N - L);
                  // vec4 R = normalize(2 * ((angle(N,L))) * N - L);
                   //vec4 R = normalize(reflect(L, N));
                  // std::cout << "test  " << R << std::endl;
                  vec4 specular_product = GLState::light_specular * material_specular *
                      pow((min(1, max(dot(R, V), 0))), object->shadingValues.Kn );

                  // +object->shadingValues.color;
                  //color = ambient_product+ diffuse_product;
                  color = (ambient_product + diffuse_product + specular_product);// +object->shadingValues.color;
                  color.x = color.x* ombre;
                  color.y = color.y * ombre;
                  color.z = color.z * ombre;
                  //std::cout << "test  " << ambient_product << " et " << diffuse_product << " et " << specular_product <<"et"<< dot(L, N) << "\n";
           //   }
            //  else
            //      color = vec4(0., 0.0, 0.0, 1.);

                  //effet miroir
                  if (object->shadingValues.Ks > 0 && depth < maxDepth){

                    vec4 miroir = normalize(2 * (min(1, max(dot(N, V), 0))) * N - V);
                     // vec4 miroir = normalize(-V - (2 * (min(1, max(dot(N, -V), 0))) * N));
                      color4 colorTemp = castRay(intersections[plusProcheObject].P, miroir, object, depth + 1);
                      colorTemp = colorTemp * max(object->shadingValues.Ks-0.1,0);
                     // color = colorTemp + color;
                  /*    color.x = max(color.x, colorTemp.x);
                      color.y = max(color.y, colorTemp.y);
                      color.z = max(color.z, colorTemp.z);
                      color.w = max(color.w, colorTemp.w);*/

                  /*    color.x = (color.x + colorTemp.x) / 2;
                      color.y = (color.y + colorTemp.y) / 2;
                      color.z = (color.z + colorTemp.z) / 2;
                      color.w = (color.w + colorTemp.w) / 2;*/
                     // color = (colorTemp + color/2);

                      color = color * max(1-(object->shadingValues.Ks - 0.1), 0.1);
                      //color =( color*(1+ 1*object->shadingValues.Ks) + colorTemp) / 2;
                      color = (color + colorTemp);
                      color.w = 1;

                  }

                  //effet transparance
                  if (object->shadingValues.Kt > 0 && depth <maxDepth) {

                      //  vec4 transparance = normalize(2 * (min(1, max(dot(-N, V), 0))) * (-N) - V);

                      //  vec4 transmission = normalize(min(max((dot(-N, V) * object->shadingValues.Kr), 1), 0) * (-N) - V);
                      
                      double n1, n2;
                      if (lastHitObject == object) {
                         n1 = object->shadingValues.Kr;
                          n2 = 1.0;
                     }
                     else {
                          n2 = object->shadingValues.Kr;
                          n1 = 1.0;
                         
                      }
                      
                      

                      double n = n1 / n2;
  
                     //double cosTeta = dot(N, V) / (length(N) * length(V));
                      double cosTeta = dot(N, E) / (length(E) * length(V));
                     if (cosTeta < 0) cosTeta = -cosTeta;
                     else
                         N = -N;

                    
                     // double c2 = sqrt(1- pow(1/ object->shadingValues.Kr,2) * pow(1 - cos(dot(E, N)),2));
                     double c2 = sqrt(1 - n*n*(1-cosTeta*cosTeta));

                      vec4 transmission = normalize(n * E + (n * cosTeta - c2) * N);
                      
                    
                     
                  
                      //////////////////
             //         double cosTeta = dot(N, V) / (length(N) * length(V));
              //        if (cosTeta < 0) {
             //             cosTeta = -cosTeta; N = -N;
               //       }
                      //double teta = acos(dot(N, V) / (length(N) * length(V)));
             //         double teta = acos(cosTeta);
             //        double teta_t = asin(n * sin(teta));
                //      vec4 M = normalize((N * dot(N, V) - V));
                    
                      ///////////////
                      
                      
                    //  vec4 transmission = normalize(sin(teta_t) * M - cos(teta_t) *N);

                    // vec4 transmission = normalize((teta_t- teta)*(-N)-V);
                    //   vec4 transmission = normalize(cos(teta_t-teta) * (-N) - V);

                       // sin(teta_t) = (n1 / n2) � sin(teta_i)

                      //  M = (N cos teta_i - V) / normalize( sin teta_i)


                        //transmission = sin teta_t M - cos  teta_t N
                     /*
                     double teta = acos(dot(V, N) / (length(N) * length(V)));


                      double racine = 1 - ((n * n) - (sin(teta) * sin(teta)));
                      vec4 transmission;
                      if (racine < 0)
                          //transmission = normalize(n * (-V) + (n * (-cos(teta)) + sqrt(-racine)) * N);
                          transmission = 0;
                      else
                           transmission = normalize(n * (-V) + (n * cos(teta) - sqrt(racine)) * N);
                      */
                          
                      color4 colorTemp = castRay(intersections[plusProcheObject].P, transmission, object, depth + 1);
                      //  colorTemp = colorTemp * object->shadingValues.Kt;
                      colorTemp = colorTemp * max(object->shadingValues.Kt - 0.1, 0);
                      // color = colorTemp + color;
                 //   color.x = max(color.x, colorTemp.x);
                 //   color.y = max(color.y, colorTemp.y);
                 //   color.z = max(color.z, colorTemp.z);
                 //   color.w = max(color.w, colorTemp.w);
     
                 // 
                 // }
                      color = color * max(1 - (object->shadingValues.Kt - 0.1), 0.1);
                      color = (color + colorTemp);
                      color.w = 1;
                  }
  


              //}
          }
    //  }
  //}

  

  //TODO: Raytracing code here

  return color;

}


    

/*
vec4 phong(Object* object) {
    //ambiante -> Ia = Isa * Ka
    vec4 material_ambient(object->shadingValues.color.x * object->shadingValues.Ka,
                            object->shadingValues.color.y * object->shadingValues.Ka,
                            object->shadingValues.color.z * object->shadingValues.Ka,
                            1.0);
    vec4 ambient_product = GLState::light_ambient * material_ambient;

    return ambient_product;
}*/



/* -------------------------------------------------------------------------- */
/* ------------  Ray trace our scene.  Output color to image and    --------- */
/* -----------   Output color to image and save to disk             --------- */
void rayTrace(){

  unsigned char *buffer = new unsigned char[GLState::window_width*GLState::window_height*4];
  printf("miaou \n");
  for(unsigned int i=0; i < GLState::window_width; i++){
    for(unsigned int j=0; j < GLState::window_height; j++){

      int idx = j*GLState::window_width+i;
      std::vector < vec4 > ray_o_dir = findRay(i,j);
      vec4 color = castRay(ray_o_dir[0], vec4(ray_o_dir[1].x, ray_o_dir[1].y, ray_o_dir[1].z,0.0), NULL, 0);
      buffer[4*idx]   = min(color.x*255, 255);
      buffer[4*idx+1] = min(color.y*255,255);
      buffer[4*idx+2] = min(color.z*255,255);
      buffer[4*idx+3] = min(color.w*255,255);
    }
  }

  write_image("output.png", buffer, GLState::window_width, GLState::window_height, 4);

  delete[] buffer;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void error_callback(int error, const char* description)
{
  fprintf(stderr, "Error: %s\n", description);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initCornellBox(){
  cameraPosition = point4( 0.0, 0.0, 6.0, 1.0 );
  lightPosition = point4( 0.0, 1.5, 0.0, 1.0 );
  lightColor = color4( 1.0, 1.0, 1.0, 1.0);

  sceneObjects.clear();

  { //Back Wall
    sceneObjects.push_back(new Square("Back Wall", Translate(0.0, 0.0, -2.0)*Scale(2.0,2.0,1.0)));
    Object::ShadingValues _shadingValues;
    _shadingValues.color = vec4(0.0,1.0,0.0,1.0);
    _shadingValues.Ka = 0.0;
    _shadingValues.Kd = 1.0;
    _shadingValues.Ks = 0.;
    _shadingValues.Kn = 16.0;
    _shadingValues.Kt = 0.0;
    _shadingValues.Kr = 0.0;
    sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
    sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }

  { //Left Wall
    sceneObjects.push_back(new Square("Left Wall", RotateY(90)*Translate(0.0, 0.0, -2.0)*Scale(2.0,2.0,1.0)));
    Object::ShadingValues _shadingValues;
    _shadingValues.color = vec4(1.0,0.0,0.0,1.0);
    _shadingValues.Ka = 0.0;
    _shadingValues.Kd = 1.0;
    _shadingValues.Ks = 0.0;
    _shadingValues.Kn = 16.0;
    _shadingValues.Kt = 0.0;
    _shadingValues.Kr = 0.0;
    sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
    sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }

  { //Right Wall
    sceneObjects.push_back(new Square("Right Wall", RotateY(-90)*Translate(0.0, 0.0, -2.0)*Scale(2.0, 2.0, 1.0 )));
    Object::ShadingValues _shadingValues;
    _shadingValues.color = vec4(0.5,0.0,0.5,1.0);
    _shadingValues.Ka = 0.0;
    _shadingValues.Kd = 1.0;
    _shadingValues.Ks = 0.0;
    _shadingValues.Kn = 16.0;
    _shadingValues.Kt = 0.0;
    _shadingValues.Kr = 0.0;
    sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
    sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }

  { //Floor
    sceneObjects.push_back(new Square("Floor", RotateX(-90)*Translate(0.0, 0.0, -2.0)*Scale(2.0, 2.0, 1.0)));
    Object::ShadingValues _shadingValues;
    _shadingValues.color = vec4(1.0,.5,.5,1.0);
    _shadingValues.Ka = 0.0;
    _shadingValues.Kd = 1.0;
    _shadingValues.Ks = 0.0;
    _shadingValues.Kn = 16.0;
    _shadingValues.Kt = 0.0;
    _shadingValues.Kr = 0.0;
    sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
    sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }

  { //Ceiling
    sceneObjects.push_back(new Square("Ceiling", RotateX(90)*Translate(0.0, 0.0, -2.0)*Scale(2.0, 2.0, 1.0)));
    Object::ShadingValues _shadingValues;
    _shadingValues.color = vec4(1.0,1.0,1.0,1.0);
    _shadingValues.Ka = 0.0;
    _shadingValues.Kd = 1.0;
    _shadingValues.Ks = 0.0;
    _shadingValues.Kn = 16.0;
    _shadingValues.Kt = 0.0;
    _shadingValues.Kr = 0.0;
    sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
    sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }

  { //Front Wall
    sceneObjects.push_back(new Square("Front Wall",RotateY(180)*Translate(0.0, 0.0, -2.0)*Scale(2.0, 2.0, 1.0)));
    Object::ShadingValues _shadingValues;
    _shadingValues.color = vec4(1.0,1.0,1.0,1.0);
    _shadingValues.Ka = 0.0;
    _shadingValues.Kd = 1.0;
    _shadingValues.Ks = 0.0;
    _shadingValues.Kn = 16.0;
    _shadingValues.Kt = 0.0;
    _shadingValues.Kr = 0.0;
    sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
    sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }


  {
  sceneObjects.push_back(new Sphere("Glass sphere", vec3(1.0, -1.2, 0.5),0.75));
 // sceneObjects.push_back(new Sphere("Glass sphere", vec3(1.0, 1, -0.5), 0.75));
  Object::ShadingValues _shadingValues;
  _shadingValues.color = vec4(1.0,0.0,0.0,1.0);
  _shadingValues.Ka = 0.0;
  _shadingValues.Kd = 1.0; // 0
  _shadingValues.Ks = 0.0;
  _shadingValues.Kn = 16.0;
  _shadingValues.Kt = 1.0;
  _shadingValues.Kr = 1.4;
  sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
  sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }

  {
  sceneObjects.push_back(new Sphere("Mirrored Sphere", vec3(-1.0, -1.2, 0.5), 0.75));
  //sceneObjects.push_back(new Sphere("Mirrored Sphere", vec3(-1.0, 1, -0.5),0.75));
  Object::ShadingValues _shadingValues;
  _shadingValues.color = vec4(1.0,1.0,1.0,1.0);
  _shadingValues.Ka = 0.0;
  _shadingValues.Kd = 0.0; //0
  _shadingValues.Ks = 1.0; //1
  _shadingValues.Kn = 16.0;
  _shadingValues.Kt = 0.0; 
  _shadingValues.Kr = 0.0; //valeur initial: 0
  sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
  sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initUnitSphere(){

    //std::cout << "test  " << GLState::light_ambient << "\n";

  cameraPosition = point4( 0.0, 0.0, 3.0, 1.0 );
  lightPosition = point4( 0.0, 0.0, 4.0, 1.0 );
  lightColor = color4( 1.0, 1.0, 1.0, 1.0);
  sceneObjects.clear();

  {
      //la sphere
  sceneObjects.push_back(new Sphere("Diffuse sphere rouge", vec3(-1., 0., 0.)));
  Object::ShadingValues _shadingValues;
  _shadingValues.color = vec4(1.0,0.0,0.0,1.0);
  _shadingValues.Ka = 0.0;
  _shadingValues.Kd = 1.0;
  _shadingValues.Ks = 0.0;
  _shadingValues.Kn = 16.0;
  _shadingValues.Kt = 0.0;
  _shadingValues.Kr = 0.0;
  sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
  sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }
  {
      //la sphere
      sceneObjects.push_back(new Sphere("Diffuse sphere verte", vec3(1., 0., 0.)));
      Object::ShadingValues _shadingValues;
      _shadingValues.color = vec4(0.0, 1.0, 0.2, 1.0);
      _shadingValues.Ka = 0.1;
      _shadingValues.Kd = 1.0;
      _shadingValues.Ks = 0.0;
      _shadingValues.Kn = 16.0;
      _shadingValues.Kt = 0.0;
      _shadingValues.Kr = 0.0;
      sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
      sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
  }





}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initUnitSquare(){
  cameraPosition = point4( 0.0, 0.0, 3.0, 1.0 );
  lightPosition = point4( 0.0, 0.0, 4.0, 1.0 );
  lightColor = color4( 1.0, 1.0, 1.0, 1.0);

  sceneObjects.clear();

  { //Back Wall
    sceneObjects.push_back(new Square("Unit Square"));
    Object::ShadingValues _shadingValues;
    _shadingValues.color = vec4(1.0,0.5,0.3,1.0);
    _shadingValues.Ka = 0.0;
    _shadingValues.Kd = 1.0;
    _shadingValues.Ks = 0.0;
    _shadingValues.Kn = 16.0;
    _shadingValues.Kt = 0.0;
    _shadingValues.Kr = 0.0;
    sceneObjects[sceneObjects.size()-1]->setShadingValues(_shadingValues);
    sceneObjects[sceneObjects.size()-1]->setModelView(mat4());
  }

}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    if (key == GLFW_KEY_1 && action == GLFW_PRESS) {

        if (scene != _SPHERE) {
            initUnitSphere();
            initGL();
            scene = _SPHERE;
        }

    }
    if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
        if (scene != _SQUARE) {
            initUnitSquare();
            initGL();
            scene = _SQUARE;
        }
    }
    if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
        if (scene != _BOX) {
            initCornellBox();
            initGL();
            scene = _BOX;
        }
    }


    if (key == GLFW_KEY_R && action == GLFW_PRESS)
        rayTrace();
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void mouseClick(GLFWwindow* window, int button, int action, int mods){

  if (GLFW_RELEASE == action){
    GLState::moving=GLState::scaling=GLState::panning=false;
    return;
  }

  if( mods & GLFW_MOD_SHIFT){
    GLState::scaling=true;
  }else if( mods & GLFW_MOD_ALT ){
    GLState::panning=true;
  }else{
    GLState::moving=true;
    TrackBall::trackball(GLState::lastquat, 0, 0, 0, 0);
  }

  double xpos, ypos;
  glfwGetCursorPos(window, &xpos, &ypos);
  GLState::beginx = xpos; GLState::beginy = ypos;

  std::vector < vec4 > ray_o_dir = findRay(xpos, ypos);
  castRayDebug(ray_o_dir[0], vec4(ray_o_dir[1].x, ray_o_dir[1].y, ray_o_dir[1].z,0.0));

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void mouseMove(GLFWwindow* window, double x, double y){

  int W, H;
  glfwGetFramebufferSize(window, &W, &H);


  float dx=(x-GLState::beginx)/(float)W;
  float dy=(GLState::beginy-y)/(float)H;

  if (GLState::panning)
    {
    GLState::ortho_x  +=dx;
    GLState::ortho_y  +=dy;

    GLState::beginx = x; GLState::beginy = y;
    return;
    }
  else if (GLState::scaling)
    {
    GLState::scalefactor *= (1.0f+dx);

    GLState::beginx = x;GLState::beginy = y;
    return;
    }
  else if (GLState::moving)
    {
    TrackBall::trackball(GLState::lastquat,
                         (2.0f * GLState::beginx - W) / W,
                         (H - 2.0f * GLState::beginy) / H,
                         (2.0f * x - W) / W,
                         (H - 2.0f * y) / H
                         );

    TrackBall::add_quats(GLState::lastquat, GLState::curquat, GLState::curquat);
    TrackBall::build_rotmatrix(GLState::curmat, GLState::curquat);

    GLState::beginx = x;GLState::beginy = y;
    return;
    }
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initGL(){

  GLState::light_ambient  = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0 );
  GLState::light_diffuse  = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0 );
  GLState::light_specular = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0 );


  std::string vshader = source_path + "/shaders/vshader.glsl";
  std::string fshader = source_path + "/shaders/fshader.glsl";

  GLchar* vertex_shader_source = readShaderSource(vshader.c_str());
  GLchar* fragment_shader_source = readShaderSource(fshader.c_str());

  GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_shader, 1, (const GLchar**) &vertex_shader_source, NULL);
  glCompileShader(vertex_shader);
  check_shader_compilation(vshader, vertex_shader);

  GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_shader, 1, (const GLchar**) &fragment_shader_source, NULL);
  glCompileShader(fragment_shader);
  check_shader_compilation(fshader, fragment_shader);

  GLState::program = glCreateProgram();
  glAttachShader(GLState::program, vertex_shader);
  glAttachShader(GLState::program, fragment_shader);

  glLinkProgram(GLState::program);
  check_program_link(GLState::program);

  glUseProgram(GLState::program);

  glBindFragDataLocation(GLState::program, 0, "fragColor");

  // set up vertex arrays
  GLState::vPosition = glGetAttribLocation( GLState::program, "vPosition" );
  GLState::vNormal = glGetAttribLocation( GLState::program, "vNormal" );

  // Retrieve transformation uniform variable locations
  GLState::ModelView = glGetUniformLocation( GLState::program, "ModelView" );
  GLState::NormalMatrix = glGetUniformLocation( GLState::program, "NormalMatrix" );
  GLState::ModelViewLight = glGetUniformLocation( GLState::program, "ModelViewLight" );
  GLState::Projection = glGetUniformLocation( GLState::program, "Projection" );

  GLState::objectVao.resize(sceneObjects.size());
  glGenVertexArrays( sceneObjects.size(), &GLState::objectVao[0] );

  GLState::objectBuffer.resize(sceneObjects.size());
  glGenBuffers( sceneObjects.size(), &GLState::objectBuffer[0] );

  for(unsigned int i=0; i < sceneObjects.size(); i++){
    glBindVertexArray( GLState::objectVao[i] );
    glBindBuffer( GL_ARRAY_BUFFER, GLState::objectBuffer[i] );
    size_t vertices_bytes = sceneObjects[i]->mesh.vertices.size()*sizeof(vec4);
    size_t normals_bytes  =sceneObjects[i]->mesh.normals.size()*sizeof(vec3);

    glBufferData( GL_ARRAY_BUFFER, vertices_bytes + normals_bytes, NULL, GL_STATIC_DRAW );
    size_t offset = 0;
    glBufferSubData( GL_ARRAY_BUFFER, offset, vertices_bytes, &sceneObjects[i]->mesh.vertices[0] );
    offset += vertices_bytes;
    glBufferSubData( GL_ARRAY_BUFFER, offset, normals_bytes,  &sceneObjects[i]->mesh.normals[0] );

    glEnableVertexAttribArray( GLState::vNormal );
    glEnableVertexAttribArray( GLState::vPosition );

    glVertexAttribPointer( GLState::vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );
    glVertexAttribPointer( GLState::vNormal, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(vertices_bytes));

  }



  glEnable( GL_DEPTH_TEST );
  glShadeModel(GL_SMOOTH);

  glClearColor( 0.8, 0.8, 1.0, 1.0 );

  //Quaternion trackball variables, you can ignore
  GLState::scaling  = 0;
  GLState::moving   = 0;
  GLState::panning  = 0;
  GLState::beginx   = 0;
  GLState::beginy   = 0;

  TrackBall::matident(GLState::curmat);
  TrackBall::trackball(GLState::curquat , 0.0f, 0.0f, 0.0f, 0.0f);
  TrackBall::trackball(GLState::lastquat, 0.0f, 0.0f, 0.0f, 0.0f);
  TrackBall::add_quats(GLState::lastquat, GLState::curquat, GLState::curquat);
  TrackBall::build_rotmatrix(GLState::curmat, GLState::curquat);

  GLState::scalefactor = 1.0;
  GLState::render_line = false;

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void drawObject(Object * object, GLuint vao, GLuint buffer){

  color4 material_ambient(object->shadingValues.color.x*object->shadingValues.Ka,
                          object->shadingValues.color.y*object->shadingValues.Ka,
                          object->shadingValues.color.z*object->shadingValues.Ka, 1.0 );
  color4 material_diffuse(object->shadingValues.color.x,
                          object->shadingValues.color.y,
                          object->shadingValues.color.z, 1.0 );
  color4 material_specular(object->shadingValues.Ks,
                           object->shadingValues.Ks,
                           object->shadingValues.Ks, 1.0 );
  float  material_shininess = object->shadingValues.Kn;

  color4 ambient_product  = GLState::light_ambient * material_ambient;
  color4 diffuse_product  = GLState::light_diffuse * material_diffuse;
  color4 specular_product = GLState::light_specular * material_specular;

  glUniform4fv( glGetUniformLocation(GLState::program, "AmbientProduct"), 1, ambient_product );
  glUniform4fv( glGetUniformLocation(GLState::program, "DiffuseProduct"), 1, diffuse_product );
  glUniform4fv( glGetUniformLocation(GLState::program, "SpecularProduct"), 1, specular_product );
  glUniform4fv( glGetUniformLocation(GLState::program, "LightPosition"), 1, lightPosition );
  glUniform1f(  glGetUniformLocation(GLState::program, "Shininess"), material_shininess );

  glBindVertexArray(vao);
  glBindBuffer( GL_ARRAY_BUFFER, buffer );
  glVertexAttribPointer( GLState::vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );
  glVertexAttribPointer( GLState::vNormal, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(object->mesh.vertices.size()*sizeof(vec4)) );

  mat4 objectModelView = GLState::sceneModelView*object->getModelView();


  glUniformMatrix4fv( GLState::ModelViewLight, 1, GL_TRUE, GLState::sceneModelView);
  glUniformMatrix3fv( GLState::NormalMatrix, 1, GL_TRUE, Normal(objectModelView));
  glUniformMatrix4fv( GLState::ModelView, 1, GL_TRUE, objectModelView);

  glDrawArrays( GL_TRIANGLES, 0, object->mesh.vertices.size() );

}


int main(void){

  GLFWwindow* window;

  glfwSetErrorCallback(error_callback);

  if (!glfwInit())
    exit(EXIT_FAILURE);

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  glfwWindowHint(GLFW_SAMPLES, 4);

  window = glfwCreateWindow(768, 768, "Raytracer", NULL, NULL);
  if (!window){
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwSetKeyCallback(window, keyCallback);
  glfwSetMouseButtonCallback(window, mouseClick);
  glfwSetCursorPosCallback(window, mouseMove);


  glfwMakeContextCurrent(window);
  gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
  glfwSwapInterval(1);

  switch(scene){
    case _SPHERE:
      initUnitSphere();
      break;
    case _SQUARE:
      initUnitSquare();
      break;
    case _BOX:
      initCornellBox();
      break;
  }

  initGL();

  while (!glfwWindowShouldClose(window)){

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    GLState::window_height = height;
    GLState::window_width  = width;

    glViewport(0, 0, width, height);


    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    mat4 track_ball =  mat4(GLState::curmat[0][0], GLState::curmat[1][0],
                            GLState::curmat[2][0], GLState::curmat[3][0],
                            GLState::curmat[0][1], GLState::curmat[1][1],
                            GLState::curmat[2][1], GLState::curmat[3][1],
                            GLState::curmat[0][2], GLState::curmat[1][2],
                            GLState::curmat[2][2], GLState::curmat[3][2],
                            GLState::curmat[0][3], GLState::curmat[1][3],
                            GLState::curmat[2][3], GLState::curmat[3][3]);

    GLState::sceneModelView  =  Translate(-cameraPosition) *   //Move Camera Back
    Translate(GLState::ortho_x, GLState::ortho_y, 0.0) *
    track_ball *                   //Rotate Camera
    Scale(GLState::scalefactor,
          GLState::scalefactor,
          GLState::scalefactor);   //User Scale

    GLfloat aspect = GLfloat(width)/height;

    switch(scene){
      case _SPHERE:
      case _SQUARE:
        GLState::projection = Perspective( 45.0, aspect, 0.01, 100.0 );
        break;
      case _BOX:
        GLState::projection = Perspective( 45.0, aspect, 4.5, 100.0 );
        break;
    }

    glUniformMatrix4fv( GLState::Projection, 1, GL_TRUE, GLState::projection);

    for(unsigned int i=0; i < sceneObjects.size(); i++){
      drawObject(sceneObjects[i], GLState::objectVao[i], GLState::objectBuffer[i]);
    }

    glfwSwapBuffers(window);
    glfwPollEvents();

  }

  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
