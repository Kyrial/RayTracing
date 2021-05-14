//////////////////////////////////////////////////////////////////////////////
//
//  --- Object.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////


#include "common.h"


vec3 vec4ToVec3(vec4 v) {
    return vec3(v.x, v.y, v.z);
}

vec4 vec3ToVec4(vec3 v) {
    return vec4(v.x, v.y, v.z, 0);
}



/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Sphere::intersect(vec4 p0, vec4 V){
  IntersectionValues result;
  //TODO: Ray-sphere setup
  //garder l'objet identifier
 

  double t=raySphereIntersection(p0, V);
  
  result.t = t; //distance
  
  // pt d'intersection = origine du rayon + distance * vect unitaire
  result.P = p0 + t * V; 
  

  //vect normal a la sphere = pt d'intersection - centre de la sphere
  result.N = result.P - vec3ToVec4(center);
  //result.N.w = 1;
  //result.P.w = 1;

  result.name = this->name;
  


  return result;
}




//p0 c'est l'origine du rayon, v la direction, 
/* -------------------------------------------------------------------------- */
/* ------ Ray = p0 + t*V  sphere at origin center and radius radius    : Find t ------- */
double Sphere::raySphereIntersection(vec4 p0, vec4 V){
  double t   = std::numeric_limits< double >::infinity();
 
  //TODO: Ray-sphere intersection;
//vec3 x ;
//length(center);
//t = (-b + sqrt(delta)) / (2 * length(V) * a);
// center - pow(radius, 2);
 //2d.(o-c) - 4 * dd * ||o - c||^2 - r^2 
  //vec4 origine = vec4(0.0, 0.0, 0.0, 1.0);

  vec3 V3 = normalize(vec4ToVec3(V)); //V direction du rayon 


  //vec3 distcentre  = vec4ToVec3(p0 - origine);
  vec3 distcentre = vec4ToVec3(p0 - center);

  double a = dot(V3, V3);
  double b = 2 * dot(distcentre, V3);
  double c = dot(distcentre, distcentre)- (radius*radius); 
  double delta = (b * b) - 4 * a * c;
  //std::cout << " delta = " << delta << std::endl;
  if (delta > EPSILON) {
     // t = (-b - sqrt(delta)) / (2.0 * length(V) * a);
      double t1= (-b - sqrt(delta)) / (2.0 * a);
      double t2 = (-b + sqrt(delta)) / (2.0 * a);
      if (t1 > 0 && t2 > EPSILON)
          t = min(t1, t2);
      if (t1 <= EPSILON)
          t = t2;
      if (t2 <= EPSILON)
          t = t1;
      if (t1 <= EPSILON && t2 <= EPSILON)
          t= std::numeric_limits< double >::infinity();
      //std::cout << "||p-c||² - r² " << length((p0 + t * V) - center) * length((p0 + t * V) - center) - radius * radius<< std::endl << std::endl;;

  }
 
  

  return t;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
Object::IntersectionValues Square::intersect(vec4 p0, vec4 V){
  IntersectionValues result;
  //TODO: Ray-square setup
  double t = raySquareIntersection(p0, V);
  result.t = t;
  result.P = p0 + t * V;
  result.N = vec3ToVec4(normal);

  result.name = this->name;
  return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
double Square::raySquareIntersection(vec4 p0, vec4 V){
  double t   = std::numeric_limits< double >::infinity();
  //TODO: Ray-square intersection;
    vec3 V3 = vec4ToVec3(V); //V direction du rayon 
    vec3 P0 = vec4ToVec3(p0); //origine du rayon
    vec3 point3 = vec4ToVec3(point);
    //      D- o . n
    //t = ------------
    //      d . n
   

    // vec3 nominateur = cross( D , on);

  //  double D = dot(point3, normal);
  //  double on = dot(P0, normal);

    //double nominateur = D - on;
    double nominateur = dot(point3, normal) - dot(P0, normal);
    double denominateur = dot(V3,normal);
    
    if (denominateur != 0) {
        double delta = nominateur / denominateur;
        // std::cout << "  dot(V3,normal) = " << dot(V3, normal) << std::endl;
        // std::cout << " nom et denom " << nominateur << " et " << denominateur << std::endl;
       //  std::cout << " delta = " << delta << std::endl;
      //   std::cout << " P0 + (V3 * delta) = " << P0 + (V3 * delta) << std::endl;
        if (delta > EPSILON) {

            vec3 pointIntersec = P0 + (V3 * delta);
            
            vec3 upRight = vec4ToVec3(mesh.vertices[1]);
            vec3 basRight = vec4ToVec3(mesh.vertices[2]);
            vec3 upLeft = vec4ToVec3(mesh.vertices[5]);
            vec3 basLeft = vec4ToVec3(mesh.vertices[3]);

     //       std::cout << " point1 " << dot(normal, cross(upRight - basRight, upRight - pointIntersec)) << std::endl;
   //         std::cout << " point2 " << dot(normal, cross(basRight - basLeft, basRight - pointIntersec)) << std::endl;// << std::endl;
    //        std::cout << " point3 " << dot(normal, cross(basLeft -upLeft, basLeft - pointIntersec)) << std::endl;
    //        std::cout << " point4 " << dot(normal, cross(upLeft - upRight, upLeft- pointIntersec)) << std::endl<<std::endl;
           /* if(dot(normal, cross(upRight - basRight, upRight - pointIntersec)) <EPSILON
                &&
                dot(normal, cross(basRight - basLeft, basRight - pointIntersec)) < EPSILON
                &&
                dot(normal, cross(basLeft - upLeft, basLeft - pointIntersec)) < EPSILON
                &&
                dot(normal, cross(upLeft - upRight, upLeft - pointIntersec)) < EPSILON                
                )*/

        //    std::cout << "(p-a).n=0 " << dot((pointIntersec - point), normal) << std::endl << std::endl;;

        //    std::cout << "br - ur " << (basRight - upRight) << std::endl;
        //    std::cout << "P - ur " << (pointIntersec - upRight) << std::endl;
         //   std::cout << "cross1 " << cross(basRight - upRight, pointIntersec - upRight) <<std::endl<< std::endl;

           double a = dot(normal, cross(basRight - upRight, pointIntersec - upRight));
           double b = dot(normal, cross(basLeft - basRight, pointIntersec - basRight));

           double c = dot(normal, cross(upLeft - basLeft, pointIntersec - basLeft));
           double d = dot(normal, cross(upRight - upLeft, pointIntersec - upLeft));

           if( (a < EPSILON && b < EPSILON && c < EPSILON&& d < EPSILON) ||
               (a > EPSILON && b > EPSILON && c > EPSILON && d > EPSILON)
               )

              // if( a* b >EPSILON)


                return delta;         
        }



//        v0v1 ^ v0P ->dot(n,) >0
    }
     
    
  return t;
}

//vec3 test = upRight - basRight;
//std::cout << " x, y " << test.x << "  " << test.y << "   " << std::endl;
//vec3 test2 = upRight - pointIntersec;
//std::cout << " x, y " << test2.x << "  " << test2.y << "   " << std::endl;
