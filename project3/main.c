#define GL_GLEXT_PROTOTYPES 1
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GL/glx.h>
#include <GL/glext.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define RENDER_BOX

#define UNPADDED_FACE_SIZE 13
#define SCALE_FACTOR 5

#define NUM_PROPERTIES 7

#define N 100
#define INTENSITY 1.0

float eye[] = {0.0, 0.0, 10.0};
float viewpt[] = {0.0, 0.0, 0.0};
float up[] = {0.0, 1.0, 0.0};

typedef struct vector_t {
    float x, y, z, nx, ny, nz;
} vector;

typedef struct face_t {
    int v1, v2, v3;
} face;

typedef struct model_t {
    GLfloat *data;
    int numVertex;
} model;

typedef struct ray_t {
  float x, y, z;
} ray;

float inputX = 0;
float inputY = 0;

float angle;
float timeStart = 0;

model *modelData = NULL;
unsigned int modelBuf = 1;
unsigned int shaderProgID;

float dotR(ray *a, ray *b) {
    return a->x * b->x +
        a->y * b->y +
        a->z * b->z;
}

ray *crossR(ray *a, ray *b) {
    ray *ret = (ray *) malloc(sizeof(ray));
    ret->x = (a->y * b->z) - (a->z * b->y);
    ret->y = (a->z * b->x) - (a->x * b->z);
    ret->z = (a->x * b->y) - (a->y * b->x);
    return ret;
}

void vect_divide(vector *vect, double scalar) {
    vect->x /= scalar;
    vect->y /= scalar;
    vect->z /= scalar;
}

void ray_divide(ray *r, double scalar) {
    r->x /= scalar;
    r->y /= scalar;
    r->z /= scalar;
}

ray *differenceR(ray *a, ray *b) {
    ray *ret = (ray *) malloc(sizeof(ray));
    ret->x = a->x - b->x;
    ret->y = a->y - b->y;
    ret->z = a->z - b->z;
    return ret;
}

double length(vector *v) {
    return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
}

double lengthR(ray *v) {
    return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
}

void normalize(vector *v) {
    vect_divide(v, length(v));
}

// moller Trumbore intersections
bool intersection(ray v0, ray v1, ray v2, ray *origin, ray *dir, ray *outVar) {

  ray *edge1 = (ray *)malloc(sizeof(ray));
  ray *edge2 = (ray *)malloc(sizeof(ray));
  ray *tvec = (ray *)malloc(sizeof(ray));
  float det, inv_det, u, v, t;
  
  // find edges sharing vertex 0
  edge1 = differenceR(&v1, &v0);
  edge2 = differenceR(&v2, &v0);
  
  // calculate determinant - det = r0 *(r1 x r2)
  ray *pvec = (ray *)malloc(sizeof(ray)); 
  pvec = crossR(dir, edge2);
  
  det = dotR(edge1, pvec);
  
  inv_det = 1.0 / det;
  
  // calculate distance from vertex0 to ray origin
  tvec = differenceR(origin, &v0);
  
  // calculate u - and test
  u = dotR(tvec, pvec) * inv_det; 
  if (u < 0.0 || u > 1.0) // intersection misses triangle
    return false;
  
  // calculate v - and test
  ray *qvec = (ray *)malloc(sizeof(ray)); 
  qvec = crossR(tvec, edge1);
  
  v = dotR(dir, qvec) * inv_det;
  if (v < 0.0 || u + v > 1.0) // intersection misses triangle
    return false;
  
  // calculate t - ray intersection with triangle 
  t = dotR(edge2, qvec) * inv_det;
  
  outVar->x = t; // return vec t,u,v
  outVar->y = u;
  outVar->z = v;
  
  free(edge1);
  free(edge2);
  free(tvec);
  free(pvec);
  free(qvec);
  
  return true;
}

model *readPlyFile(char const *filename) {
    model *ret = (model *)malloc(sizeof(model));
    float sumX = 0.0, sumY = 0.0, sumZ = 0.0, maxX = 0, maxY = 0, maxZ = 0;

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Invalid file %s\n", filename);
        exit(-1);
    }
    // Skip the 'ply\n'
    fseek(fp, 4, SEEK_SET);

    int endHeader = 0;
    int properties = 0;
    int numVertex = 0, numFaces = 0;
    while (!feof(fp) && !endHeader) {
        char *line = NULL;
        size_t len = 0;
        getline(&line, &len, fp);
        switch(line[0]) {
            case 'p':
                properties += 1;
                break;
            case 'e':
                if(line[1] == 'n') {
                    endHeader = 1;
                }
                else {
                    char *type = (char *)malloc(7);
                    int size = 0;
                    sscanf(line, "element %s %d", type, &size);
                    if (type[0] == 'v') {
                        numVertex = size;
                    }
                    else {
                        numFaces = size;
                    }
                }
        }
    }
    if (properties != NUM_PROPERTIES) {
        fprintf(stderr, "Error, ply file does not have exactly %d properties, it has %d properties\n", NUM_PROPERTIES, properties);
        exit(-1);
    }

    vector *vertices = (vector *)calloc(numVertex, sizeof(vector));
	  for(int i = 0; i < numVertex; i++) {
		  //read in a line
		  char one[10], two[10], three[10], four[10], five[10], six[10];
		  fscanf(fp, "%s %s %s %s %s %s", one, two, three, four, five, six);

		  //convert to float
		  //put into vertices
		  vertices[i].x = strtof(one,NULL);
		  vertices[i].y = strtof(two,NULL);
		  vertices[i].z = strtof(three,NULL);
		  vertices[i].nx = strtof(four,NULL);
		  vertices[i].ny = strtof(five,NULL);
		  vertices[i].nz = strtof(six,NULL);
	  }


    face *faces = (face *)calloc(numFaces, sizeof(face));
	  for(int i = 0; i < numFaces; i++)
	  {
		  //read in a line
		  char one[10], two[10], three[10], four[10];
		  fscanf(fp, "%s %s %s %s",one, two, three, four);
		
		  //convert to int
		  //put into faces
		  faces[i].v1 = atoi(two);
		  faces[i].v2 = atoi(three);
		  faces[i].v3 = atoi(four);
	  }    
    // centering and scaling
    for (int i = 0; i < numVertex; i++) {

      sumX += vertices[i].x;
      sumY += vertices[i].y;
      sumZ += vertices[i].z;
      if (vertices[i].x > maxX) maxX = vertices[i].x;
      if (vertices[i].y > maxY) maxY = vertices[i].y;
      if (vertices[i].z > maxZ) maxZ = vertices[i].z;
    }
    sumX /= numVertex;
    sumY /= numVertex;
    sumZ /= numVertex;
	
    float scale = maxX;
    if (maxY > scale) scale = maxY;
    if (maxZ > scale) scale = maxZ;
    scale /= SCALE_FACTOR;

    for (int i = 0; i < numVertex; i++) {
        vertices[i].x -= sumX;
        vertices[i].y -= sumY;
        vertices[i].z -= sumZ;
        vect_divide(&vertices[i], scale);
    }
    // end centering and scaling
    GLfloat *facesAndNormals = (GLfloat *) calloc(numFaces * 6, sizeof(GLfloat) * 3);   
    int ndx = 0;
    int ndx2 = numFaces * 9;
    for (int i = 0; i < numFaces; i++) {

	      vector *v1 = &vertices[faces[i].v1];
        vector *v2 = &vertices[faces[i].v2];
        vector *v3 = &vertices[faces[i].v3];

        facesAndNormals[ndx++] = v1->x;
        facesAndNormals[ndx++] = v1->y;
        facesAndNormals[ndx++] = v1->z;
        facesAndNormals[ndx++] = v2->x;
        facesAndNormals[ndx++] = v2->y;
        facesAndNormals[ndx++] = v2->z;
        facesAndNormals[ndx++] = v3->x;
        facesAndNormals[ndx++] = v3->y;
        facesAndNormals[ndx++] = v3->z;
        facesAndNormals[ndx2++] = v1->nx;
        facesAndNormals[ndx2++] = v1->ny;
        facesAndNormals[ndx2++] = v1->nz;
        facesAndNormals[ndx2++] = v2->nx;
        facesAndNormals[ndx2++] = v2->ny;
        facesAndNormals[ndx2++] = v2->nz;
        facesAndNormals[ndx2++] = v3->nx;
        facesAndNormals[ndx2++] = v3->ny;
        facesAndNormals[ndx2++] = v3->nz;

    }

    ret->data = facesAndNormals;
    ret->numVertex = numFaces * 3;

    free(vertices);
    free(faces);

    return ret;
}

char *read_shader_program(char const *filename) {
    FILE *fp;
    char *content = NULL;
    int fd, count;
    fd = open(filename,O_RDONLY);
    if (fd == 0) {
        fprintf(stderr, "Shader source `%s` not found\n", filename);
        exit(-1);
    }
    count = lseek(fd,0,SEEK_END);
    close(fd);
    content = (char *)calloc(1,(count+1));
    fp = fopen(filename,"r");
    count = fread(content,sizeof(char),count,fp);
    content[count] = '\0';
    fclose(fp);
    return content;
}

void view_volume() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0,1.0,1.0,27.5);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye[0],eye[1],eye[2],viewpt[0],viewpt[1],viewpt[2],up[0],up[1],up[2]);
}

void set_uniform(float *diffuse, float *specular, float shininess, int programID) {
  
  int location;
  location = glGetUniformLocation(programID,"diffuse");
  glUniform4fv(location, 1, diffuse); // 4 elems - floats - type vector : for sending color  
  location = glGetUniformLocation(programID,"specular");  
  glUniform4fv(location, 1, specular); // 4 elems - floats - type vector : for sending color
  location = glGetUniformLocation(programID,"shininess");
  glUniform1f(location, shininess); // 1 elem
}

void set_color(GLfloat *destColor, GLfloat r, GLfloat g, GLfloat b, GLfloat intensity) {

  destColor[0] = r;
  destColor[1] = g;
  destColor[2] = b;
  destColor[3] = intensity;
}

void set_vertex(ray *destVert, GLfloat x, GLfloat y, GLfloat z) {

  destVert->x = x;
  destVert->y = y;
  destVert->z = z;
}

void draw_scene(float defColor[], GLfloat colorData[], GLfloat wallData[], GLfloat wallNorm[]) {      
  glClear(GL_DEPTH_BUFFER_BIT);
      // draw the bunny
  glPushMatrix();
      set_uniform(defColor, defColor, 5.0, shaderProgID);  
	  glDrawArrays(GL_TRIANGLES, 0, modelData->numVertex);
  glPopMatrix();

  for (int jj = 0; jj < 10; jj++){
    
    float color[3] = { colorData[(jj * 3) + 0], colorData[(jj * 3) + 1], colorData[(jj * 3) + 2]};
      if((jj * 3) >= 12 && (jj * 3) <= 24) { 
        set_uniform(color, color, 1.0, shaderProgID);
      } else {
        set_uniform(color, color, 8.0, shaderProgID);
      }
      
    glBegin(GL_TRIANGLES);
        glNormal3f(wallNorm[(jj * 9) + 0], wallNorm[(jj * 9) + 1], wallNorm[(jj * 9) + 2]);
        glVertex3f(wallData[(jj * 9) + 0], wallData[(jj * 9) + 1], wallData[(jj * 9) + 2]);
        glNormal3f(wallNorm[(jj * 9) + 3], wallNorm[(jj * 9) + 4], wallNorm[(jj * 9) + 5]);
        glVertex3f(wallData[(jj * 9) + 3], wallData[(jj * 9) + 4], wallData[(jj * 9) + 5]);
        glNormal3f(wallNorm[(jj * 9) + 6], wallNorm[(jj * 9) + 7], wallNorm[(jj * 9) + 8]);
        glVertex3f(wallData[(jj * 9) + 6], wallData[(jj * 9) + 7], wallData[(jj * 9) + 8]);
    glEnd();

    glBegin(GL_QUADS); // floor
      glVertex3f(-4.0, -1.6, -4.0);
      glVertex3f(4.0, -1.6, -4.0);
      glVertex3f(4.0, -1.6, 7.0);
      glVertex3f(-4.0, -1.6, 7.0); 
    glEnd();

      GLfloat lightColor[] = {1.0, 1.0, 1.0};
      set_uniform(lightColor, lightColor, 8.0, shaderProgID);
    glBegin(GL_QUADS); // area light
      glVertex3f(-1.5, 3.98, -1.5);
      glVertex3f(1.5, 3.98, -1.5);
      glVertex3f(1.5, 3.98, 1.5);
      glVertex3f(-1.5, 3.98, 1.5); 
    glEnd();

      set_uniform(defColor, defColor, 8.0, shaderProgID);
    glBegin(GL_QUADS); // ceiling
      glVertex3f(-4.0, 3.98, -4.0);
      glVertex3f(4.0, 3.98, -4.0);
      glVertex3f(4.0, 3.98, 6.0);
      glVertex3f(-4.0, 3.98, 6.0); 
    glEnd();

  } 
}

void renderScene() {

  glClearColor(0.35,0.35,0.35,1.0);
  glClearAccum(0.0,0.0,0.0,0.0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT); 

	glLoadIdentity();
	gluLookAt(eye[0],eye[1],eye[2],viewpt[0],viewpt[1],viewpt[2],up[0],up[1],up[2]);

  #ifdef RENDER_BOX
                   
  GLfloat colorData[] = {0.99, 0.9, 0.79, 0.99, 0.9, 0.79,  //backWall
                         0.99, 0.9, 0.79, 0.99, 0.9, 0.79,  //floor
                         1.0, 0.0, 1.0, 1.0, 0.0, 1.0,      //rightWall
                         0.0, 1.0, 0.0, 0.0, 1.0, 0.0,      //leftWall
                         1.0, 1.0, 1.0, 1.0, 1.0, 1.0};     //roof  
                
  GLfloat wallData[] = {-4, -1.7, -4, -4,  4,   -4,  4,  4,  -4,
		                    -4, -1.7, -4,  4, -1.7, -4,  4,  4,  -4,  //backWall
		                    -4, -1.7, -4, -4, -1.7,  7,  4, -1.7, 7,
		                    -4, -1.7, -4,  4, -1.7, -4,  4, -1.7, 7,  //floor
		                     4,  4,   -4,  4, -1.7, -4,  4, -1.7, 6,
		                     4,  4,   -4,  4,  4,    6,  4, -1.7, 6,  //rightWall
		                    -4,  4,   -4, -4, -1.7, -4, -4, -1.7, 6,
		                    -4,  4,   -4, -4,  4,    6, -4, -1.7, 6,  //leftWall
		                    -4,  4,   -4, -4,  4,    6,  4,  4,   6,
		                    -4,  4,   -4,  4,  4,   -4,  4,  4,   6}; //roof

  GLfloat wallNorm[] = { 0,  0, 1,  0,  0, 1,  0,  0, 1,
	                       0,  0, 1,  0,  0, 1,  0,  0, 1,  //backWall
	                       0,  1, 0,  0,  1, 0,  0,  1, 0,
	                       0,  1, 0,  0,  1, 0,  0,  1, 0,  //floor
                        -1,  0, 0, -1,  0, 0, -1,  0, 0,
	                      -1,  0, 0, -1,  0, 0, -1,  0, 0,  //rightWall
	                       1,  0, 0,  1,  0, 0,  1,  0, 0,
	                       1,  0, 0,  1,  0, 0,  1,  0, 0,  //leftWall
	                       0, -1, 0,  0, -1, 0,  0, -1, 0,
	                       0, -1, 0,  0, -1, 0,  0, -1, 0}; //roof

  float bunnyColor[] = {1.0, .75, .6, 1.0};
  // allocating big structure to hold all our data
  // in mem: [ wall data | vertices of faces]
  GLfloat *totalData = (GLfloat *) calloc((modelData->numVertex) + 30, sizeof(GLfloat) * 3); 
  memcpy(totalData, wallData, 90 * sizeof(GLfloat));
  memcpy(totalData + 90, modelData->data, (modelData->numVertex) * 3 * sizeof(GLfloat)); 

  // allocating big structure to hold all our normals
  // in mem: [ wall data | vertices of faces]
  GLfloat *totalNorm = (GLfloat *) calloc((modelData->numVertex) + 30, sizeof(GLfloat) * 3); 
  memcpy(totalNorm, wallNorm, 90 * sizeof(GLfloat));
  memcpy(totalNorm + 90, modelData->data + ((modelData->numVertex) * 3), (modelData->numVertex) * 3 * sizeof(GLfloat)); 
  #endif

  // --------------------- instant radiosity ------------------------------
  
  ray *randRay = (ray *)calloc(1, sizeof(ray)); // dir vector
  ray *hitPt = (ray *)calloc(1, sizeof(ray));
  ray *origin = (ray*)calloc(1, sizeof(ray));
  ray *dir = (ray*)calloc(1, sizeof(ray));
  ray *finalRay = (ray *)calloc(1, sizeof(ray));
  ray *intersectRay = (ray *) calloc(1, sizeof(ray));
  ray *reflectionRay = (ray *) calloc(1, sizeof(ray)); 
  
  srand(1000);
  // pick a random ray
  for (int i = 0; i < N; i++) {

    // pick random ray from light - area light 1 x 1
    float startPtX = (((float)rand()/(float)(RAND_MAX)) * 2.0) - 1; 
    float startPtZ = (((float)rand()/(float)(RAND_MAX)) * 2.0) - 1; 

    // add gl_light pointing in dir of ray
    float rayNew_pos[] = { startPtX, 4.0, startPtZ };
    float rayNew_color[] = {1.0, 1.0, 1.0, INTENSITY};

      // create hemisphere [0, 2pi] off of light
    double azimuth = (((double)rand()/(double)(RAND_MAX)) * (2*M_PI)); // random [0, 2pi]
    double elevation = asin(((double)rand()/(double)(RAND_MAX) * 1)); // arcsin(r)
    
    randRay->x = (-sin(azimuth) * cos(elevation));
    randRay->y = -1 * (sin(elevation));
    randRay->z = (cos(azimuth) * cos(elevation));
   
    float rayNew_dirn[] = {randRay->x, randRay->y, randRay->z, 1.0}; 

    glLightfv(GL_LIGHT0,GL_POSITION, rayNew_pos);
    glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION, rayNew_dirn);
    glLightfv(GL_LIGHT0,GL_DIFFUSE, rayNew_color);
    glLightfv(GL_LIGHT0,GL_SPECULAR, rayNew_color);
    glLightf(GL_LIGHT0,GL_CONSTANT_ATTENUATION, 1.0);
    glLightf(GL_LIGHT0,GL_LINEAR_ATTENUATION, .1);
    glLightf(GL_LIGHT0,GL_QUADRATIC_ATTENUATION, .01);
    
    finalRay->x = rayNew_pos[0] + randRay->x;
    finalRay->y = rayNew_pos[1] + randRay->y;
    finalRay->z = rayNew_pos[2] + randRay->z;     
    
    origin->x = startPtX;
    origin->y = 4.0;
    origin->z = startPtZ; 
    
    ray norm;
      
    float minimumHit = INFINITY;
    // add loop to go over ALL tri's of wall (two tris) and all tris of bunny
    for (int k = 0; k < ((modelData->numVertex*3)+90); k+=9) { 
      // calculate intersections                                          
      dir->x = randRay->x;
      dir->y = randRay->y;
      dir->z = randRay->z;
         
      // grabing triangle verticies for checking intersections
      ray v0;
        set_vertex(&v0, totalData[k], totalData[k+1], totalData[k+2]);
      ray v1;
        set_vertex(&v1, totalData[k+3], totalData[k+4], totalData[k+5]);
      ray v2;
        set_vertex(&v2, totalData[k+6], totalData[k+7], totalData[k+8]);
      ray v0n;
        set_vertex(&v0n, totalNorm[k], totalNorm[k+1], totalNorm[k+2]);;
      ray v1n;
        set_vertex(&v1n, totalNorm[k+3], totalNorm[k+4], totalNorm[k+5]);
      ray v2n;
        set_vertex(&v2n, totalNorm[k+6], totalNorm[k+7], totalNorm[k+8]);
      
      norm.x = (v0n.x + v1n.x + v2n.x)/3;
      norm.y = (v0n.y + v1n.y + v2n.y)/3;
      norm.z = (v0n.z + v1n.z + v2n.z)/3;
       
       // - check t - ray v0, v1, v2 are triangle coordinates (verts)
       if( intersection(v0, v1, v2, origin, dir, intersectRay) ) { // hit  
            // hit store lowest t value - indicating our first intersection
          if(minimumHit > intersectRay->x && intersectRay->x > 1) {
            minimumHit = intersectRay->x;
          }
       }  // end of if                  
      } // end of inner for loop
      
      if(minimumHit < INFINITY) {
              
        hitPt->x = rayNew_pos[0]*(1-minimumHit) + (finalRay->x *minimumHit);
        hitPt->y = rayNew_pos[1]*(1-minimumHit) + (finalRay->y *minimumHit);
        hitPt->z = rayNew_pos[2]*(1-minimumHit) + (finalRay->z *minimumHit);               
                            
        reflectionRay->x = finalRay->x - (2 * dotR(finalRay, &norm) * norm.x );
        reflectionRay->y = finalRay->y - (2 * dotR(finalRay, &norm) * norm.y );
        reflectionRay->z = finalRay->z - (2 * dotR(finalRay, &norm) * norm.z );

        // attach opGL spotlight in dir of reflected ray   
        float lightnew_pos[] = { hitPt->x, hitPt->y, hitPt->z };
        float lightnew_dir[] = { reflectionRay->x, reflectionRay->y, reflectionRay->z, 1.0 };
        glLightfv(GL_LIGHT0,GL_POSITION, lightnew_pos);
        glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION, lightnew_dir);                                   
        glLightf(GL_LIGHT0,GL_CONSTANT_ATTENUATION, 1.0);
        glLightf(GL_LIGHT0,GL_LINEAR_ATTENUATION, .1);
        glLightf(GL_LIGHT0,GL_QUADRATIC_ATTENUATION, .01);
             
        GLfloat newColor[4];       
                  
        // test x, y, z constant values with our insterect ray to find wall color
        if(hitPt->x <= (-4 + .0005) && hitPt->x >= (-4 - .0005)) {   
          set_color(newColor, colorData[18], colorData[19], colorData[20], INTENSITY);   
        }           
        else if(hitPt->x <= (4 + .0005) && hitPt->x >= (4 - .0005)) {
          set_color(newColor, colorData[12], colorData[13], colorData[14], INTENSITY);
        }           
        else if(hitPt->z <= (-4 + .0005) && hitPt->z >= (-4 - .0005)) {
          set_color(newColor, colorData[0], colorData[1], colorData[2], INTENSITY);
        }            
        else if(hitPt->y <= (-1.7 + .0005) && hitPt->y >= (-1.7 - .0005)) { 
          set_color(newColor, colorData[6], colorData[7], colorData[8], INTENSITY);
        }           
        else {
          set_color(newColor, bunnyColor[0], bunnyColor[1], bunnyColor[2], INTENSITY);
        }
        
        glLightfv(GL_LIGHT0,GL_DIFFUSE, newColor);
        glLightfv(GL_LIGHT0,GL_SPECULAR, newColor);

        rayNew_color[0] = newColor[0];
        rayNew_color[1] = newColor[1];
        rayNew_color[2] = newColor[2];
         
        draw_scene(bunnyColor, colorData, wallData, wallNorm);
        glAccum(GL_ACCUM, 1.0/(N/4));
      
      } // end of if 
   //} // end of outer for loop - do two bounces  
  }
  
  glAccum(GL_RETURN,1.0);
  glutSwapBuffers();
      
  free(hitPt);
  free(reflectionRay);
  free(intersectRay);
  free(dir);
  free(randRay);
  free(finalRay);
  free(origin);
  free(totalData);
  free(totalNorm);
  sleep(2);
}

int compileShader(GLuint shaderID, char const *shaderName) {
    int result;
    glCompileShader(shaderID);
    glGetShaderiv(shaderID, GL_COMPILE_STATUS, &result);
    if(result == GL_FALSE) {
        GLint logSize = 0;
        glGetShaderiv(shaderID, GL_INFO_LOG_LENGTH, &logSize);
        char *compileLog = (char *)malloc(logSize);
        GLsizei written = 0;
        glGetShaderInfoLog(shaderID, logSize, &written, compileLog);
        fprintf(stderr, "Error compiling %s shader:\n%s\n", shaderName, compileLog);
        glDeleteShader(shaderID);
        return GL_FALSE;
    }
    return GL_TRUE;
}

unsigned int set_shaders() {
    char *vs, *fs;
    GLuint v, f, p;

    v = glCreateShader(GL_VERTEX_SHADER);
    f = glCreateShader(GL_FRAGMENT_SHADER);
    vs = read_shader_program("main.vert");
    fs = read_shader_program("main.frag");
    glShaderSource(v,1,(const char **)&vs,NULL);
    glShaderSource(f,1,(const char **)&fs,NULL);
    free(vs);
    free(fs);
    if (compileShader(v, "vertex") == GL_FALSE) {
        exit(-1);
    }
    if (compileShader(f, "fragment") == GL_FALSE) {
        exit(-1);
    }
    p = glCreateProgram();
    glAttachShader(p,f);
    glAttachShader(p,v);
    glLinkProgram(p);
    glUseProgram(p);
    return p;
}

void HandleKeyDown(unsigned char key, int x, int y) {
    switch(key)
	{
		case 'w':
			inputY = -0.25;
			break;
        case 'W':
            inputY = -1.0;
            break;
		case 's':
			inputY = 0.25;
			break;
        case 'S':
            inputY = 1.0;
            break;
		case 'a':
			inputX = -0.25;
			break;
        case 'A':
            inputX = -1.0;
            break;
		case 'd':
			inputX = 0.25;
			break;
        case 'D':
            inputX = 1.0;
            break;
        case 'p':
            printf("Current camera position: %lf, %lf, %lf\n", eye[0], eye[1], eye[2]);
            printf("Current view point: %lf, %lf, %lf\n", viewpt[0], viewpt[1], viewpt[2]);
            break;
		case 'q':
            glDeleteBuffers(1, &modelBuf);
			exit(1);
	}
}

void HandleKeyUp(unsigned char key, int x, int y) {
    switch(key)
	{
		case 'w':
        case 'W':
			inputY = (inputY < 0) ? 0 : inputY;
			break;
		case 's':
        case 'S':
			inputY = (inputY > 0) ? 0 : inputY;
			break;
		case 'a':
        case 'A':
			inputX = (inputX < 0) ? 0 : inputX;
			break;
		case 'd':
        case 'D':
			inputX = (inputX > 0) ? 0 : inputX;
			break;
		case 'q':
            glDeleteBuffers(1, &modelBuf);
			exit(1);
	}
}

void send_model(model *model) {
    glBindBuffer(GL_ARRAY_BUFFER, modelBuf); 
    glBufferData(GL_ARRAY_BUFFER, sizeof(vector) * model->numVertex, model->data, GL_STATIC_DRAW); 
    glVertexPointer(3, GL_FLOAT, 3 * sizeof(GLfloat), 0);
    glNormalPointer(GL_FLOAT, 3 * sizeof(GLfloat), (void *)(3 * model->numVertex * sizeof(GLfloat)));
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
}

int main(int argc, char **argv) {
    modelData = readPlyFile("bunny.ply");
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_RGBA | GLUT_ACCUM |GLUT_DOUBLE);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(768, 768);
    glutCreateWindow("Bunny");
    
    glEnable(GL_DEPTH_TEST);
    
    send_model(modelData);
    view_volume();
    shaderProgID = set_shaders();
    
    glutDisplayFunc(renderScene);
	  glutIdleFunc(renderScene); // comment out to see one ray
    glutKeyboardFunc(HandleKeyDown);
	glutKeyboardUpFunc(HandleKeyUp);
    glutMainLoop();
   
	return 0;
}
