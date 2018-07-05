#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>

typedef struct vec3f_s vec3f_t;
typedef struct matrix44f_s matrix44f_t;
typedef struct ray_s ray_t;
typedef struct hit_s hit_t;
typedef struct poly_s poly_t;
typedef struct object_s object_t;
typedef struct bbox_s bbox_t;
typedef struct scene_s scene_t;

struct vec3f_s {
    float v[3];
};

struct matrix44f_s {
    float m[4][4];
};

struct ray_s {
    vec3f_t origin;
    vec3f_t direction;
};

struct hit_s {
    vec3f_t point;
    // TODO angle
    // TODO object/poly that was hit
};

struct poly_s {
    int n;
    vec3f_t p[4];
    vec3f_t min;
    vec3f_t max;
};

struct object_s {
    // TODO int type;
    float x;
    float y;
    float z;
    float hy;
    float wx;
    float lz;
    // TODO float hy[4];
    // TODO int fold;
    // TODO int r,g,b;
    // TODO other objects
};

struct bbox_s {
    vec3f_t min;
    vec3f_t max;
    poly_t *polys;
    int polys_len;
    int polys_cap;
    bbox_t *a;
    bbox_t *b;
};

struct scene_s {
    SDL_Renderer *renderer;
    SDL_Window *window;
    int view_w;
    int view_h;
    int view_scale;
    bbox_t bvh;
    vec3f_t dir_light;
    int objects_len;
    object_t *objects;
};

static void init_sdl();
static void deinit_sdl();
static void make_objects();
static void make_polys();
static void add_poly(poly_t *poly);
static void poly_find_min_max(poly_t *poly);
static void bbox_find_min_max(bbox_t* bbox);
static void construct_bvh(bbox_t *bbox);
static int construct_bvh_axis(bbox_t *bbox, int axis, int alloc);
static int construct_bvh_axis_sort(const void *va, const void *vb, void *arg);
static void draw_objects();
static int ray_cast(ray_t *ray, hit_t *ret);
static void ray_init_shadow(ray_t *ray, hit_t *hit);
static void ray_init_ortho_from_cam_pixel(ray_t *ray, matrix44f_t *cam_to_world, int view_x, int view_y);
static void wait_for_quit();
static float vec3f_mag(vec3f_t *v);
static void vec3f_norm(vec3f_t *v);
static void vec3f_cross(vec3f_t *result, vec3f_t *v1, vec3f_t *v2);
static void vec3f_mult_scalar(vec3f_t *v, float f);
static void vec3f_sub(vec3f_t *v, vec3f_t *o);
static void vec3f_set(vec3f_t *v, float x, float y, float z);
static void matrix44f_calc_transform(matrix44f_t *ret, vec3f_t *from, vec3f_t *to);
static void vec3f_mult_matrix_point(vec3f_t *v, matrix44f_t *m);
static void vec3f_mult_matrix_vector(vec3f_t *v, matrix44f_t *m);
static float randf(float min, float max);

scene_t scene;

int main(int argc, char **argv) {
    memset(&scene, 0, sizeof(scene_t));
    scene.view_w = 256;
    scene.view_h = 240;
    scene.view_scale = 2;
    init_sdl();
    make_objects();
    make_polys();
    construct_bvh(&scene.bvh);
    draw_objects();
    wait_for_quit();
    deinit_sdl();
    return 0;
}

static void init_sdl() {
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(
        scene.view_w * scene.view_scale,
        scene.view_h * scene.view_scale,
        0,
        &scene.window,
        &scene.renderer
    );
    SDL_RenderSetLogicalSize(scene.renderer, scene.view_w, scene.view_h);
    SDL_SetRenderDrawColor(scene.renderer, 0, 0, 0, 0);
    SDL_RenderClear(scene.renderer);
    //SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
}

static void deinit_sdl() {
    SDL_DestroyRenderer(scene.renderer);
    SDL_DestroyWindow(scene.window);
    SDL_Quit();
}

static void make_objects() {
    int i, j;
    float size;
    size = 16.0;
    object_t *object;
    // TODO read level from config file
    scene.objects_len = 100;
    scene.objects = calloc(100, sizeof(object_t)); // TODO free
    for (i = 0; i < 10; i++) {
        for (j = 0; j < 10; j++) {
            // TODO different object types
            object = &scene.objects[i*10+j];
            object->x = i * size;
            object->y = 0;
            object->z = j * size;
            object->wx = size;
            object->hy = randf(1.0, 10.0);
            object->lz = size;
        }
    }
}

static void make_polys() {
    int i;
    object_t *object;
    poly_t poly;
    vec3f_t v[8];

    memset(&poly, 0, sizeof(poly_t));
    poly.n = 4;
    for (i = 0; i < scene.objects_len; i++) {
        // TODO different object types
        object = &scene.objects[i];
        //    f-------g
        //   /|      /|
        //  / |     / |
        // b--|----c  |
        // |  e----|--h
        // | /     | /
        // a-------d
        //
        // origin: a
        // min:    a
        // max:    g
        // x-axis: ae
        // y-axis: ab
        // z-axis: ad
        vec3f_set(&v[0], object->x + 0,          object->y + 0,          object->z + 0         ); // a
        vec3f_set(&v[1], object->x + 0,          object->y + object->hy, object->z + 0         ); // b
        vec3f_set(&v[2], object->x + 0,          object->y + object->hy, object->z + object->lz); // c
        vec3f_set(&v[3], object->x + 0,          object->y + 0,          object->z + object->lz); // d
        vec3f_set(&v[4], object->x + object->wx, object->y + 0,          object->z + 0         ); // e
        vec3f_set(&v[5], object->x + object->wx, object->y + object->hy, object->z + 0         ); // f
        vec3f_set(&v[6], object->x + object->wx, object->y + object->hy, object->z + object->lz); // g
        vec3f_set(&v[7], object->x + object->wx, object->y + 0,          object->z + object->lz); // h

        // abcd
        poly.p[0] = v[0]; poly.p[1] = v[1]; poly.p[2] = v[2]; poly.p[3] = v[3]; add_poly(&poly);
        // efgh
        poly.p[0] = v[4]; poly.p[1] = v[5]; poly.p[2] = v[6]; poly.p[3] = v[7]; add_poly(&poly);
        // abfe
        poly.p[0] = v[0]; poly.p[1] = v[1]; poly.p[2] = v[5]; poly.p[3] = v[4]; add_poly(&poly);
        // dcgh
        poly.p[0] = v[3]; poly.p[1] = v[2]; poly.p[2] = v[6]; poly.p[3] = v[7]; add_poly(&poly);
        // aehd
        poly.p[0] = v[0]; poly.p[1] = v[4]; poly.p[2] = v[7]; poly.p[3] = v[3]; add_poly(&poly);
        // bfgc
        poly.p[0] = v[1]; poly.p[1] = v[5]; poly.p[2] = v[6]; poly.p[3] = v[2]; add_poly(&poly);
    }
}

static void add_poly(poly_t *poly) {
    // Ensure capacity
    if (scene.bvh.polys_len <= scene.bvh.polys_cap) {
        scene.bvh.polys_cap += 128;
        scene.bvh.polys = realloc(scene.bvh.polys, sizeof(poly_t) * scene.bvh.polys_cap); // TODO free
    }

    // Find min, max
    poly_find_min_max(poly);

    // Add poly
    scene.bvh.polys[scene.bvh.polys_len] = *poly;
    scene.bvh.polys_len += 1;
}

static void poly_find_min_max(poly_t *poly) {
    int axis, n;
    for (axis = 0; axis < 3; axis++) {
        poly->min.v[axis] = FLT_MAX;
        poly->max.v[axis] = FLT_MIN;
        for (n = 0; n < poly->n; n++) {
            if (poly->p[n].v[axis] < poly->min.v[axis]) poly->min.v[axis] = poly->p[n].v[axis];
            if (poly->p[n].v[axis] > poly->max.v[axis]) poly->max.v[axis] = poly->p[n].v[axis];
        }
    }
}

static void bbox_find_min_max(bbox_t* bbox) {
    int axis, i;
    poly_t *poly;
    for (int axis = 0; axis < 3; axis++) {
        bbox->min.v[axis] = FLT_MAX;
        bbox->max.v[axis] = FLT_MIN;
        for (i = 0; i < bbox->polys_len; i++) {
            poly = &bbox->polys[i];
            if (poly->min.v[axis] < bbox->min.v[axis]) bbox->min.v[axis] = poly->min.v[axis];
            if (poly->max.v[axis] > bbox->max.v[axis]) bbox->max.v[axis] = poly->max.v[axis];
        }
    }
}

static void construct_bvh(bbox_t *bbox) {
    int axis, i, n;
    int split[3];
    poly_t *poly;

    // Bail if 2 or less polys
    if (bbox->polys_len <= 2) return;

    // Find min, max of bbox
    bbox_find_min_max(bbox);

    // TODO construct_bvh_axis on each axis

    // TODO recurse construct_bvh(bbox->a/b)
}

static int construct_bvh_axis(bbox_t *bbox, int axis, int alloc) {
    int i, j, a_len, b_len, rank;
    if (bbox->polys_len < 2) {
        return -1;
    }
    qsort_r(
        bbox->polys,
        bbox->polys_len,
        sizeof(poly_t),
        construct_bvh_axis_sort,
        (void*)&axis
    );
    rank = -1;
    for (i = 0; i < bbox->polys_len; i++) {
        for (j = 0; j < bbox->polys_len; j++) {
        }
    }
    a_len = bbox->polys_len / 2;
    b_len = bbox->polys_len - a_len;
    if (alloc) {
        bbox->a = calloc(1, sizeof(bbox_t)); // TODO free
        bbox->a->polys = bbox->polys;
        bbox->a->polys_len = a_len;
        bbox->b = calloc(1, sizeof(bbox_t)); // TODO free
        bbox->b->polys = bbox->polys + a_len;
        bbox->b->polys_len = b_len;
    }

    rank = a_len - (bbox->polys_len / 2);
    if (rank < 0) {
        rank *= -1;
    }
    return rank;
}

static int construct_bvh_axis_sort(const void *va, const void *vb, void *arg) {
    poly_t *a, *b;
    int axis;
    a = (poly_t*)va;
    b = (poly_t*)vb;
    axis = *((int*)arg);
    if (a->min.v[axis] < b->min.v[axis]) {
        return -1;
    } else if (a->min.v[axis] > b->min.v[axis]) {
        return 1;
    }
    return 0;
}

static void draw_objects() {
    matrix44f_t cam_to_world;
    vec3f_t origin;
    vec3f_t cam;
    ray_t ray;
    hit_t hit;
    int view_x, view_y;

    // Set origin and cam coorinates
    vec3f_set(&origin, 0, 0, 0);
    vec3f_set(&cam, 0, 1, 9);

    // Calc transformation matrix from camera to world coorinate system
    matrix44f_calc_transform(&cam_to_world, &cam, &origin);

    // Draw view pixel by pixel
    for (view_x = 0; view_x < scene.view_w; view_x++) {
        for (view_y = 0; view_y < scene.view_h; view_y++) {
            ray_init_ortho_from_cam_pixel(&ray, &cam_to_world, view_x, view_y);
            if (ray_cast(&ray, &hit)) {
                ray_init_shadow(&ray, &hit);
                if (ray_cast(&ray, &hit)) {
                    // TODO draw in shadow
                } else {
                    // TODO draw in light
                }
            } else {
                // TODO draw bg color
            }
        }
    }
    SDL_RenderPresent(scene.renderer);
}

static int ray_cast(ray_t *ray, hit_t *ret) {
    // TODO
    return 0;
}

static void ray_init_shadow(ray_t *ray, hit_t *hit) {
    ray->direction = scene.dir_light;
    vec3f_mult_scalar(&ray->direction, -1.0);
    ray->origin = hit->point;
}

static void ray_init_ortho_from_cam_pixel(ray_t *ray, matrix44f_t *cam_to_world, int view_x, int view_y) {
    vec3f_set(&ray->direction, 0, 0, -1);
    vec3f_set(
        &ray->origin,
        scene.bvh.min.v[0] + ((scene.bvh.max.v[0] - scene.bvh.min.v[0]) * (view_x + 0.5)) / scene.view_w,
        scene.bvh.min.v[1] + ((scene.bvh.max.v[1] - scene.bvh.min.v[1]) * (view_y + 0.5)) / scene.view_h,
        0
    );
    vec3f_mult_matrix_point(&ray->origin, cam_to_world);
    vec3f_mult_matrix_vector(&ray->direction, cam_to_world);
}

static void wait_for_quit() {
    SDL_Event event;
    while (1) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            break;
    }
}

static float vec3f_mag(vec3f_t *v) {
    return sqrtf(pow(v->v[0], 2.0) + pow(v->v[1], 2.0) + pow(v->v[2], 2.0));
}

static void vec3f_norm(vec3f_t *v) {
    float mag;
    mag = vec3f_mag(v);
    if (mag == 0.0) {
        v->v[0] /= mag;
        v->v[1] /= mag;
        v->v[2] /= mag;
    }
}

static void vec3f_cross(vec3f_t *result, vec3f_t *v1, vec3f_t *v2) {
    result->v[0] = (v1->v[1] * v2->v[2]) - (v1->v[2] * v2->v[1]);
    result->v[1] = (v1->v[0] * v2->v[2]) - (v1->v[2] * v2->v[0]);
    result->v[2] = (v1->v[0] * v2->v[1]) - (v1->v[1] * v2->v[0]);
}

static void vec3f_mult_scalar(vec3f_t *v, float f) {
    v->v[0] *= f;
    v->v[1] *= f;
    v->v[2] *= f;
}

static void vec3f_sub(vec3f_t *v, vec3f_t *o) {
    v->v[0] -= o->v[0];
    v->v[1] -= o->v[1];
    v->v[2] -= o->v[2];
}

static void vec3f_set(vec3f_t *v, float x, float y, float z) {
    v->v[0] = x;
    v->v[1] = y;
    v->v[2] = z;
}

static void matrix44f_calc_transform(matrix44f_t *ret, vec3f_t *from, vec3f_t *to) {
    vec3f_t up, x, y, z;

    z = *from;
    vec3f_sub(&z, to);
    vec3f_norm(&z);

    vec3f_set(&up, 0, 1, 0);
    vec3f_cross(&x, &up, &z);

    vec3f_cross(&y, &z, &x);

    ret->m[0][0] = x.v[0];
    ret->m[0][1] = x.v[1];
    ret->m[0][2] = x.v[2];

    ret->m[1][0] = y.v[0];
    ret->m[1][1] = y.v[1];
    ret->m[1][2] = y.v[2];

    ret->m[2][0] = z.v[0];
    ret->m[2][1] = z.v[1];
    ret->m[2][2] = z.v[2];

    ret->m[0][3] = 0;
    ret->m[1][3] = 0;
    ret->m[2][3] = 0;
    ret->m[3][3] = 1;
}

static void vec3f_mult_matrix_point(vec3f_t *v, matrix44f_t *m) {
    float a, b, c, w;
    a = v->v[0] * m->m[0][0] + v->v[1] * m->m[1][0] + v->v[2] * m->m[2][0] + m->m[3][0];
    b = v->v[0] * m->m[0][1] + v->v[1] * m->m[1][1] + v->v[2] * m->m[2][1] + m->m[3][1];
    c = v->v[0] * m->m[0][2] + v->v[1] * m->m[1][2] + v->v[2] * m->m[2][2] + m->m[3][2];
    w = v->v[0] * m->m[0][3] + v->v[1] * m->m[1][3] + v->v[2] * m->m[2][3] + m->m[3][3];
    v->v[0] = a / w;
    v->v[1] = b / w;
    v->v[2] = c / w;
}

static void vec3f_mult_matrix_vector(vec3f_t *v, matrix44f_t *m) {
    float a, b, c;
    a = v->v[0] * m->m[0][0] + v->v[1] * m->m[1][0] + v->v[2] * m->m[2][0];
    b = v->v[0] * m->m[0][1] + v->v[1] * m->m[1][1] + v->v[2] * m->m[2][1];
    c = v->v[0] * m->m[0][2] + v->v[1] * m->m[1][2] + v->v[2] * m->m[2][2];
    v->v[0] = a;
    v->v[1] = b;
    v->v[2] = c;
}

static float randf(float min, float max) {
    float scale = rand() / (float) RAND_MAX;
    return min + scale * ( max - min );
}

//                SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
//                SDL_RenderDrawPoint(renderer,
