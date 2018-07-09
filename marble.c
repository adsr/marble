#define _GNU_SOURCE
#include <stdlib.h>
#include <time.h>
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
    vec3f_t orig;
    vec3f_t dir;
};

struct hit_s {
    vec3f_t point;
    float t;
    object_t *object;
    // TODO angle
    // TODO object/poly that was hit
};

struct poly_s {
    int n;
    vec3f_t p[4];
    vec3f_t min;
    vec3f_t max;
    object_t *object;
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
static int construct_bvh_axis_min_sort(const void *va, const void *vb, void *arg);
static void draw_objects();
static int ray_cast(ray_t *ray, bbox_t *bbox, hit_t *ret, int depth);
static int ray_intersects_bbox(ray_t *ray, bbox_t *bbox);
static int ray_intersects_rect(ray_t *ray, poly_t* poly, hit_t *ret);
static int ray_intersects_tri(ray_t *ray, poly_t *poly, hit_t *ret);
static void ray_init_shadow(ray_t *ray, hit_t *hit);
static void ray_init_ortho_from_cam_pixel(ray_t *ray, matrix44f_t *cam_to_world, int view_x, int view_y);
static void wait_for_quit();
static float vec3f_mag(vec3f_t *v);
static void vec3f_norm(vec3f_t *v);
static void vec3f_cross(vec3f_t *result, vec3f_t *v1, vec3f_t *v2);
static float vec3f_dot(vec3f_t *v1, vec3f_t *v2);
static void vec3f_mult_scalar(vec3f_t *v, float f);
static void vec3f_sub(vec3f_t *v, vec3f_t *o);
static void vec3f_add(vec3f_t *v, vec3f_t *o);
static void vec3f_set(vec3f_t *v, float x, float y, float z);
static void matrix44f_calc_transform(matrix44f_t *ret, vec3f_t *from, vec3f_t *to);
static void vec3f_mult_matrix_point(vec3f_t *v, matrix44f_t *m);
static void vec3f_mult_matrix_vector(vec3f_t *v, matrix44f_t *m);
static float randf(float min, float max);

scene_t scene;

int main(int argc, char **argv) {
    srand(time(NULL));
    memset(&scene, 0, sizeof(scene_t));
    scene.view_w = 320/4;
    scene.view_h = 480/4;
    scene.view_scale = 2;
    vec3f_set(&scene.dir_light, 1, 0, 0);
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
    int i, j, n;
    float size;
    size = 16.f;
    n = 3;
    object_t *object;
    // TODO read level from config file
    scene.objects_len = n*n;
    scene.objects = calloc(n*n, sizeof(object_t)); // TODO free
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            // TODO different object types
            object = &scene.objects[i*n+j];
            object->x = i * size;
            object->y = 0;
            object->z = j * size;
            object->wx = 8.f;
            object->hy = 32.f * (((i*n+j) + 1)/(1.f*n*n));
            object->lz = 8.f;
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
        poly.object = object;
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

    int i;
    for (i=0;i<poly->n;i++)
        printf("%f,%f,%f,\n", poly->p[i].v[0], poly->p[i].v[1], poly->p[i].v[2]);

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
    int rank[3];
    poly_t *poly;

    // Bail if 2 or less polys
    if (bbox->polys_len <= 2) {
        return;
    }

    // Find min, max of bbox
    bbox_find_min_max(bbox);

    // Rank split on each axis
    for (axis = 0; axis < 3; axis++) {
        rank[axis] = construct_bvh_axis(bbox, axis, 0);
    }

    // Pick best split
    if (       rank[0] != -1 && rank[0] <= rank[1] && rank[0] <= rank[2]) {
        axis = 0;
    } else if (rank[1] != -1 && rank[1] <= rank[0] && rank[1] <= rank[2]) {
        axis = 1;
    } else if (rank[2] != -1 && rank[2] <= rank[0] && rank[2] <= rank[1]) {
        axis = 2;
    } else {
        return; // Could not split
    }

    // Recurse
    construct_bvh_axis(bbox, axis, 1);
    construct_bvh(bbox->a);
    construct_bvh(bbox->b);
}

static int construct_bvh_axis(bbox_t *bbox, int axis, int alloc) {
    int ideal_i, d, s, i, j, clean_partition;

    // Bail if less than 2 ranges
    if (bbox->polys_len < 2) {
        return -1;
    }

    // Sort polys by min
    qsort_r(
        bbox->polys,
        bbox->polys_len,
        sizeof(poly_t),
        construct_bvh_axis_min_sort,
        (void*)&axis
    );

    // Find best split index
    ideal_i = bbox->polys_len / 2;
    for (d = 0; d <= ideal_i; d++) {
        for (s = -1; s <= 1; s += 2) {
            i = ideal_i + (d * s);
            if (i < 1 || i >= bbox->polys_len) {
                continue;
            }
            clean_partition = 1;
            for (j = 0; j < bbox->polys_len; j++) {
                if (i == j) {
                    continue;
                } else if (j < i && bbox->polys[j].max.v[axis] > bbox->polys[i].max.v[axis]) {
                    clean_partition = 0;
                    break;
                } else if (j > i && bbox->polys[j].min.v[axis] < bbox->polys[i].max.v[axis]) {
                    clean_partition = 0;
                    break;
                }
            }
            if (clean_partition) {
                goto construct_bvh_axis_found;
            }
            if (d == 0) {
                break;
            }
        }
    }

    // Not found
    return -1;

construct_bvh_axis_found:
    if (alloc) {
        bbox->a = calloc(1, sizeof(bbox_t)); // TODO free
        bbox->a->polys = bbox->polys;
        bbox->a->polys_len = i;
        bbox->b = calloc(1, sizeof(bbox_t)); // TODO free
        bbox->b->polys = bbox->polys + i;
        bbox->b->polys_len = bbox->polys_len - i;
    }

    return d;
}

static int construct_bvh_axis_min_sort(const void *va, const void *vb, void *arg) {
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
    ray_t ray;
    hit_t hit;
    int view_x, view_y;
    vec3f_t tmp;

    // Draw view pixel by pixel
    int shadow,light,bg;
    shadow = light = bg = 0;
    for (view_y = 0; view_y < scene.view_h; view_y++) {
        for (view_x = 0; view_x < scene.view_w; view_x++) {
            vec3f_set(
                &ray.orig,
                (1.f-(1.f*view_x / scene.view_w)) * (scene.bvh.max.v[0]) * -2.f,
                ((1.f*view_y / scene.view_h)) * (scene.bvh.max.v[0]) * 2.f + scene.bvh.max.v[1],
                ((1.f*view_x / scene.view_w)) * (scene.bvh.max.v[2]) * -2.f
            );
            vec3f_set(&ray.dir, 1.f, -1.f, 1.f);
//if (view_x %10 ==0 &&view_y %10==0) printf("%f,%f,%f,\n", ray.orig.v[0], ray.orig.v[1], ray.orig.v[2]);
//tmp = ray.dir;
//vec3f_mult_scalar(&tmp, 50.f);
//vec3f_add(&tmp, &ray.orig);
//if (view_x %10 ==0 &&view_y %10==0) printf("%f,%f,%f,\n", tmp.v[0], tmp.v[1], tmp.v[2]);


            if (ray_cast(&ray, &scene.bvh, &hit, 0)) {
                // hit
                ++light;

                SDL_SetRenderDrawColor(scene.renderer, 30+hit.object->x, 30+hit.object->z, 30+hit.object->hy, 255);
                SDL_RenderDrawPoint(scene.renderer, view_x, view_y);
            } else {
                // background
                ++bg;
                SDL_SetRenderDrawColor(scene.renderer, 0, 0, 0, 255);
                SDL_RenderDrawPoint(scene.renderer, view_x, view_y);
            }
        }
    }
    printf("shadow=%d light=%d bg=%d\n", shadow, light, bg);
    SDL_RenderPresent(scene.renderer);
}

static int ray_cast(ray_t *ray, bbox_t *bbox, hit_t *ret, int depth) {
    int i, n;
    hit_t tmp = {0};
    float tmin;
    if (!ray_intersects_bbox(ray, bbox)) return 0;
    if (bbox->a && ray_cast(ray, bbox->a, ret, depth+1)) return 1;
    if (bbox->b && ray_cast(ray, bbox->b, ret, depth+1)) return 1;
    n = 0;
    tmin = 0.f;
    for (i = 0; i < bbox->polys_len; i++) {
        // TODO different poly types
        if (ray_intersects_rect(ray, &bbox->polys[i], &tmp)) {
            n += 1;
            if (tmp.t == tmp.t && (tmin == 0.f || tmp.t < tmin)) {
                *ret = tmp;
                tmin = tmp.t;
            }
        }
    }
    return n;
}

static int ray_intersects_bbox(ray_t *ray, bbox_t *bbox) {
    vec3f_t reflect_dir;
    float t[6], tmin, tmax;

    vec3f_set(
        &reflect_dir,
        ray->dir.v[0] == 0.f ? FLT_MIN : 1.f / ray->dir.v[0],
        ray->dir.v[1] == 0.f ? FLT_MIN : 1.f / ray->dir.v[1],
        ray->dir.v[2] == 0.f ? FLT_MIN : 1.f / ray->dir.v[2]
    );

    t[0] = (bbox->min.v[0] - ray->orig.v[0]) * reflect_dir.v[0];
    t[1] = (bbox->max.v[0] - ray->orig.v[0]) * reflect_dir.v[0];
    t[2] = (bbox->min.v[1] - ray->orig.v[1]) * reflect_dir.v[1];
    t[3] = (bbox->max.v[1] - ray->orig.v[1]) * reflect_dir.v[1];
    t[4] = (bbox->min.v[2] - ray->orig.v[2]) * reflect_dir.v[2];
    t[5] = (bbox->max.v[2] - ray->orig.v[2]) * reflect_dir.v[2];

    tmin = fmax(fmax(fmin(t[0], t[1]), fmin(t[2], t[3])), fmin(t[4], t[5]));
    tmax = fmin(fmin(fmax(t[0], t[1]), fmax(t[2], t[3])), fmax(t[4], t[5]));

    if (tmax < 0) return 0;
    if (tmin > tmax) return 0;

    // t = tmin;
    return 1;
}

static int ray_intersects_rect(ray_t *ray, poly_t* poly, hit_t *ret) {
    poly_t tmp;

    tmp = *poly;
    tmp.n = 3;

    tmp.p[0] = poly->p[0];
    tmp.p[1] = poly->p[1];
    tmp.p[2] = poly->p[2];
    if (ray_intersects_tri(ray, &tmp, ret)) return 1;

    tmp.p[0] = poly->p[0];
    tmp.p[1] = poly->p[2];
    tmp.p[2] = poly->p[3];
    return ray_intersects_tri(ray, &tmp, ret);
}

static int ray_intersects_tri(ray_t *ray, poly_t *poly, hit_t *ret) {
    vec3f_t edge1, edge2, h, s, q;
    vec3f_t *vertex;
    float a, f, u, v, t;

    vertex = poly->p;

    printf("ray.orig(%3.1f %3.1f %3.1f) tri(%3.1f,%3.1f,%3.1f %3.1f,%3.1f,%3.1f %3.1f,%3.1f,%3.1f) ", 
        ray->orig.v[0], ray->orig.v[1], ray->orig.v[2],
        vertex[0].v[0], vertex[0].v[1], vertex[0].v[2],
        vertex[1].v[0], vertex[1].v[1], vertex[1].v[2],
        vertex[2].v[0], vertex[2].v[1], vertex[2].v[2]
    );

    edge1 = vertex[1];
    vec3f_sub(&edge1, &vertex[0]);
    edge2 = vertex[2];
    vec3f_sub(&edge2, &vertex[0]);

    vec3f_cross(&h, &ray->dir, &edge2);
    a = vec3f_dot(&edge1, &h);

    if (fabs(a) < 1e-8) {
        printf("parallel\n");
        return 0;
    }

    f = 1.f / a;

    s = ray->orig;
    vec3f_sub(&s, &vertex[0]);

    u = f * vec3f_dot(&s, &h);
    if (u < 0.f || u > 1.f) {
        printf("u<0 || u>1 (f=%f u=%f)\n", f, u);
        return 0;
    }

    vec3f_cross(&q, &s, &edge1);

    v = f * vec3f_dot(&ray->dir, &q);
    if (v < 0.0 || u + v > 1.0) {
        printf("v<0 || u+v>1\n");
        return 0;
    }

    t = f * vec3f_dot(&edge2, &q);
    if (t > 1e-8) {
        ret->point = ray->dir;
        vec3f_mult_scalar(&ret->point, t);
        vec3f_add(&ret->point, &ray->orig);
        ret->t = t;
        ret->object = poly->object;
        printf("YA\n");
        return 1;
    }

    printf("nah\n");
    return 0;
}

static void wait_for_quit() {
    SDL_Event event;
    while (1) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            break;
    }
}

static float vec3f_mag(vec3f_t *v) {
    return sqrtf(pow(v->v[0], 2.f) + pow(v->v[1], 2.f) + pow(v->v[2], 2.f));
}

static void vec3f_norm(vec3f_t *v) {
    float mag;
    mag = vec3f_mag(v);
    if (mag != 0.f) {
        v->v[0] /= mag;
        v->v[1] /= mag;
        v->v[2] /= mag;
    }
}

static void vec3f_cross(vec3f_t *result, vec3f_t *v1, vec3f_t *v2) {
    result->v[0] = (v1->v[0] * v2->v[2]) - (v1->v[2] * v2->v[1]);
    result->v[1] = (v1->v[2] * v2->v[0]) - (v1->v[0] * v2->v[2]);
    result->v[2] = (v1->v[0] * v2->v[1]) - (v1->v[1] * v2->v[0]);
}

static float vec3f_dot(vec3f_t *a, vec3f_t *b) {
    return (a->v[0] * b->v[0] + a->v[1] * b->v[1] + a->v[2] * b->v[2]);
}

static void vec3f_mult_scalar(vec3f_t *v, float f) {
    v->v[0] *= f;
    v->v[1] *= f;
    v->v[2] *= f;
}

static void vec3f_add(vec3f_t *v, vec3f_t *o) {
    v->v[0] += o->v[0];
    v->v[1] += o->v[1];
    v->v[2] += o->v[2];
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

static float randf(float min, float max) {
    float scale = rand() / (float) RAND_MAX;
    return min + scale * ( max - min );
}
