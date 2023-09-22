#include "common.h"
#include "vec3.h"
#include "ray.h"
#include "plane.h"
#include "sphere.h"
#include "camera.h"
#include "scene_objects.h"
#include "material.h"
#include "color.h"

#include <iostream>

using namespace std;

// {x0, x1, z0, z1, k}
double area_light[5] = { 1, 5, 1, 7, 5 };



color ray_color(const ray& r, const color& background, const hittable& world, int depth) {
    hit_record rec;
    //r.dir change 50/50 0,1,2
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
}

color path_ray_color(const ray& r, const color& background, const hittable& world, int depth) {
    hit_record object;
    // bounced enough times -> passed ray bounce limit (max)
    if (depth <= 0)
        return color(0, 0, 0);

    // if nothing was hit
    if (!world.hit(r, 0.001, infinity, object))
        return background;

    // Material material = ray.thingHit->material;
    // Color emittance = material.emittance;
    // equivalent to:
    color emitted = object.mat_ptr->emitted(object.u, object.v, object.p);

    //r.dir change 50/50 0,1,2

    ray scattered; // new ray
    scattered.origin = object.p;

    color attenuation; // albedo/light intensity


    // scatter calculates a random unit vector in the hemisphere of the normal where the obj was hit
    if (!object.mat_ptr->scatter(r, object, attenuation, scattered))
        return emitted;

    // phong light
    scattered.phongReflect(scattered, object.get_point(), area_light[4]);
    vec3 specular;

    //dot product should be between scattered.dir and object.normal
    double dotProduct = 0.0; // between scatter.dir and r.dir
    dotProduct += scattered.get_direction().vec[0] * object.get_normal().vec[0];
    dotProduct += scattered.get_direction().vec[1] * object.get_normal().vec[1];
    dotProduct += scattered.get_direction().vec[2] * object.get_normal().vec[2];

    double k_e = 5.0; // e of a shiny object
    double k_d = 0.6; // k_d of a shiny object

    double max_res = max(0.0, dotProduct);
    specular = attenuation * pow(max_res, k_e);

    int check = random_int(0, 1);
    //create shadow dir
    if (check == 0) {

        scattered.changeDirection(scattered, area_light[4]);
    }

    return emitted * k_d + specular * ray_color(scattered, background, world, depth - 1);
    //return emitted + attenuation * ray_color(scattered, background, world, depth - 1);

}

color dist_ray_color(const ray& r, const color& background, const hittable& world, int depth) {
    hit_record object;
    // bounced enough times -> passed ray bounce limit (max)
    if (depth <= 0)
        return color(0, 0, 0);

    // if nothing was hit
    if (!world.hit(r, 0.001, infinity, object))
        return background;

    // Material material = ray.thingHit->material;
    // Color emittance = material.emittance;
    // equivalent to:
    color emitted = object.mat_ptr->emitted(object.u, object.v, object.p);

    //r.dir change 50/50 0,1,2

    ray scattered; // new ray
    //scattered.changeDirection(scattered, area_light[4]);
    scattered.origin = object.p;

    color attenuation; // albedo

    // scatter calculates a random unit vector in the hemisphere of the normal where the obj was hit
    // return emitted (ambient light) because shadow/no scatter
    if (!object.mat_ptr->scatter(r, object, attenuation, scattered))
        return emitted;

    //* phong light - maybe put in check if statement
    //scattered.phongReflect(scattered, object.get_point(), area_light[4]);
    vec3 specular;

    //dot product should be between scattered.dir and object.normal
    double dotProduct = 0.0; // between scatter.dir and r.dir
    dotProduct += scattered.get_direction().vec[0] * object.get_normal().vec[0];
    dotProduct += scattered.get_direction().vec[1] * object.get_normal().vec[1];
    dotProduct += scattered.get_direction().vec[2] * object.get_normal().vec[2];

    double k_e = 5.0; // e of a shiny object
    double k_d = 0.6; // k_d of a shiny object

    double max_res = max(0.0, dotProduct);
    specular = attenuation * pow(max_res, k_e); //*/

    int check = random_int(0, 1);
    /*create shadow dir
    if (check == 0) {

        scattered.changeDirection(scattered, area_light[4]);
    } */

    return emitted * k_d + specular * ray_color(scattered, background, world, depth - 1);
    //return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
}

hittable_list scene() {
    hittable_list world;

    auto material_ground = make_shared<lambertian>(color(0.1, 0.1, 0.1));
    auto material_center = make_shared<lambertian>(color(0.24, 0.70, 0.54));
    auto material_right = make_shared<lambertian>(color(0.55, 0.57, 0.55));
    auto material_left = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    auto diffused_light = make_shared<diffuse_light>(color(4, 4, 4));


    world.add(make_shared<sphere>(point3(0.0, -100.5, 0.0), 100.0, material_ground));

    world.add(make_shared<sphere>(point3(-1.0, 1.0, -2.0), 1.5, material_right));
    world.add(make_shared<sphere>(point3(-4.0, 1.0, 1.0), 1.5, material_center));
    world.add(make_shared<sphere>(point3(0.0, 1.0, 4.0), 1.5, material_left));
    world.add(make_shared<xz_rect>(area_light[0], area_light[1], area_light[2], area_light[3], area_light[4], diffused_light));

    return world;
}

int main() {

    // Image

    auto aspect_ratio = 3.0 / 2.0;
    int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    int samples_per_pixel = 25; // 1000 - 500 - 100 - 50
    int max_depth = 5; // 50 - 50 - 25 - 15

    // World
    auto world = scene();

    // Camera
    point3 lookfrom(16, 4, 4); //(13, 2, 3)
    point3 lookat(0, 1, 0); //(0,1,-1)
    const vec3 vup(0, 1, 0);
    const auto dist_to_focus = 10.0;
    auto vfov = 30.0; // field of view angle
    auto aperture = 0.1;
    color background(0.0, 0.0, 0.0);

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                //*generate a random ray - path tracing
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);

                r.changeOrigin(r, lookfrom); // path tracing start from eye

                pixel_color += path_ray_color(r, background, world, max_depth); //*/


                /* DISTRIBUTED RAY TRACING
                for (int k = 0; k < 5; k++) {
                    auto u = (i + random_double()) / (image_width - 1); // k/i+1
                    auto v = (j + random_double()) / (image_height - 1); // l/j+1
                    ray r = cam.get_ray(u, v);
                    r.changeDirection(r, area_light[4]); // distributed ray tracing
                    pixel_color += dist_ray_color(r, background, world, max_depth);
                } //*/



            }
            //pixel_color = pixel_color / samples_per_pixel;
            print_pixel_colour(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}