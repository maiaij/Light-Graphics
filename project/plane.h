#include "common.h"
#include "bbox.h"
#include "hittable.h"

class xz_rect : public hittable {
public:
    xz_rect() {}

    xz_rect(
        double _x0, double _x1, double _z0, double _z1, double _k, shared_ptr<material> mat
    ) : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(box& output_box) const override {
        // The bounding box must have non-zero width in each dimension, so pad the Y
        // dimension a small amount.
        output_box = box(point3(x0, k - 0.0001, z0), point3(x1, k + 0.0001, z1));
        return true;
    }

public:
    shared_ptr<material> mp;
    double x0, x1, z0, z1, k;
};



bool xz_rect::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    auto t = (k - r.get_origin().y()) / r.get_direction().y();
    if (t < t_min || t > t_max)
        return false;

    auto x = r.get_origin().x() + t * r.get_direction().x();
    auto z = r.get_origin().z() + t * r.get_direction().z();
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;

    rec.u = (x - x0) / (x1 - x0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;
    auto outward_normal = vec3(0, 1, 0);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mp;

    return true;
}
