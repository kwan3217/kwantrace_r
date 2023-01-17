use crate::vector::{Ray, Vector};

pub(crate) trait Render {
    fn intersectLocal(&self, rv:&Ray)->Option<f64>;
}

pub struct Sphere {

}

impl Render for Sphere {
    /// Intersect a ray and the unit sphere
    ///
    /// The sphere is defined by x^2+y^2+z^2=1, and the
    /// ray is defined by x=x0+vx*t, y=y0+vy*t, z=z0+zy*t
    ///
    ///  * (x0+vx*t)^2+(y0+vy*t)^2+(z0+vz*t)^2-1=0
    ///  * x0^2+2*x0*vx*t+vx^2*t^2 + y0^2+2*y0*vy*t+vy^2*t^2 + z0^2+2*z0*vz*t+vz^2*t^2 - 1 =0
    ///  * t^2(vx^2+vy^2+vz^2)+ t  (2*x0*vx+2*y0*vy+2*z0*vz)+    (x0^2+y0^2+z0^2-1)=0
    ///  * a=vx^2+vy^2+vz^2=v.v
    ///  * b=2*(x0*vx+y0*vy+z0*vz)=2*(r0.v)
    ///  * c=(x0^2+y0^2+z0^2-1)=r0.r0-1
    fn intersectLocal(&self, rv: &Ray) -> Option<f64> {
        let a=rv.v.dot(&rv.v);
        let b=2.0*rv.r0.dot(&rv.v);
        let c=rv.r0.dot(&rv.r0)-1.0;
        let d=b*b-4.0*a*c;
        if d<0.0 {
            None
        } else {
            let tp=(-b+d.sqrt())/(2.0*a);
            let tm=(-b-d.sqrt())/(2.0*a);
            if tp<0.0 && tm<0.0 {
                // Both intersections behind camera
                None
            } else if tp<0.0 {
                // tp is behind camera, return tm
                Some(tm)
            } else if tm<0.0 {
                // tm is behind camera, return tp
                Some(tp)
            } else if tp<tm {
                // both are in front, tp is closer
                Some(tp)
            } else {
                // both are in front, tm is closer
                Some(tm)
            }
        }
    }
}