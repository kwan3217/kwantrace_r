use crate::vector::{Direction, HMatrix, Matrix3x3, Ray, Vector};
use crate::transform::{Transform, TransformList, Translate};

pub trait Render {
    fn intersect_local(&self, rv_b:&Ray) ->Option<f64>;
    fn M_rb(&self)->&HMatrix;
    fn M_br(&self)->&HMatrix;
    fn translate(&mut self,x:f64,y:f64,z:f64) {
        let T=Translate::make(x,y,z);
        self.get_transforms().push(Box::new(T));
    }
    fn intersect(&self, rv_r:&Ray)->Option<f64> {
        self.intersect_local(&self.M_br().mul_ray(rv_r))
    }
    fn get_transforms(&mut self)->&mut TransformList;
    fn set_M_rb(&mut self,M_rb:HMatrix);
    fn set_M_br(&mut self,M_br:HMatrix);
    fn prepare_render(&mut self) {
        let M_rb=self.get_transforms().get_matrix();
        let M_br=M_rb.inv();
        self.set_M_rb(M_rb);
        self.set_M_br(M_br);
    }
}

pub struct Sphere {
    _m_rb:HMatrix,
    _m_br:HMatrix,
    _transforms:TransformList,
}

impl Sphere {
    pub fn make()->Sphere {
        Sphere{_m_rb:HMatrix::identity(), _m_br: HMatrix::identity(), _transforms: vec![] }
    }
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
    fn intersect_local(&self, rv: &Ray) -> Option<f64> {
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

    fn M_rb(&self) -> &HMatrix {
        &self._m_rb
    }

    fn M_br(&self) -> &HMatrix {
        &self._m_br
    }

    fn set_M_rb(&mut self,M_rb:HMatrix) {
        self._m_rb=M_rb;
    }

    fn set_M_br(&mut self,M_br:HMatrix) {
        self._m_br=M_br;
    }

    fn get_transforms(&mut self) -> &mut TransformList {
        &mut self._transforms
    }
}

pub struct Union {
    pub itemList:Vec<Box<dyn Render>>,
    _transforms:TransformList,
    _m_rb:HMatrix,
    _m_br:HMatrix,
}

impl Union {
    pub fn make() -> Union {
        Union{_m_rb:HMatrix::identity(), _m_br: HMatrix::identity(), _transforms: vec![],itemList:vec![] }
    }
}

impl Render for Union {
    fn prepare_render(&mut self) {
        let M_rb=self.get_transforms().get_matrix();
        let M_br=M_rb.inv();
        self.set_M_rb(M_rb);
        self.set_M_br(M_br);
        for this_render in &mut self.itemList {
            this_render.prepare_render();
        }
    }

    fn intersect_local(&self, rv_b: &Ray) -> Option<f64> {
        /* In many languages, we would keep track of the closest
           valid parameter, by tracking it. The initial value
            is either very large or literally infinity, so that
            any valid parameter is less than it.

            Here instead we use the Option enum. We check
            the first thing and keep it as a Some(t) or None. We
            then iterate through the other things. If this one
            is better (has an intersection while we don't yet,
            or intersection is closer than current best) we keep
            the best intersection. When we are done, we return
            the best intersection, without having to check for
            thinks like are we still pointing at infinity.

            Note that we are in intersect_local() for the union,
            but are calling intersect() for the children. Each
            child performs its own transform and therefore what
            it considers to be M_rb is actually M_ib where i is
            the intermediate frame (the Union body frame). This requires
            one matrix-ray transform for the Union, and one for each
            child. Maybe later we will concatenate the
            reference-from-intermediate and intermediate-from-body
            transformations in prepare_render(). */
        let mut result=self.itemList[0].intersect(rv_b);
        for this_render in &self.itemList[1..] {
            let this_result= this_render.intersect(rv_b);
            match this_result {
                Some(_) => {
                    if result.is_none() {
                        result=this_result;
                    } else {
                        if this_result<result {
                            result=this_result;
                        }
                    }
                },
                None    => ()
            };

        }
        result
    }

    fn M_rb(&self) -> &HMatrix {
        &self._m_rb
    }

    fn M_br(&self) -> &HMatrix {
        &self._m_br
    }

    fn set_M_rb(&mut self,M_rb:HMatrix) {
        self._m_rb=M_rb;
    }

    fn set_M_br(&mut self,M_br:HMatrix) {
        self._m_br=M_br;
    }

    fn get_transforms(&mut self) -> &mut TransformList {
        &mut self._transforms
    }

}