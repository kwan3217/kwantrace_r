//! Find the intersection between a ray and a sphere

use std::fs::File;
use std::io::Write;
// Linear algebra stuff from
// https://rust-lang-nursery.github.io/rust-cookbook/science/mathematics/linear_algebra.html
use ndarray::{array,Array,Array1,ArrayView1};
// Vectors are 1D and functions like dot know how to treat them as rows or columns -- IE
// you do M@v (where v is a compatible column vector) as M.dot(v),
//    and v@M (where v is a compatible row    vector) as v.dot(M)
//    and v.u (where both vectors are just vectors, don't care row or column) as v.dot(u)

/// Intersect a ray and the unit sphere
///
pub fn ray_sphere(r0:ArrayView1<f64>,v:ArrayView1<f64>)->Option<f64> {
    /* The sphere is defined by x^2+y^2+z^2=1, and the
     ray is defined by x=x0+vx*t, y=y0+vy*t, z=z0+zy*t

     (x0+vx*t)^2+(y0+vy*t)^2+(z0+vz*t)^2-1=0
     x0^2+2*x0*vx*t+vx^2*t^2 +
     y0^2+2*y0*vy*t+vy^2*t^2 +
     z0^2+2*z0*vz*t+vz^2*t^2 - 1 =0
     t^2(vx^2+vy^2+vz^2)+
     t  (2*x0*vx+2*y0*vy+2*z0*vz)+
        (x0^2+y0^2+z0^2-1)=0
     a=vx^2+vy^2+vz^2=v.v
     b=2*(x0*vx+y0*vy+z0*vz)=2*(r0.v)
     c=(x0^2+y0^2+z0^2-1)=r0.r0-1
    */
    let a=v.dot(&v);
    let b=2.0*r0.dot(&v);
    let c=r0.dot(&r0)-1.0;
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


fn main()->std::io::Result<()> {
    const n_rows:usize=120;
    const n_cols:usize=160;
    let mut img=Array::<u8,_>::zeros((n_rows,n_cols));
    let v =array![ 0.0, 0.0, 1.0];
    for i_row in 0..n_rows {
        for i_col in 0..n_cols {
            let r0=array![ ((i_col as f64)/(n_cols as f64)-0.5)*4.0, ((i_row as f64)/(n_rows as f64)-0.5)*3.0,-2.0];
            match ray_sphere(r0.view(),v.view()) {
                Some(t) => {
                    println!("{i_row},{i_col} Intersect: t={t}");
                    img[[i_row,i_col]]=(t*128.0) as u8;
                },
                None    => println!("No intersect")
            }
        }
    }
    let mut file=File::create("out.ppm")?;
    write!(file,"P5 {n_cols} {n_rows} 255\n")?;
    for i_row in 0..n_rows {
        for i_col in 0..n_cols {
            file.write(&[img[[i_row, i_col]]])?;
        }
    }
    Ok(())
}
