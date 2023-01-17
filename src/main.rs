#![warn(missing_docs)]
//! Implementation of a ray-tracer in Rust

use std::fs::File;
use std::io::Write;
use vector::{Direction, HMatrix, Matrix3x3, Position, Vector};

mod vector;

mod transform;

struct Ray {
    r0:Position,
    v:Direction,
}

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
fn ray_sphere(rv:&Ray)->Option<f64> {
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


fn main()->std::io::Result<()> {
    const n_cols:usize=1920;
    const n_rows:usize=1080;
    let mut img=[[0;n_cols];n_rows];
    let mut rv=Ray{r0:Position{x:0.0,y:0.0,z:0.0},v:Direction{x:0.0,y:0.0,z:1.0}};
    for (i_row,row) in img.iter_mut().enumerate() {
        for (i_col,pix) in row.iter_mut().enumerate() {
            rv.r0.x=((i_col as f64)/(n_cols as f64)-0.5)*16.0/4.0;
            rv.r0.y=((i_row as f64)/(n_rows as f64)-0.5)* 9.0/4.0;
            rv.r0.z=-2.0;
            match ray_sphere(&rv) {
                Some(t) => {
                    *pix=(t * 128.0) as u8;
                },
                None    => ()
            }
        }
        print!(".");
    }
    let mut file=File::create("out.ppm")?;
    write!(file,"P5 {n_cols} {n_rows} 255\n")?;
    for row in img.iter() {
        file.write(row)?;
    }
    Ok(())
}
