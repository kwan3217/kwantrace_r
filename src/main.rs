#![warn(missing_docs)]
//! Implementation of a ray-tracer in Rust

use std::fs::File;
use std::io::Write;
use vector::{Direction, HMatrix, Matrix3x3, Position, Ray, Vector};
use crate::render::{Render,Sphere};

mod vector;
mod transform;
mod render;

fn main()->std::io::Result<()> {
    const n_cols:usize=1920;
    const n_rows:usize=1080;
    let mut img=[[0;n_cols];n_rows];
    let mut rv=Ray{r0:Position{x:0.0,y:0.0,z:0.0},v:Direction{x:0.0,y:0.0,z:1.0}};
    let sphere=Sphere{};
    for (i_row,row) in img.iter_mut().enumerate() {
        for (i_col,pix) in row.iter_mut().enumerate() {
            rv.r0.x=((i_col as f64)/(n_cols as f64)-0.5)*16.0/4.0;
            rv.r0.y=((i_row as f64)/(n_rows as f64)-0.5)* 9.0/4.0;
            rv.r0.z=-2.0;
            match sphere.intersectLocal(&rv) {
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
