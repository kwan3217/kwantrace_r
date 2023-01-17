#![warn(missing_docs)]
//! Implementation of a ray-tracer in Rust

use std::fs::File;
use std::io::Write;

/// Position vector, a 3-element vector which participates in translation.
///
///When thinking about the vector as a homogeneous vector, all
///vectors of this type are considered to have an implicit w=1 component.
///This is handled in the matrix multiplication of a position vector.
struct Position {
    x:f64,
    y:f64,
    z:f64,
}

/// Direction vector, a 3-element vector which does not participate in translation.
///
///When thinking about the vector as a homogeneous vector, all
///vectors of this type are considered to have an implicit w=1 component.
///This is handled in the matrix multiplication of a position vector.
struct Direction {
    x:f64,
    y:f64,
    z:f64,
}

///Operations which are common to both Position and Direction, and operations
///between them.
trait Vector {
    /// Accessor for x component of a vector, since traits can't require fields,
    /// and therefore we can't require
    fn _x(&self)->f64;
    fn _y(&self)->f64;
    fn _z(&self)->f64;
    /// Dot product between two vectors. This ignores the implied w component
    /// of either vector
    fn dot(&self,other:&impl Vector)->f64 {
        self._x()*other._x()+
        self._y()*other._y()+
        self._z()*other._z()
    }
}

impl Vector for Position {
    fn _x(&self)->f64 {self.x}
    fn _y(&self)->f64 {self.y}
    fn _z(&self)->f64 {self.z}
}

impl Vector for Direction {
    fn _x(&self)->f64 {self.x}
    fn _y(&self)->f64 {self.y}
    fn _z(&self)->f64 {self.z}
}

impl Position {
    fn add_dir(&self,rhs:&Direction) -> Position {
        Position{x:self.x+rhs.x,
            y:self.y+rhs.y,
            z:self.z+rhs.z}
    }
}

impl Direction {
    fn add_dir(&self,rhs:&Direction) -> Direction {
        Direction{x:self.x+rhs.x,
                  y:self.y+rhs.y,
                  z:self.z+rhs.z}
    }
}


/// Represents a 3x3 matrix
///
/// The operations that we need:
/// * (/)  Mp - Matrix-Position multiplication
/// * (/)  Md - Matrix-Direction multiplication
/// * (/) -Md - Negated Matrix-Direction multiplication
/// * (/)  M^-1 - Matrix inverse
/// * ( )  M^-T - matrix transpose inverse=(M^-1)^T. This is used in transforming normal vectors,
///       which are direction vectors. As such, we only do this operation on the upper 3x3 matrix
/// The inverse of a homogeneous matrix is also homogeneous.
struct Matrix3x3 {
    e:[[f64;3];3],
}

impl Matrix3x3 {
    fn mul(&self,rhs:&Matrix3x3)->Matrix3x3 {
        // looped implementation
        /*
        let mut e: [[f64; 3]; 3] = [[0.0; 3]; 3];
        for m in 0..3 {
            for p in 0..3 {
                e[m][p] = self.e[m][0] * rhs.e[0][p];
                for n in 1..3 {
                    e[m][p] += self.e[m][n] * rhs.e[n][p];
                }
            }
        Matrix3x3 { e: e }

         */
        // manually unrolled implementation, I expect it to be right on the edge of readability.
        // Each cell is the dot product between the corresponding row from left and col from right.
        // In each cell, the inner two indices are the same for each term and increment in the cell,
        //   the first a index is the row index for the cell and will be the same in all terms in a row,
        //   the second b index is the column index for the cell and will be the same for all terms in a column
        let a=&self.e; //shorten the names so we can fit on the page better
        let b=&rhs.e;
        Matrix3x3{e:[[a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0],a[0][0]*b[0][1]+a[0][1]*b[1][1]+a[0][2]*b[2][1],a[0][0]*b[0][2]+a[0][1]*b[1][2]+a[0][2]*b[2][2]],
                     [a[1][0]*b[0][0]+a[1][1]*b[1][0]+a[1][2]*b[2][0],a[1][0]*b[0][1]+a[1][1]*b[1][1]+a[1][2]*b[2][1],a[1][0]*b[0][2]+a[1][1]*b[1][2]+a[1][2]*b[2][2]],
                     [a[2][0]*b[0][0]+a[2][1]*b[1][0]+a[2][2]*b[2][0],a[2][0]*b[0][1]+a[2][1]*b[1][1]+a[2][2]*b[2][1],a[2][0]*b[0][2]+a[2][1]*b[1][2]+a[2][2]*b[2][2]]]}
    }
    fn mul_pos(&self, rhs: &Position) -> Position {
        Position{x:self.e[0][0]*rhs._x()+self.e[0][1]*rhs._y()+self.e[0][2]*rhs._z(),
                 y:self.e[1][0]*rhs._x()+self.e[1][1]*rhs._y()+self.e[1][2]*rhs._z(),
                 z:self.e[2][0]*rhs._x()+self.e[2][1]*rhs._y()+self.e[2][2]*rhs._z()}
    }
    fn mul_dir(&self, rhs: &Direction) -> Direction {
        Direction{x:self.e[0][0]*rhs._x()+self.e[0][1]*rhs._y()+self.e[0][2]*rhs._z(),
                  y:self.e[1][0]*rhs._x()+self.e[1][1]*rhs._y()+self.e[1][2]*rhs._z(),
                  z:self.e[2][0]*rhs._x()+self.e[2][1]*rhs._y()+self.e[2][2]*rhs._z()}
    }
    fn mul_neg_dir(&self, rhs: &Direction) -> Direction {
        Direction{x:-(self.e[0][0]*rhs._x()+self.e[0][1]*rhs._y()+self.e[0][2]*rhs._z()),
                  y:-(self.e[1][0]*rhs._x()+self.e[1][1]*rhs._y()+self.e[1][2]*rhs._z()),
                  z:-(self.e[2][0]*rhs._x()+self.e[2][1]*rhs._y()+self.e[2][2]*rhs._z())}
    }
    /// Explicit formula for 3x3 inverse from
    /// <https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices>.
    /// This effectively does Cramer's rule and claims to be "efficient", which I interpret
    /// as that any improvement will be a different, substantially more difficult algorithm.
    ///
    /// It has 18 multiplies and 9 subtracts to calculate
    /// the capital letters, three multiplies, two adds, and a divide for the deteminant inverse,
    /// and 9 more multiplies by the determinant inverse, for a total of 30 multiplies, 11 adds/subtracts,
    /// and 1 divide.
    fn inv(&self) -> Matrix3x3 {
        let a=self.e[0][0]; let b=self.e[0][1]; let c=self.e[0][2];
        let d=self.e[1][0]; let e=self.e[1][1]; let f=self.e[1][2];
        let g=self.e[2][0]; let h=self.e[2][1]; let i=self.e[2][2];
        //Note that the order of these capital letter elements is transposed. Note that terms for
        //D, B, H, and F are reversed from original reference in order to get rid of the extra minus sign
        let A=e*i-f*h;let D=c*h-b*i;let G=b*f-c*e;
        let B=f*g-d*i;let E=a*i-c*g;let H=c*d-a*f;
        let C=d*h-e*g;let F=b*g-a*h;let I=a*e-b*d;
        let detinv=1.0/(a*A+b*B+c*C);
        Matrix3x3{e:[[A*detinv,D*detinv,G*detinv],
                     [B*detinv,E*detinv,H*detinv],
                     [C*detinv,F*detinv,I*detinv]]}
    }
}

struct HMatrix {
    M:Matrix3x3,
    T:Direction,
}

impl HMatrix {
    fn mul(&self, rhs: &HMatrix) -> HMatrix {
        HMatrix {M:self.M.mul(&self.M),
                 T:self.M.mul_dir(&rhs.T).add_dir(&self.T)}
    }
    fn mul_dir(&self, rhs: &Direction) -> Direction {
        self.M.mul_dir(rhs)
    }
    fn mul_pos(&self, rhs: &Position) -> Position {
        self.M.mul_pos(rhs).add_dir(&self.T)
    }
    fn inv(&self) -> HMatrix {
        let Am1=self.M.inv();
        HMatrix{T:Am1.mul_neg_dir(&self.T),M:Am1}
    }
}

trait Transform {
    fn get_matrix(&self)->HMatrix;
}

struct Translate {
    T:Direction,
}

impl Transform for Translate {
    fn get_matrix(&self) -> HMatrix {
        HMatrix{M:Matrix3x3{e:[[1.0,0.0,0.0],
                               [0.0,1.0,0.0],
                               [0.0,0.0,1.0]]},
                T:Direction{x:self.T.x,y:self.T.y,z:self.T.z}
        }
    }
}

struct UniformScale {
    S:f64,
}

impl Transform for UniformScale {
    fn get_matrix(&self) -> HMatrix {
        HMatrix{M:Matrix3x3{e:[[self.S,0.0,0.0],
            [0.0,self.S,0.0],
            [0.0,0.0,self.S]]},
            T:Direction{x:0.0,y:0.0,z:0.0}
        }
    }
}

struct Scale {
    S:Direction,
}

impl Transform for Scale {
    fn get_matrix(&self) -> HMatrix {
        HMatrix{M:Matrix3x3{e:[[self.S.x,0.0,0.0],
            [0.0,self.S.y,0.0],
            [0.0,0.0,self.S.z]]},
            T:Direction{x:0.0,y:0.0,z:0.0}
        }
    }
}

struct RotateX {
    S:f64,
}

impl Transform for RotateX {
    fn get_matrix(&self) -> HMatrix {
        let c=self.S.cos();
        let s=self.S.sin();
        HMatrix{M:Matrix3x3{e:[[1.0,0.0,0.0],
            [0.0,c,-s],
            [0.0,s,c]]},
            T:Direction{x:0.0,y:0.0,z:0.0}
        }
    }
}

struct RotateY {
    S:f64,
}

impl Transform for RotateY {
    fn get_matrix(&self) -> HMatrix {
        let c=self.S.cos();
        let s=self.S.sin();
        HMatrix{M:Matrix3x3{e:[[c,0.0,  s],
            [0.0,1.0,0.0],
            [ -s,0.0,c]]},
            T:Direction{x:0.0,y:0.0,z:0.0}
        }
    }
}

struct RotateZ {
    S:f64,
}

impl Transform for RotateZ {
    fn get_matrix(&self) -> HMatrix {
        let c=self.S.cos();
        let s=self.S.sin();
        HMatrix{M:Matrix3x3{e:[[c,-s,0.0],
            [ s,c,0.0],
            [0.0,0.0,1.0]]},
            T:Direction{x:0.0,y:0.0,z:0.0}
        }
    }
}

/// Get a single transformation matrix from a vector of transforms.
///
/// This is something that we can't do in C++: Add a function to an existing
/// class, without subclassing it etc.
impl Transform for Vec<Box<dyn Transform>> {
    fn get_matrix(&self) -> HMatrix {
        let this_transform=&self[0];
        let mut result =this_transform.get_matrix();
        for subtrans in self.iter() {
            result=subtrans.get_matrix().mul(&result);
        }
        result
    }
}

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
