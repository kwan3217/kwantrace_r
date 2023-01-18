use crate::vector::{Direction, HMatrix, Matrix3x3};

pub trait Transform {
    fn get_matrix(&self)->HMatrix;
}

pub struct Translate {
    T:Direction,
}

impl Translate {
    pub(crate) fn make(x:f64, y:f64, z:f64) ->Translate {
        Translate{T:Direction{x,y,z}}
    }
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

pub type TransformList=Vec<Box<dyn Transform>>;

/// Get a single transformation matrix from a vector of transforms.
///
/// This is something that we can't do in C++: Add a function to an existing
/// class, without subclassing it etc.
impl Transform for TransformList {
    fn get_matrix(&self) -> HMatrix {
        let this_transform=&self[0];
        let mut result =HMatrix::identity();
        for subtrans in self.iter() {
            result=subtrans.get_matrix().mul(&result);
        }
        result
    }
}
