/// Position vector, a 3-element vector which participates in translation.
///
///When thinking about the vector as a homogeneous vector, all
///vectors of this type are considered to have an implicit w=1 component.
///This is handled in the matrix multiplication of a position vector.
pub struct Position {
    pub x:f64,
    pub y:f64,
    pub z:f64,
}

/// Direction vector, a 3-element vector which does not participate in translation.
///
///When thinking about the vector as a homogeneous vector, all
///vectors of this type are considered to have an implicit w=1 component.
///This is handled in the matrix multiplication of a position vector.
pub struct Direction {
    pub x:f64,
    pub y:f64,
    pub z:f64,
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
pub struct Matrix3x3 {
    pub(crate) e:[[f64;3];3],
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

pub struct HMatrix {
    pub(crate) M:Matrix3x3,
    pub(crate) T:Direction,
}

impl HMatrix {
    pub(crate) fn mul(&self, rhs: &HMatrix) -> HMatrix {
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

///Operations which are common to both Position and Direction, and operations
///between them.
pub(crate) trait Vector {
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

pub struct Ray {
    pub(crate) r0:Position,
    pub(crate) v:Direction,
}
