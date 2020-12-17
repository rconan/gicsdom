//!
//! # Analytic Ray Tracing
//!

use std::fmt;

pub type Vector = [f64; 3];

trait Arithmetic {
    fn dot(&self, other: &[f64]) -> f64;
    fn norm_square(&self) -> f64;
    fn norm(&self) -> f64;
    fn normalize(&mut self) -> Self;
}
impl Arithmetic for Vector {
    fn dot(&self, other: &[f64]) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }
    fn norm_square(&self) -> f64 {
        self.dot(self)
    }
    fn norm(&self) -> f64 {
        self.norm_square().sqrt()
    }
    fn normalize(&mut self) -> Self {
        let n = self.norm();
        self[0] /= n;
        self[1] /= n;
        self[2] /= n;
        *self
    }
}

/// # Ray definition
///
/// A ray is defined with:
///  - a point of origin: $\vec p = [x,y,z]$,
///  - a direction vector: $\vec u = [k,l,m]$ such as $\| \vec u \|=1$.
///
/// The ray tracing equation is given by: $$\vec{p^\prime} = \vec p + s \vec u,$$ where $s$ is the optical path length.
pub struct Ray {
    /// Ray point of origin
    pub p: Vector,
    /// Ray direction vector
    pub u: Vector,
}
/// # Ray builder
///
/// Build a new [`Ray`](crate::analytic::Ray)
pub struct NewRay {
    /// Ray point of origin
    pub p: Vector,
    /// Ray direction vector
    pub u: Vector,
}
impl Default for NewRay {
    fn default() -> Self {
        Self {
            p: [0f64; 3],
            u: [0f64, 0f64, -1f64],
        }
    }
}
impl NewRay {
    /// Build the [`Ray`](crate::analytic::Ray)
    pub fn build(self) -> Ray {
        Ray {
            p: self.p,
            u: self.u,
        }
    }
    /// Set the [`Ray`](crate::analytic::Ray) point of origin
    pub fn point_of_origin(self, p: Vector) -> Self {
        Self { p, ..self }
    }
    /// Set the [`Ray`](crate::analytic::Ray) direction vector
    pub fn direction_vector(self, u: Vector) -> Self {
        Self { u, ..self }
    }
    /// Set the [`Ray`](crate::analytic::Ray) direction vector from polar coordinates
    pub fn polar_direction_vector(self, z: f64, a: f64) -> Self {
        let ca = a.cos();
        let sa = a.sin();
        let sz = z.sin();
        let cz = z.cos();
        let u = [sz * ca, sz * sa, -cz].normalize();
        self.direction_vector(u)
    }
}
/// Create a [`NewRay`](crate::analytic::NewRay) at the origin propagate downward (z<0)
pub fn new_ray() -> NewRay {
    NewRay::default()
}
impl Ray {
    /// Compute the distance $s$ from the ray current location to [`Conic`](crate::analytic::Conic)
    /// We find the distance $s$ from:
    /// $$s=\frac{mR-\vec\alpha\cdot\vec\beta + \sqrt{(\vec\alpha\cdot\vec\beta-zR)^2-\|\vec\beta\|^2(\|\vec\alpha\|^2-2zR)}}{\|\vec\alpha\|^2}$$
    /// $$\vec\alpha = [x,y,z\sqrt{\kappa+1}]$$
    /// $$\vec\beta = [k,l,m\sqrt{\kappa+1}]$$
    pub fn distance_to(&self, conic: &Conic) -> f64 {
        let q = (conic.constant + 1f64).sqrt();
        let p: Vector = [
            self.p[0] - conic.origin[0],
            self.p[1] - conic.origin[1],
            self.p[2] - conic.origin[2],
        ];
        let alpha: Vector = [p[0], p[1], p[2] * q];
        let beta: Vector = [self.u[0], self.u[1], self.u[2] * q];
        let a = beta.norm_square();
        let b = 2f64 * (alpha.dot(&beta) - self.u[2] * conic.radius);
        let c = alpha.norm_square() - 2f64 * p[2] * conic.radius;
        0.5 * (-b + (b * b - 4f64 * a * c).sqrt()) / a
    }
    /// Trace ray from ray current position to [`Conic`](crate::analytic::Conic)
    pub fn trace_to(&mut self, conic: &Conic) {
        let s = self.distance_to(conic);
        self.p[0] += self.u[0] * s;
        self.p[1] += self.u[1] * s;
        self.p[2] += self.u[2] * s;
    }
    /// Solve ray tracing equation for $z$ given $x$ and $y$
    pub fn solve_for_z(&self, x: f64, y: f64) -> f64 {
        let x = x - self.p[0];
        let y = y - self.p[1];
        let num = x * x + y * y;
        let denom = self.u[0] * self.u[0] + self.u[1] * self.u[1];
        if denom < 1e-30 {
            std::f64::INFINITY
        } else {
            self.p[2] + self.u[2] * (num / denom).sqrt()
        }
    }
}
impl fmt::Display for Ray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "P: [{:+15.9},{:+15.9},{:+15.9}] ; U: [{:+.6},{:+.6},{:+.6}]",
            self.p[0], self.p[1], self.p[2], self.u[0], self.u[1], self.u[2],
        )
    }
}

/// # Conic definition
///
/// A conic surface is defined by the set of coordinates $(x,y,z)$ that satisfies
/// $$ F(x,y,z) = r^2 - 2zR + z^2(\kappa+1)=0,$$
/// where $r^2=x^2+y^2$, $R$ is the radius of curvature and $\kappa$ is the conic constant.
pub struct Conic {
    /// Conic constant $\kappa$
    pub constant: f64,
    /// Radius of curvature
    pub radius: f64,
    /// Origin vector
    pub origin: Vector,
}
impl Conic {
    /// Creates a new `Conic`
    pub fn new(constant: f64, radius: f64) -> Self {
        Self {
            constant,
            radius,
            origin: [0f64; 3],
        }
    }
    /// Creates a new `Conic` with GMT M1 prescription:
    ///  - $\kappa=-0.9982857$
    ///  - $R=36$
    pub fn gmt_m1() -> Self {
        Self::new(-0.9982857, 36.0)
    }
    /// Creates a new `Conic` with GMT M2 prescription:
    ///  - $\kappa=-0.71692784$
    ///  - $R=-4.1639009$
    pub fn gmt_m2() -> Self {
        Self {
            constant: -0.71692784,
            radius: -4.1639009,
            origin: [0f64, 0f64, 20.26247614],
        }
    }
    /// Solve conic surface $F(x,y,z)$ for z :
    pub fn height_at(&self, v: Vector) -> Vector {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        _v[2] = self.c() * r2 / (1f64 + self.sqrt_(r2));
        _v
    }
    /// Conic $x$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial x} = -\frac{x}{R+\sqrt{R^2-(\kappa+1)r^2}}
    /// $$
    ///
    pub fn x_partial_at(&self, v: Vector) -> f64 {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        -self.c() * _v[0] / self.sqrt_(r2)
    }
    /// Conic $y$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial y} = -\frac{y}{R+\sqrt{R^2-(\kappa+1)r^2}}
    /// $$
    ///
    pub fn y_partial_at(&self, v: Vector) -> f64 {
        let mut _v = v;
        _v[2] = 0f64;
        let r2 = _v.norm_square();
        -self.c() * _v[1] / self.sqrt_(r2)
    }
    /// Conic $z$ partial derivative
    /// $$
    /// \frac{\partial{F}}{\partial z} = 1
    /// $$
    ///
    pub fn z_partial_at(&self, _v: Vector) -> f64 {
        1f64
    }
    /// normal vector to conic surface
    pub fn normal_at(&self, v: Vector) -> Vector {
        [
            self.x_partial_at(v),
            self.y_partial_at(v),
            self.z_partial_at(v),
        ]
        .normalize()
    }
    /// Reflect ray from conic surface
    pub fn reflect(&self, ray: &mut Ray) {
        let n = self.normal_at(ray.p);
        let q = 2f64 * ray.u.dot(&n);
        ray.u[0] -= q * n[0];
        ray.u[1] -= q * n[1];
        ray.u[2] -= q * n[02];
        ray.u = ray.u.normalize();
    }
    fn kp1(&self) -> f64 {
        self.constant + 1f64
    }
    fn c(&self) -> f64 {
        1f64 / self.radius
    }
    fn sqrt_(&self, r2: f64) -> f64 {
        (1f64 - self.kp1() * self.c() * self.c() * r2).sqrt()
    }
}
