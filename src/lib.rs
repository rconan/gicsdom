pub mod ceo;
pub mod astrotools;
pub mod agws;


pub struct Rotation {
    pub o: f64,
    pub axis: u8,
}
impl Rotation {
    pub fn new(o: f64, axis: u8) -> Rotation {
        Rotation { o, axis }
    }
    pub fn apply(&self, v: Vec<f64>) -> Vec<f64> {
        let (s, c) = self.o.sin_cos();
        let x = v[0];
        let y = v[1];
        let z = v[2];
        let rot_v: Vec<f64> = match self.axis {
            0 => vec![x, c * y + s * z, -s * y + c * z],
            1 => vec![c * x - s * z, y, s * x + c * z],
            2 => vec![c * x + s * y, -s * x + c * y, z],
            _ => {
                println!("Rotation axis must be one of: 0 for x, 1 for y or 2 for z!");
                vec![0.0; 3]
            }
        };
        rot_v
    }
}


