use super::{Centroiding, Gmt, Source};
use crate::agws::probe::LensletArray;
use ndarray::{Axis,Array2, ArrayView, ShapeBuilder,stack,Ix2};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::time::Instant;

pub enum Mirror {
    M1,
    M2,
}

pub enum RigidBodyMotion {
    Txyz(usize, f64),
    Rxyz(usize, f64),
}

pub struct Calibration {
    gmt: Gmt,
    src: Source,
    pub cog: Centroiding,
    n_side_lenslet: i32,
    lenslet_size: f64,
    m1_rbm: Vec<Vec<f64>>,
    m2_rbm: Vec<Vec<f64>>,
    pub n_data: u32,
    pub n_mode: u32,
}

impl Calibration {
    pub fn new(optics: LensletArray, n_px_lenslet: Option<i32>) -> Calibration {
        Calibration {
            gmt: Gmt::new(),
            src: Source::new(
                1,
                optics.n_side_lenslet as f64 * optics.lenslet_size,
                optics.n_side_lenslet * n_px_lenslet.or(Some(16)).unwrap() + 1,
            ),
            cog: Centroiding::new(),
            n_side_lenslet: optics.n_side_lenslet,
            lenslet_size: optics.lenslet_size,
            m1_rbm: vec![vec![0.; 6]; 7],
            m2_rbm: vec![vec![0.; 6]; 7],
            n_data: 0,
            n_mode: 0,
        }
    }
    pub fn build(
        &mut self,
        valid_lenslets: &Vec<i8>,
        m1_n_mode: Option<u64>,
        m2_n_mode: Option<u64>,
    ) -> &mut Self {
        self.gmt.build(m1_n_mode.unwrap(), m2_n_mode);
        self.src.build("V", vec![0f32], vec![0f32], vec![0f32]);
        self.cog.build(self.n_side_lenslet as u32, None);
        //        self.cog.process(detector, None);
        self.n_data = self
            .cog
            .set_valid_lenslets(None, Some(valid_lenslets.clone()));
        //println!("# valid lenslets: {}", self.n_data);
        self.n_data *= 2;
        self
    }
    pub fn calibrate(&mut self, mirror: Vec<Mirror>, t_or_r: Vec<RigidBodyMotion>) -> Vec<f32> {
        self.n_mode = 7 * t_or_r.len() as u32;
        let mut calibration: Vec<f32> =
            Vec::with_capacity((self.n_mode * 2 * self.cog.n_valid_lenslet) as usize);
        for k in 0..7 {
            for m in mirror.iter() {
                for rbm in t_or_r.iter() {
                    calibration.extend::<Vec<f32>>(self.sample(k, m, rbm));
                }
            }
        }
        calibration
    }
    pub fn sample(&mut self, sid: usize, mirror: &Mirror, t_or_r: &RigidBodyMotion) -> Vec<f32> {
        let (idx, stroke) = match t_or_r {
            RigidBodyMotion::Txyz(idx, stroke) => (*idx, *stroke),
            RigidBodyMotion::Rxyz(idx, stroke) => (*idx + 3, *stroke),
        };
        // PUSH
        match mirror {
            Mirror::M1 => {
                self.m1_rbm[sid][idx] = stroke;
            }
            Mirror::M2 => {
                self.m2_rbm[sid][idx] = stroke;
            }
        }
        self.gmt.update(Some(&self.m1_rbm), Some(&self.m2_rbm));
        self.src.through(&mut self.gmt).xpupil().lenslet_gradients(
            self.n_side_lenslet,
            self.lenslet_size,
            &mut self.cog,
        );
        let c_push = self.cog.grab().valids(None);
        // PULL
        match mirror {
            Mirror::M1 => {
                self.m1_rbm[sid][idx] = -stroke;
            }
            Mirror::M2 => {
                self.m2_rbm[sid][idx] = -stroke;
            }
        }
        self.gmt.update(Some(&self.m1_rbm), Some(&self.m2_rbm));
        self.src.through(&mut self.gmt).xpupil().lenslet_gradients(
            self.n_side_lenslet,
            self.lenslet_size,
            &mut self.cog,
        );
        let c_pull = self.cog.grab().valids(None);
        // RESET
        match mirror {
            Mirror::M1 => {
                self.m1_rbm[sid][idx] = 0f64;
            }
            Mirror::M2 => {
                self.m2_rbm[sid][idx] = 0f64;
            }
        }
        // OUT
        c_push
            .iter()
            .zip(c_pull)
            .map(|x| 0.5 * (x.0 - x.1) / stroke as f32)
            .collect()
    }
}
pub fn pseudo_inverse(calibration: Vec<Vec<f32>>, n_mode: usize) -> Array2<f32> {
    let mut a_view: Vec<ArrayView<f32,Ix2>> = Vec::new();
    for c in calibration.iter() {
        let n_data = c.len()/n_mode;
        a_view.push(ArrayView::from_shape((n_data, n_mode).strides((1, n_data)), c).unwrap());
    }
    let mut d = stack(Axis(0),&a_view).unwrap();
    let (u, sig, v_t) = d.svddc_inplace(UVTFlag::Some).unwrap();
    //println!("eigen values: {}", sig);

    let i_sig = sig.mapv(|x| 1.0 / x);

    let l_sv = Array2::from_diag(&i_sig);
    //print!("Computing the pseudo-inverse");
    let now = Instant::now();
    let m: Array2<f32> = v_t.unwrap().t().dot(&l_sv.dot(&u.unwrap().t()));
    //println!(" in {}ms", now.elapsed().as_millis());
    //self.reconstructor = m.into_dimensionality::<Ix1>().unwrap().to_vec();
    //println!("reconstructor shape: {:?}",m.shape());
    m.to_owned()
}
/*
pub fn dot(&self, data: Vec<f32>) -> Array2<f32> {
    let n = 2 * self.cog.n_valid_lenslet as usize;
    let m = Array2::from_shape_vec((n,self.n_modes as usize),self.reconstructor.clone());
    let v = Array::from_shape_vec((data.len(), 1), data).unwrap();
    let c = m.to_owned().unwrap().dot(&v);
    c
}
*/
