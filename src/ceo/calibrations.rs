use super::{Centroiding, Gmt, LensletArray, Source};
use ndarray::{stack, Array2, ArrayView, Axis, Ix2, ShapeBuilder};
use ndarray_linalg::svddc::{SVDDCInplace, UVTFlag};
use std::ops::Range;
//use std::time::Instant;

#[derive(Clone)]
pub enum Mirror {
    M1,
    M1MODES,
    M2,
}

#[derive(Clone)]
pub enum Segment {
    Txyz(f64, Option<Range<usize>>),
    Rxyz(f64, Option<Range<usize>>),
    Modes(f64, Range<usize>),
}
impl Segment {
    pub fn n_mode(&self) -> usize {
        match self {
            Segment::Txyz(_, i) => match i {
                Some(i) => i.end - i.start,
                None => 3,
            },
            Segment::Rxyz(_, i) => match i {
                Some(i) => i.end - i.start,
                None => 3,
            },
            Segment::Modes(_, n) => n.end - n.start,
        }
    }
    pub fn strip(&self) -> (f64, Range<usize>) {
        match self {
            Segment::Txyz(stroke, idx) => (
                *stroke,
                match idx.as_ref() {
                    Some(x) => Range {
                        start: x.start,
                        end: x.end,
                    },
                    None => 0..3,
                },
            ),
            Segment::Rxyz(stroke, idx) => (
                *stroke,
                match idx.as_ref() {
                    Some(x) => Range {
                        start: x.start + 3,
                        end: x.end + 3,
                    },
                    None => 3..6,
                },
            ),
            Segment::Modes(m, n) => (
                *m,
                Range {
                    start: n.start,
                    end: n.end,
                },
            ),
        }
    }
    pub fn range(&self) -> Range<usize> {
        self.strip().1
    }
    pub fn stroke(&self) -> f64 {
        self.strip().0
    }
}

pub struct Calibration {
    gmt: Gmt,
    src: Source,
    pub cog: Centroiding,
    n_side_lenslet: i32,
    lenslet_size: f64,
    m1_rbm: Vec<Vec<f64>>,
    m1_mode: Vec<Vec<f64>>,
    m2_rbm: Vec<Vec<f64>>,
    pub n_data: u32,
    pub n_mode: usize,
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
            m1_mode: vec![vec![]; 7],
            m2_rbm: vec![vec![0.; 6]; 7],
            n_data: 0,
            n_mode: 0,
        }
    }
    pub fn build(
        &mut self,
        zen: f32,
        azi: f32,
        valid_lenslets: &Vec<i8>,
        m1_n_mode: Option<usize>,
        m2_n_mode: Option<usize>,
    ) -> &mut Self {
        self.gmt.build(m1_n_mode.unwrap(), m2_n_mode);
        self.m1_mode = vec![vec![0.; m1_n_mode.or(Some(1)).unwrap()]; 7];
        self.src.build("R+I", vec![zen], vec![azi], vec![0f32]);
        self.cog.build(self.n_side_lenslet as u32, None);
        //        self.cog.process(detector, None);
        self.n_data = self
            .cog
            .set_valid_lenslets(None, Some(valid_lenslets.clone()));
        //println!("# valid lenslets: {}", self.n_data);
        self.n_data *= 2;
        self
    }
    pub fn calibrate(&mut self, mirror: Vec<Mirror>, segments: Vec<Vec<Segment>>) -> Vec<f32> {
        self.n_mode = 0; //14; //7 * t_or_r.len() as u32;
        let mut calibration: Vec<f32> =
            Vec::with_capacity((self.n_mode * 2 * self.cog.n_valid_lenslet as usize) as usize);
        for (k, segment) in segments.iter().enumerate() {
            for m in mirror.iter() {
                for rbm in segment.iter() {
                    let (stroke, idx) = rbm.strip();
                    for l in idx {
                        self.n_mode += 1;
                        calibration.extend::<Vec<f32>>(self.sample(k, m, l, stroke));
                    }
                }
            }
        }
        calibration
    }
    pub fn sample(&mut self, sid: usize, mirror: &Mirror, k: usize, stroke: f64) -> Vec<f32> {
        // PUSH
        match mirror {
            Mirror::M1 => {
                self.m1_rbm[sid][k] = stroke;
                self.gmt.update(Some(&self.m1_rbm), None, None);
            }
            Mirror::M1MODES => {
                self.m1_mode[sid][k] = stroke;
                self.gmt.update(None, None, Some(&self.m1_mode));
            }
            Mirror::M2 => {
                self.m2_rbm[sid][k] = stroke;
                self.gmt.update(None, Some(&self.m2_rbm), None);
            }
        }
        self.src.through(&mut self.gmt).xpupil().lenslet_gradients(
            self.n_side_lenslet,
            self.lenslet_size,
            &mut self.cog,
        );
        let c_push = self.cog.grab().valids(None);
        // PULL
        match mirror {
            Mirror::M1 => {
                self.m1_rbm[sid][k] = -stroke;
                self.gmt.update(Some(&self.m1_rbm), None, None);
            }
            Mirror::M1MODES => {
                self.m1_mode[sid][k] = -stroke;
                self.gmt.update(None, None, Some(&self.m1_mode));
            }
            Mirror::M2 => {
                self.m2_rbm[sid][k] = -stroke;
                self.gmt.update(None, Some(&self.m2_rbm), None);
            }
        }
        self.src.through(&mut self.gmt).xpupil().lenslet_gradients(
            self.n_side_lenslet,
            self.lenslet_size,
            &mut self.cog,
        );
        let c_pull = self.cog.grab().valids(None);
        // RESET
        match mirror {
            Mirror::M1 => {
                self.m1_rbm[sid][k] = 0f64;
                self.gmt.update(Some(&self.m1_rbm), None, None);
            }
            Mirror::M1MODES => {
                self.m1_mode[sid][k] = 0f64;
                self.gmt.update(None, None, Some(&self.m1_mode));
            }
            Mirror::M2 => {
                self.m2_rbm[sid][k] = 0f64;
                self.gmt.update(None, Some(&self.m2_rbm), None);
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
pub fn composite(    calibration: Vec<Vec<f32>>,
                     n_mode_vec: Vec<usize>,
                     composite_axis: Option<usize>,
) -> Array2<f32> {
    let mut a_view: Vec<ArrayView<f32, Ix2>> = Vec::new();
    for (c,n_mode) in calibration.iter().zip(n_mode_vec) {
        let n_data = c.len() / n_mode;
        //println!("Compositing IM: [{};{}]",n_data,n_mode);
        a_view.push(ArrayView::from_shape((n_data, n_mode).strides((1, n_data)), c).unwrap());
    }
    stack(Axis(composite_axis.or(Some(0)).unwrap()), &a_view).unwrap()
}
pub fn pseudo_inverse(
    d: &mut Array2<f32>,
    n_mode_vec: Vec<usize>,
    sig_n_thresholded: Option<usize>,
    composite_axis: Option<usize>,
) -> Array2<f32> {
    /*
    let mut a_view: Vec<ArrayView<f32, Ix2>> = Vec::new();
    for (c,n_mode) in calibration.iter().zip(n_mode_vec) {
        let n_data = c.len() / n_mode;
        println!("Compositing IM: [{};{}]",n_data,n_mode);
        a_view.push(ArrayView::from_shape((n_data, n_mode).strides((1, n_data)), c).unwrap());
    }
    let mut d = stack(Axis(composite_axis.or(Some(0)).unwrap()), &a_view).unwrap();
    */
    let (u, sig, v_t) = d.svddc_inplace(UVTFlag::Some).unwrap();
    //println!("eigen values: {}", sig);

    let mut i_sig = sig.mapv(|x| 1.0 / x);
    let n = i_sig.len();
    if sig_n_thresholded.is_some() {
        for k in 0..sig_n_thresholded.unwrap() {
            i_sig[n - 1 - k] = 0.0;
        }
    }

    let l_sv = Array2::from_diag(&i_sig);
    //print!("Computing the pseudo-inverse");
    //    let now = Instant::now();
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
