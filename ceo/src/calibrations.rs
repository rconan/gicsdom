use super::{
    cu::Single, shackhartmann::Geometric, Builder, Cu, Gmt, ShackHartmann, Source, GMT, SOURCE
};
use log;
use std::ops::Range;
//use std::time::Instant;

#[derive(Clone, Debug)]
/// GMT mirror functions
pub enum Mirror {
    /// M1 rigid body motion
    M1,
    /// M1 modal surface
    M1MODES,
    /// M2 rigid body motion
    M2,
    /// M1 modal surface
    M2MODES,
}

/// GMT segment functions
#[derive(Clone, Debug)]
pub enum Segment {
    /// Rigid body translations (stroke[m],range(`None`: 0..3))
    Txyz(f64, Option<Range<usize>>),
    /// Rigid body rotations (stroke[rd],range(`None`: 3..6))
    Rxyz(f64, Option<Range<usize>>),
    /// Modal surface coefficients
    Modes(f64, Range<usize>),
}
impl Segment {
    /// returns the number of rigid body motions based on the specified ranges
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
    /// returns the stroke and range
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
    /// returns the range
    pub fn range(&self) -> Range<usize> {
        self.strip().1
    }
    /// returns the stroke
    pub fn stroke(&self) -> f64 {
        self.strip().0
    }
}

/// GMT segment rigid body motion and surface figure calibration
///
/// `Calibration` creates its own GMT simulation with a `Gmt`, a `Source` and a `Builder`.
/// The calibration is performed by estimating the geometric centroids associated with the calibrated functions.
pub struct Calibration<T: Clone +  Builder<Component = ShackHartmann<Geometric>>> {
    gmt_blueprint: GMT,
    src_blueprint: SOURCE,
    wfs_blueprint: T,
    pub n_data: usize,
    pub n_mode: usize,
    pub poke: Cu<Single>,
    poke_qr: Cu<Single>,
}
impl<T: Clone + Builder<Component = ShackHartmann<Geometric>>> Calibration<T> {
    /// Creates a new `Calibration` with the blueprints of the `Gmt`, the `Source` and the `Builder`
    pub fn new(gmt: &Gmt, src: &Source, wfs_blueprint: T) -> Calibration<T>
    {
        Calibration {
            gmt_blueprint: gmt.into(),
            src_blueprint: src.into(),
            wfs_blueprint: wfs_blueprint,
            n_data: 0,
            n_mode: 0,
            poke: Cu::new(),
            poke_qr: Cu::new(),
        }
    }
    /// Performs the calibration of a single `Segment` function for a single `Mirror`
    pub fn sample(
        gmt: &mut Gmt,
        src: &mut Source,
        wfs: &mut ShackHartmann<Geometric>,
        sid: usize,
        mirror: &Mirror,
        k: usize,
        stroke: f64,
        valids: bool,
    ) -> Vec<f32> {
        let mut m1_rbm = vec![vec![0.; 6]; 7];
        let mut m1_mode = vec![vec![0.; gmt.m1_n_mode]; 7];
        let mut m2_rbm = vec![vec![0.; 6]; 7];
        let mut m2_mode = vec![vec![0.; gmt.m2_n_mode]; 7];
        let mut update = |x: f64, gmt: &mut Gmt| {
            match mirror {
                Mirror::M1 => {
                    m1_rbm[sid][k] = x;
                    gmt.update(Some(&m1_rbm), None, None, None);
                }
                Mirror::M1MODES => {
                    m1_mode[sid][k] = x;
                    gmt.update(None, None, Some(&m1_mode), None);
                }
                Mirror::M2 => {
                    m2_rbm[sid][k] = x;
                    gmt.update(None, Some(&m2_rbm), None, None);
                }
                Mirror::M2MODES => {
                    m2_mode[sid][k] = x;
                    gmt.update(None, None, None, Some(&m2_mode));
                }
            };
        };
        // PUSH
        update(stroke, gmt);
        wfs.reset();
        src.through(gmt).xpupil().through(wfs);
        wfs.process();
        let c_push = if valids {
            wfs.get_data().from_dev()
        } else {
            wfs.centroids.from_dev()
        };
        //println!("c push sum: {:?}", c_push.iter().sum::<f32>());
        // PULL
        update(-stroke, gmt);
        wfs.reset();
        src.through(gmt).xpupil().through(wfs);
        wfs.process();
        let c_pull = if valids {
            wfs.get_data().from_dev()
        } else {
            wfs.centroids.from_dev()
        };
        //println!("c pull sum: {:?}", c_pull.iter().sum::<f32>());
        // RESET
        update(0f64, gmt);
        // OUT
        c_push
            .iter()
            .zip(c_pull)
            .map(|x| 0.5 * (x.0 - x.1) / stroke as f32)
            .collect()
    }
    /// Calibrates the given mirror and segment functions:
    ///
    /// * `mirror`: `Vec` of `Mirror` functions
    /// * `segments`: a `Vec` the same size as the number of segment in the `mirror` with `Vec` elements of `Segment` functions
    pub fn calibrate(
        &mut self,
        mirror: Vec<Mirror>,
        segments: Vec<Vec<Segment>>,
        wfs_intensity_threshold: Option<f64>,
    ) {
        self.n_mode = segments
            .iter()
            .flat_map(|x| x.iter().map(|y| y.n_mode()))
            .sum::<usize>()
            * mirror.len();
        let mut calibration: Vec<f32> = vec![];
        let mut nnz = 0_usize;
        let wfs_intensity_threshold = wfs_intensity_threshold.unwrap_or(0.0);
        for (k, segment) in segments.iter().enumerate() {
            for m in mirror.iter() {
                for rbm in segment.iter() {
                    log::trace!("segment: {:?}", rbm);
                    let (stroke, idx) = rbm.strip();
                    let mut gmt = self.gmt_blueprint.clone().build();
                    let mut wfs = self.wfs_blueprint.clone().build();
                    let mut src = self.src_blueprint.clone().build();
                    src.through(&mut gmt).xpupil();
                    wfs.calibrate(&mut src, wfs_intensity_threshold);
                    nnz = wfs.n_valid_lenslet();
                    //println!("# valid lenslet: {}", wfs.n_valid_lenslet());
                    for l in idx {
                        calibration.extend::<Vec<f32>>(Calibration::<T>::sample(
                            &mut gmt,
                            &mut src,
                            &mut wfs,
                            k,
                            m,
                            l,
                            stroke,
                            wfs_intensity_threshold > 0f64,
                        ));
                    }
                }
            }
        }
        self.n_data = nnz * 2;
        log::info!("Calibration matrix: [{}x{}]", self.n_data, self.n_mode);
        self.poke = Cu::array(nnz * 2, self.n_mode);
        self.poke.to_dev(&mut calibration);
    }
    pub fn qr(&mut self) -> &mut Self {
        self.poke_qr = self.poke.clone();
        self.poke_qr.qr();
        self
    }
    pub fn solve(&mut self, data: &mut Cu<Single>) -> Cu<Single> {
        self.poke_qr.qr_solve(data)
    }
}
/*
/// reshapes vectors in matrices and stack the matrices together:
///
/// * `calibration` - `Vec` of vectors to reshape
/// * `n_mode_vec` - `Vec` of the matrices number of columns
/// * `composite_axis` - the axis along which to stack the matrices: row wise (0) or column wise (1)
pub fn composite(
    calibration: Vec<Vec<f32>>,
    n_mode_vec: Vec<usize>,
    composite_axis: Option<usize>,
) -> Array2<f32> {
    let mut a_view: Vec<ArrayView<f32, Ix2>> = Vec::new();
    for (c, n_mode) in calibration.iter().zip(n_mode_vec) {
        let n_data = c.len() / n_mode;
        //println!("Compositing IM: [{};{}]",n_data,n_mode);
        a_view.push(ArrayView::from_shape((n_data, n_mode).strides((1, n_data)), c).unwrap());
    }
    stack(Axis(composite_axis.or(Some(0)).unwrap()), &a_view).unwrap()
}
/// returns the pseudo-inverse but computing the singular value decomposition
///
/// * `d` - the matrix to inverse
/// * `sig_n_thresholded` - the number of filtered eigen values by increasing order
pub fn pseudo_inverse(d: &mut Array2<f32>, sig_n_thresholded: Option<usize>) -> Array2<f32> {
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
 */

/*
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gmt_calib_m2tt() {
        let pupil_size = 25.5f64;
        let n_lenslet = 48i32;
        let lenslet_size = pupil_size / n_lenslet as f64;
        let mut gmt = Gmt::new();
        gmt.build(0, None);
        let mut src = Source::new(1, pupil_size, n_lenslet * 16 + 1);
        src.build("V", vec![0.0], vec![0.0], vec![0.0]);
        src.set_fwhm(4.0);

        src.through(&mut gmt).xpupil();

        let mut cog = Centroiding::new();
        cog.build(n_lenslet as u32, None)
            .set_valid_lenslets(None, Some(src.masklet(n_lenslet as usize, 0.9)));
        src.lenslet_gradients(n_lenslet, lenslet_size, &mut cog);
        let s0 = cog.grab().valids(None);
        assert!(!s0.iter().any(|x| x.is_nan()));

        let lenslets = LensletArray {
            n_side_lenslet: n_lenslet,
            lenslet_size: lenslet_size,
        };
        let mut calib = Calibration::new(lenslets, None);
        calib.build(0.0, 0.0, &cog.valid_lenslets, Some(0), None);
        let mirror = vec![Mirror::M2];
        let segments = vec![vec![Segment::Rxyz(1e-6, Some(0..2))]; 7];
        let calibration = calib.calibrate(mirror, segments);

        let mut d = composite(vec![calibration], vec![14], Some(0));
        //        println!("{:?}",d.shape());
        let m = pseudo_inverse(&mut d, None);

        let rt = vec![vec![0f64, 0f64, 0f64, 1e-6, 1e-6, 0f64]; 7];
        gmt.update(None, Some(&rt), None);
        src.through(&mut gmt)
            .xpupil()
            .lenslet_gradients(n_lenslet, lenslet_size, &mut cog);
        let s = cog.grab().valids(None);
        //        println!("{:?}",src.segment_piston_10e(-9));
        assert!(!src.phase().iter().any(|x| x.is_nan()));

        let ds = s
            .iter()
            .zip(s0.iter())
            .map(|x| x.0 - x.1)
            .collect::<Vec<f32>>();
        let slopes = Array::from_shape_vec((ds.len(), 1), ds).unwrap();
        let m2_tt = 1e6 * m.dot(&slopes);
        //        println!("{:?}",m2_tt.into_shape((7,2)).unwrap());
        assert!(
            m2_tt.iter().all(|x| (x - 1.0).abs() < 1e-3),
            format!("{:?}", m2_tt.into_shape((7, 2)).unwrap())
        );
    }
}
*/
