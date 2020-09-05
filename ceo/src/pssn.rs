use super::{Cu,Source};
use std::{fmt, mem};
use super::ceo_bindings::pssn as ceo_pssn;
use serde::ser::{Serialize, Serializer, SerializeStruct};

/// Wrapper for CEO PSSn
///
/// # Examples
///
/// ```
/// use gicsdom::ceo;
/// let mut src = ceo::Source::new(1,25.5,401);
/// src.build("V",vec![0.0],vec![0.0],vec![0.0]);
/// let mut gmt = ceo::Gmt::new();
/// gmt.build(27,None);
/// src.through(&mut gmt).xpupil();
/// println!("WFE RMS: {:.3}nm",src.wfe_rms_10e(-9)[0]);
/// let mut pssn = ceo::PSSn::new();
/// pssn.build(&mut src);
/// println!("PSSn: {:?}",pssn.reset(&mut src).estimates);
/// ```
///
// NEW PSSN
pub struct TelescopeError;
pub struct AtmosphereTelescopeError;
pub struct PSSn<S> {
    _c_: ceo_pssn,
    pub r0_at_zenith: f32,
    pub oscale: f32,
    pub zenith_angle: f32,
    /// PSSn estimates
    pub estimates: Vec<f32>,
    mode: std::marker::PhantomData<S>,
    otf: Vec<f32>,
}
impl<S> PSSn<S> {
    /// Creates a new `PSSn` with r0=16cm at zenith, L0=25m a zenith distance of 30 degrees
    pub fn new() -> PSSn<S> {
        PSSn {
            _c_: unsafe { mem::zeroed() },
            r0_at_zenith: 0.16,
            oscale: 25.0,
            zenith_angle: 30_f32.to_radians(),
            estimates: vec![],
            mode: std::marker::PhantomData,
            otf: Vec::new(),
        }
    }
    /// Creates a new `PSSn` from r0 at zenith and L0 a zenith distance of 30 degrees
    pub fn from_r0_and_outerscale(r0_at_zenith: f32, oscale: f32) -> PSSn<S> {
        PSSn {
            _c_: unsafe { mem::zeroed() },
            r0_at_zenith: r0_at_zenith,
            oscale: oscale,
            zenith_angle: 30_f32.to_radians(),
            estimates: vec![],
            mode: std::marker::PhantomData,
            otf: Vec::new(),
        }
    }
    /// Initializes PSSn atmosphere and telescope transfer function from a `Source` object
    pub fn build(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.setup(&mut src._c_, self.r0(), self.oscale);
        }
        self.estimates = vec![0.0; self._c_.N as usize];
        self
    }
    /// Integrates the `Source` optical transfer function
    pub fn accumulate(&mut self, src: &mut Source) {
        unsafe {
            self._c_.otf(&mut src._c_);
        }
    }
    /// Integrates the `Source` optical transfer function
    pub fn integrate(&mut self, src: &mut Source) {
        unsafe {
            self._c_.otf(&mut src._c_);
        }
    }
    /// Resets the `Source` optical transfer function to its initial value
    pub fn reset(&mut self) -> &mut Self {
        self._c_.N_O = 0;
        self
    }
    /// Computes `PSSn` spatial uniformity
    pub fn spatial_uniformity(&mut self) -> f32 {
        let mut pssn_values = self.estimates.clone();
        pssn_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        100. * ((pssn_values.len() as f32)
            * (*pssn_values.last().unwrap() - *pssn_values.first().unwrap()))
            / pssn_values.iter().sum::<f32>()
    }
    pub fn r0(&self) -> f32 {
        (self.r0_at_zenith.powf(-5_f32 / 3_f32) / self.zenith_angle.cos()).powf(-3_f32 / 5_f32)
    }
    pub fn xotf(&mut self) -> &Self {
        let mut d_otf  = Cu::vector(2*self._c_.NN as usize);
        d_otf.malloc();
        unsafe {
            self._c_.xotf(d_otf.as_ptr());
        }
        self.otf = d_otf.from_dev();
        self
    }
}
impl<T> Serialize for PSSn<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("PSSn", 1)?;
        state.serialize_field("r0",&self.r0())?;
        state.serialize_field("L0",&self.oscale)?;
        state.serialize_field("values",&self.estimates)?;
        state.serialize_field("otf",&self.otf)?;
        state.end()
    }
}
impl PSSn<TelescopeError> {
    /// Estimates the `PSSn` values
    pub fn peek(&mut self) -> &mut Self {
        unsafe { self._c_.eval1(self.estimates.as_mut_ptr()) }
        self
    }
}
impl PSSn<AtmosphereTelescopeError> {
    /// Estimates the `PSSn` values
    pub fn peek(&mut self) -> &mut Self {
        unsafe { self._c_.oeval1(self.estimates.as_mut_ptr()) }
        self
    }
}
impl<S> Drop for PSSn<S> {
    /// Frees CEO memory before dropping `PSSn`
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl<S> fmt::Display for PSSn<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}]",
            self.estimates
                .iter()
                .map(|x| format!("{:.4}", x))
                .collect::<Vec<String>>()
                .as_slice()
                .join(",")
        )
    }
}

