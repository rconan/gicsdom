use std::{f32, mem};

use super::ceo_bindings::{dev2host, imaging};
use super::Propagation;
use super::{Centroiding, Gmt, Source,Conversion};

#[derive(Copy, Clone)]
/// A square lenslet array
pub struct LensletArray {
    /// The number of lenslet per side
    pub n_side_lenslet: i32,
    /// Dimension [m] of one lenslet
    pub lenslet_size: f64,
}
#[derive(Copy, Clone)]
/// Detector noise specifications
pub struct NoiseDataSheet {
    /// Read-out noise rms
    pub rms_read_out_noise: f64,
    /// Number of background photons
    pub n_background_photon: f64,
    /// Noise factor
    pub noise_factor: f64,
}
impl NoiseDataSheet {
    fn rms_read_out_noise(rms_read_out_noise: f64) -> Self {
        NoiseDataSheet {
            rms_read_out_noise,
            n_background_photon: 0f64,
            noise_factor: 1f64,
        }
    }
}
impl Default for NoiseDataSheet {
    fn default() -> Self {
        NoiseDataSheet {
            rms_read_out_noise: 0f64,
            n_background_photon: 0f64,
            noise_factor: 1f64,
        }
    }
}

/// An optical imager with a detector
///
/// The optical imager is a square lenslet array which focal plane lies on the detector.
/// The detector continuously integrates the images formed on the detector until it is explicitly reset
pub struct Imaging {
    _c_: imaging,
    dft_osf: i32,
}
impl Imaging {
    /// Creates a new `Imaging`
    pub fn new() -> Imaging {
        Imaging {
            _c_: unsafe { mem::zeroed() },
            dft_osf: 1,
        }
    }
    /// Set `Imaging` parameters
    ///
    /// * `n_sensor` - the number of `Imaging` sensor
    /// * `n_side_lenslet` - the size of the square lenslet array
    /// * `n_px_lenslet` - the number of pixel per lenslet, for a total resolution of (`n_side_lenslet`X`n_px_lenslet`+1)^2
    /// * `dft_osf` - the discrete Fourier transform oversampling factor
    /// * `n_px_imagelet` - the sampling of a lenslet focal plane image
    /// * `binning` - binning factor of a imagelet
    pub fn build(
        &mut self,
        n_sensor: i32,
        n_side_lenslet: i32,
        n_px_lenslet: i32,
        dft_osf: i32,
        n_px_imagelet: i32,
        binning: i32,
    ) -> &mut Self {
        unsafe {
            self._c_.setup3(
                n_px_lenslet,
                n_side_lenslet,
                dft_osf,
                n_px_imagelet,
                binning,
                n_sensor,
            );
        }
        self.dft_osf = dft_osf;
        self
    }
    pub fn __ceo__(&self) -> &imaging {
        &self._c_
    }
    /// Returns the frame from the GPU
    pub fn frame_transfer(&mut self, frame: &mut Vec<f32>) -> &mut Self {
        unsafe {
            dev2host(
                frame.as_mut_ptr(),
                self._c_.d__frame,
                self.resolution() * self.resolution() * self._c_.N_SOURCE,
            );
        }
        self
    }
    /// Resets the detector frame to zero
    pub fn reset(&mut self) -> &mut Self {
        unsafe {
            self._c_.reset();
        }
        self
    }
    /// Returns the detector resolution
    pub fn resolution(&self) -> i32 {
        self._c_.N_PX_CAMERA * self._c_.N_SIDE_LENSLET
    }
    /// Return the number of frames that have been summed since the last reset
    pub fn n_frame(&self) -> u32 {
        self._c_.N_FRAME as u32
    }
    /// Reads out the detector by adding noise to the frame if a `NoiseDataSheet` is passed and the intensity is scaled according to the detector `exposure` time [s]
    pub fn readout(
        &mut self,
        exposure: f64,
        detector_noise_properties: Option<NoiseDataSheet>,
    ) -> &mut Self {
        detector_noise_properties.map_or_else(
            || (),
            |p| unsafe {
                self._c_.readout1(
                    exposure as f32,
                    p.rms_read_out_noise as f32,
                    p.n_background_photon as f32,
                    p.noise_factor as f32,
                );
            },
        );
        self
    }
    /// Sets the pixel scale
    pub fn set_pixel_scale(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.pixel_scale = (src.wavelength() / src.pupil_size) as f32
                * (self._c_.N_SIDE_LENSLET * self._c_.BIN_IMAGE / self.dft_osf) as f32;
        }
        self
    }
    /// Returns the detector pixel scale
    pub fn pixel_scale(&mut self, src: &mut Source) -> f32 {
        (src.wavelength() / src.pupil_size) as f32
            * (self._c_.N_SIDE_LENSLET * self._c_.BIN_IMAGE / self.dft_osf) as f32
    }
    /// Sets the detector pointing direction
    pub fn set_pointing(&mut self, mut zen: Vec<f32>, mut azi: Vec<f32>) -> &mut Self {
        unsafe {
            self._c_.absolute_pointing = 1;
            self._c_
                .set_pointing_direction(zen.as_mut_ptr(), azi.as_mut_ptr());
        }
        self
    }
}
impl Drop for Imaging {
    /// Frees CEO memory before dropping `Imaging`
    fn drop(&mut self) {
        unsafe {
            self._c_.cleanup();
        }
    }
}
impl Propagation for Imaging {
    /// Fourier propagates the wavefront to the focal plane onto the detector
    fn propagate(&mut self, src: &mut Source) -> &mut Self {
        unsafe {
            self._c_.propagate(&mut src._c_);
        }
        self
    }
    fn time_propagate(&mut self, _secs: f64, src: &mut Source) -> &mut Self {
        self.propagate(src)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// Test the intensity per lenslet
    fn imaging_flux() {
        let pupil_size = 25.5f64;
        let n_side_lenslet = 20;
        let n_px_lenslet = 32;
        let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;
        let lenslet_size = (pupil_size / n_side_lenslet as f64) as f32;
        let mut gmt = Gmt::new();
        gmt.build(0, None);
        let mut src = Source::new(1, pupil_size, pupil_sampling);
        src.build("V", vec![0f32], vec![0f32], vec![18f32]);
        let mut sensor = Imaging::new();
        sensor.build(1, n_side_lenslet, n_px_lenslet, 2, 2 * n_px_lenslet, 1);
        let mut cog = Centroiding::new();
        cog.build(n_side_lenslet as u32, None);

        src.through(&mut gmt).xpupil().through(sensor.reset());
        cog.process(&sensor, None);
        let fluxlet = cog
            .lenslet_flux()
            .iter()
            .cloned()
            .fold(-f32::INFINITY, f32::max);
        println!("Sensor lenslet flux: {}", fluxlet);
        let fluxlet_expected = src.n_photon()[0] * lenslet_size * lenslet_size;
        println!("Lenslet expected flux: {}", fluxlet_expected);
        assert!((fluxlet - fluxlet_expected).abs() / fluxlet_expected < 1e-1);
    }

    #[test]
    fn imaging_noise_photon() {
        let pupil_size = 25.5f64;
        let n_side_lenslet = 40;
        let n_px_lenslet = 16;
        let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;
        let lenslet_size = (pupil_size / n_side_lenslet as f64) as f32;
        let mut gmt = Gmt::new();
        gmt.build(0, None);
        let mut src = Source::new(1, pupil_size, pupil_sampling);
        src.build("V", vec![0f32], vec![0f32], vec![18f32]);
        let mut sensor = Imaging::new();
        sensor.build(1, n_side_lenslet, n_px_lenslet, 2, 2*n_px_lenslet, 1);
        let p = sensor.pixel_scale(&mut src) as f64;

        let mut cog0 = Centroiding::new();
        cog0.build(n_side_lenslet as u32, None);
        src.through(&mut gmt).xpupil().through(&mut sensor);
        let nv = cog0
            .process(&sensor, None)
            .set_valid_lenslets(Some(0.9), None);
        println!("Valid lenslet #: {}",nv);

        let mut cog = Centroiding::new();
        cog.build(n_side_lenslet as u32, Some(p));
        src.through(&mut gmt).xpupil().through(&mut sensor);
        src.through(&mut gmt).xpupil().through(&mut sensor);
        sensor.readout(1f64, Some(NoiseDataSheet::default()) );
        let s = cog
            .process(&sensor, Some(&cog0))
            .grab()
            .valids(Some(&cog0.valid_lenslets));
        println!("Pixel scale: {}mas",p.to_mas());
        let m = s.iter().sum::<f32>() / nv as f32;
        let v = s.iter().map(|x| (x - m).powi(2)).sum::<f32>() / nv as f32;
        println!("Centroid rms error: {}",(v.sqrt() as f64).to_mas());

        let fluxlet = cog
            .lenslet_flux()
            .iter()
            .cloned()
            .fold(-f32::INFINITY, f32::max);
        let fluxlet_expected = src.n_photon()[0] * lenslet_size * lenslet_size;
        println!("flux ratio: {}",fluxlet/fluxlet_expected);
        let fwhm = 1.03*src.wavelength()/lenslet_size as f64;
        println!("FWHM: {}mas",fwhm.to_mas());
        let v_expected = fwhm.powi(2)/(2f64*2f64.ln()*fluxlet_expected as f64);
        println!("Expected centroid rms error: {}",v_expected.sqrt().to_mas());

        assert!((v as f64-v_expected).abs()/v_expected<1e-1);
    }

    #[test]
    fn imaging_noise_readout() {
        let n_side_lenslet = 40;
        let n_px_lenslet = 16;

        let mut sensor = Imaging::new();
        sensor.build(1, n_side_lenslet, n_px_lenslet, 2, 2 * n_px_lenslet, 1);

        sensor
            .reset()
            .readout(1f64, Some(NoiseDataSheet::rms_read_out_noise(1f64)));
        let n = sensor.resolution().pow(2);
        let mut frame = vec![0f32; n as usize];
        sensor.frame_transfer(&mut frame);

        let m = frame.iter().sum::<f32>() / n as f32;
        let v = frame.iter().map(|x| (x - m).powi(2)).sum::<f32>() / n as f32;
        println!("RON: {}", v.sqrt());
        assert!((1f32 - v.sqrt()).abs() < 1e-2);
    }

    #[test]
    fn imaging_pointing() {
        let pupil_size = 25.5f64;
        let n_side_lenslet = 1;
        let n_px_lenslet = 511;
        let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;
//        let lenslet_size = (pupil_size / n_side_lenslet as f64) as f32;
        let mut gmt = Gmt::new();
        gmt.build(0, None);
        let mut src = Source::new(1, pupil_size, pupil_sampling);
        src.build("V", vec![0f32], vec![0f32], vec![18f32]);
        let mut sensor = Imaging::new();
        sensor.build(1, n_side_lenslet, n_px_lenslet, 2, 128, 1);
        let p = src.wavelength()/pupil_size/2f64;
        println!("Pixel scale: {}mas",p.to_mas());

        let mut cog0 = Centroiding::new();
        cog0.build(n_side_lenslet as u32, None);
        src.through(&mut gmt).xpupil().through(&mut sensor);
        let nv = cog0
            .process(&sensor, None)
            .set_valid_lenslets(Some(0.9), None);
        println!("Valid lenslet #: {}",nv);

        let mut cog = Centroiding::new();
        cog.build(n_side_lenslet as u32, Some(p));

        sensor.set_pixel_scale(&mut src);
        sensor.set_pointing(vec![1e-9], vec![0.0]);
        src.through(&mut gmt).xpupil().through(&mut sensor);

        let s = cog
            .process(&sensor, Some(&cog0))
            .grab()
            .valids(Some(&cog0.valid_lenslets));
        println!("s: {:?}",s);

    }

    #[test]
    fn imaging_exposure() {
        let pupil_size = 25.5f64;
        let n_side_lenslet = 1;
        let n_px_lenslet = 511;
        let pupil_sampling = n_side_lenslet * n_px_lenslet + 1;
        //        let lenslet_size = (pupil_size / n_side_lenslet as f64) as f32;
        let mut gmt = Gmt::new();
        gmt.build(0, None);
        let mut src = Source::new(1, pupil_size, pupil_sampling);
        src.build("V", vec![0f32], vec![0f32], vec![18f32]);
        let mut sensor = Imaging::new();
        sensor.build(1, n_side_lenslet, n_px_lenslet, 2, 128, 1);

        let mut frame = vec![0f32;sensor.resolution().pow(2) as usize];
        src.through(&mut gmt).xpupil().through(&mut sensor);
        sensor.frame_transfer(&mut frame);
        let f0 =  frame.iter().sum::<f32>();

        src.through(&mut gmt).xpupil().through(&mut sensor);
        sensor.frame_transfer(&mut frame);
        let f =  frame.iter().sum::<f32>();
        println!("Flux ratio: {}",f/f0);
        assert_eq!((f/f0)as usize, 2);

        src.through(&mut gmt).xpupil().through(sensor.reset());
        sensor.frame_transfer(&mut frame);
        let f =  frame.iter().sum::<f32>();
        println!("Flux ratio: {}",f/f0);
        assert_eq!((f/f0)as usize, 1);

        sensor.reset();
        for _ in 0..10 {
            src.through(&mut gmt).xpupil().through(&mut sensor);
        }
        sensor.frame_transfer(&mut frame);
        let f =  frame.iter().sum::<f32>();
        println!("Flux ratio: {}",(f/f0) as usize);
        assert_eq!((f/f0)as usize, 10);
    }
}
