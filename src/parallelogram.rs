#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize};

/// An parallelogram defined by two vectors. Degenerate parallelograms
/// are not guaranteed to work properly.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::prelude::*;
/// 
/// let gram1 = Parallelogram::new((0.0, 0.0), (10.0, 0.0), (0.0, 5.0));
/// let gram2 = Parallelogram::new((-1.0, -0.5), (2.0, 1.0), (-1.0, 1.0));
/// let gram3 = Parallelogram::new((8.0, 2.0), (4.0, 0.0), (1.0, 2.0));
/// 
/// let circle = Circle::new((-2.0, -2.0), 3.0);
/// 
/// assert!(sat_overlap(&gram1, &circle));
/// assert!(sat_overlap(&gram2, &gram1));
/// 
/// assert!(!sat_overlap(&gram3, &circle));
/// ```
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Parallelogram
{

    pub position: (f32, f32),
    pub u: (f32, f32),
    pub v: (f32, f32)
    
}

impl Parallelogram
{

    /// Create a new AABB at the given position with the given width and height.
    pub fn new(position: (f32, f32), u: (f32, f32), v: (f32, f32)) -> Parallelogram
    {

        return Parallelogram { position, u, v };

    }

    /// Creates a rectangular parallelogram so that it can be rotated later.
    pub fn rectangle(position: (f32, f32), width: f32, height: f32) -> Parallelogram
    {

        return Parallelogram { position, u: (width, 0.0), v: (0.0, height) };

    }

    pub fn points(&self) -> [(f32, f32); 4]
    {

        //TODO: Determine if this needs to be optimized or if the compiler does it for us
        return 
        [
            
            (0.0, 0.0),
            self.u,
            (self.u.0 + self.v.0, self.u.1 + self.v.1),
            self.v

        ];

    }

}

impl crate::Shape for Parallelogram
{

    fn position(&self) -> (f32, f32)
    {

        return self.position;

    }

    fn set_position(&mut self, position: (f32, f32))
    {

        self.position = position;

    }
    
    fn num_axes(&self) -> usize
    {

        return 2;

    }

    fn get_axis(&self, index: usize, _target: (f32, f32)) -> (f32, f32)
    {

        return match index
        {

            0 => (-self.u.1, self.u.0),
            _ => (-self.v.1, self.v.0)
            
        };

    }

    fn project(&self, axis: (f32, f32), _normalize: bool) -> (f32, f32)
    {

        return crate::project(self.position, axis, &self.points());

    }

    fn needs_closest(&self, _index: usize) -> bool
    {

        return false;

    }

    fn get_closest(&self, target: (f32, f32)) -> (f32, f32)
    {

        return crate::closest(self.position, target, &self.points());

    }

    fn point(&self, _index: usize) -> (f32, f32)
    {

        return self.position;

    }

}

impl crate::Rotate for Parallelogram
{

    fn rotate(&mut self, angle: f32)
    {

        let sin = f32::sin(angle);
        let cos = f32::cos(angle);

        self.rotate_sincos(sin, cos);

    }

    fn rotate_sincos(&mut self, sin: f32, cos: f32)
    {

        self.u = crate::rotate!(sin, cos, self.u);
        self.v = crate::rotate!(sin, cos, self.v);

    }

}

#[cfg(test)]
mod paralellogram_tests
{

    use super::*;
    use crate::{float_equal, Shape};

    #[test]
    fn test_num_axes()
    {

        let gram = Parallelogram::new((1.0, 2.0), (3.0, 1.0), (2.0, -2.0));

        assert_eq!(gram.num_axes(), 2);

    }

    #[test]
    fn test_get_axis()
    {

        let gram = Parallelogram::new((4.0, 10.0), (1.0, -2.0), (2.0, 3.0));

        let axis1 = gram.get_axis(0, (1.0, 0.0));
        let axis2 = gram.get_axis(1, (13.0, 20.0));

        assert!(float_equal(axis1.0, 2.0));
        assert!(float_equal(axis1.1, 1.0));
        assert!(float_equal(axis2.0, -3.0));
        assert!(float_equal(axis2.1, 2.0));

    }

    #[test]
    fn test_project()
    {

        let gram = Parallelogram::new((1.0, 2.0), (2.0, 0.0), (-1.0, 1.0));

        let axis1 = (1.0, 0.0);
        let axis2 = (1.0, -1.0);

        let projection1 = gram.project(axis1, true);
        let projection2 = gram.project(axis2, false);

        assert!(float_equal(projection1.0, 0.0));
        assert!(float_equal(projection1.1, 3.0));
        assert!(float_equal(projection2.0, -3.0));
        assert!(float_equal(projection2.1, 1.0));

    }

    #[test]
    fn test_needs_closest()
    {

        let gram = Parallelogram::new((1.0, 2.0), (3.0, 3.0), (2.0, 0.0));

        assert!(!gram.needs_closest(0));

    }

}