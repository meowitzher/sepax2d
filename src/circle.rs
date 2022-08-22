#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize};

/// A struct representing a circle via a position and radius.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::prelude::*;
/// 
/// let circle1 = Circle::new((0.0, 0.0), 2.0);
/// let circle2 = Circle::new((2.0, 2.0), 2.0);
/// 
/// let resolution = sat_collision(&circle1, &circle2);
/// let difference = 2.0 * (f32::sqrt(2.0) - 1.0); // 2 root 2 - 2
/// assert!(resolution.0 - difference < f32::EPSILON && resolution.0 - difference > -f32::EPSILON);
/// assert!(resolution.1 - difference < f32::EPSILON && resolution.1 - difference > -f32::EPSILON);
/// 
/// let polygon = Polygon::from_vertices((0.0, 3.0), vec![(0.0, 0.0), (1.0, 0.0), (1.0, -1.5), (0.0, -1.5)]);
/// assert!(sat_overlap(&polygon, &circle1));
/// assert!(sat_overlap(&circle1, &polygon));
/// ```
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Circle
{

    pub position: (f32, f32),
    pub radius: f32

}

impl Circle
{

    /// Create a new circle with given position and radius.
    pub fn new(position: (f32, f32), radius: f32) -> Circle
    {

        return Circle { position, radius };

    }

}

impl crate::Shape for Circle
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

        return 1;

    }

    fn get_axis(&self, _index: usize, target: (f32, f32)) -> (f32, f32)
    {

        return (target.0 - self.position.0, target.1 - self.position.1);

    }

    fn project(&self, axis: (f32, f32), normalize: bool) -> (f32, f32)
    {

        let projection = (self.position.0 * axis.0) + (self.position.1 * axis.1);
        
        //The projection of the circle along the axis is found by simply adding distance r
        //along it. However, if the vector is not a unit vector, then adding r equates to
        //adding r times the length.
        let mut magnitude = 1.0;
        if !normalize
        {

            magnitude = f32::sqrt((axis.0 * axis.0) + (axis.1 * axis.1));

        }

        return (projection - (self.radius * magnitude), projection + (self.radius * magnitude));

    }

    fn needs_closest(&self, _index: usize) -> bool
    {

        return true;

    }

    fn get_closest(&self, _target: (f32, f32)) -> (f32, f32)
    {

        return self.position;

    }

    fn point(&self, _index: usize) -> (f32, f32)
    {

        return self.position;

    }

}

impl crate::Rotate for Circle
{

    fn rotate(&mut self, _angle: f32)
    {

    }

    fn rotate_sincos(&mut self, _sin: f32, _cos: f32)
    {

    }

}

#[cfg(test)]
mod circle_tests
{

    use super::*;
    use crate::{float_equal, Shape};

    #[test]
    fn test_num_axes()
    {

        let circle = Circle::new((1.0, 2.0), 3.0);

        assert_eq!(circle.num_axes(), 1);

    }

    #[test]
    fn test_get_axis()
    {

        let circle = Circle::new((2.0, 3.0), 2.0);

        let axis1 = circle.get_axis(0, (1.0, 0.0));
        let axis2 = circle.get_axis(32, (13.0, 20.0));

        assert!(float_equal(axis1.0, -1.0));
        assert!(float_equal(axis1.1, -3.0));
        assert!(float_equal(axis2.0, 11.0));
        assert!(float_equal(axis2.1, 17.0));

    }

    #[test]
    fn test_project()
    {

        let circle = Circle::new((1.0, 2.0), 2.0);

        let axis1 = (1.0, 0.0);
        let axis2 = (1.0, 1.0);

        let projection1 = circle.project(axis1, true);
        let projection2 = circle.project(axis2, false);

        assert!(float_equal(projection1.0, -1.0));
        assert!(float_equal(projection1.1, 3.0));
        assert!(float_equal(projection2.0, 3.0 - 2.0 * f32::sqrt(2.0)));
        assert!(float_equal(projection2.1, 3.0 + 2.0 * f32::sqrt(2.0)));

    }

    #[test]
    fn test_needs_closest()
    {

        let circle = Circle::new((1.0, 2.0), 3.);

        assert!(circle.needs_closest(3));

    }

}