/// An axis-aligned bounding box, that is a rectangle aligned
/// along the Cartesian coordinate system.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::aabb::AABB;
/// use sepax2d::circle::Circle;
/// use sepax2d::{sat_collision, sat_overlap};
/// 
/// let box1 = AABB::new((0.0, 0.0), 10.0, 5.0);
/// let box2 = AABB::new((-1.0, -0.5), 2.0, 1.0);
/// let box3 = AABB::new((8.0, 2.0), 4.0, 2.0);
/// 
/// let circle = Circle::new((-2.0, -2.0), 3.0);
/// 
/// assert!(sat_overlap(&box1, &circle));
/// assert!(sat_overlap(&box2, &box1));
/// 
/// let resolution = sat_collision(&box1, &box3);
/// assert!(resolution.0 - 2.0 < f64::EPSILON && resolution.0 - 2.0 > -f64::EPSILON);
/// assert!(resolution.1 - 0.0 < f64::EPSILON && resolution.1 - 0.0 > -f64::EPSILON);
/// ```
pub struct AABB
{

    pub position: (f64, f64),
    pub width: f64,
    pub height: f64
    

}

impl AABB
{

    /// Create a new AABB at the given position with the given width and height.
    pub fn new(position: (f64, f64), width: f64, height: f64) -> AABB
    {

        return AABB { position, width, height };

    }

    fn points(&self) -> [(f64, f64); 4]
    {

        //TODO: Determine if this needs to be optimized or if the compiler does it for us
        return 
        [
            
            (0.0, 0.0),
            (self.width, 0.0),
            (self.width, self.height),
            (0.0, self.height)

        ];

    }

}

impl crate::Shape for AABB
{

    fn position(&self) -> (f64, f64)
    {

        return self.position;

    }

    fn num_axes(&self) -> usize
    {

        return 2;

    }

    fn get_axis(&self, index: usize, _target: (f64, f64)) -> (f64, f64)
    {

        return match index
        {

            0 => (1.0, 0.0),
            _ => (0.0, 1.0)
            
        };

    }

    fn project(&self, axis: (f64, f64), _normalize: bool) -> (f64, f64)
    {

        return crate::project(self.position, axis, &self.points());

    }

    fn needs_closest(&self, _index: usize) -> bool
    {

        return false;

    }

    fn get_closest(&self, target: (f64, f64)) -> (f64, f64)
    {

        return crate::closest(self.position, target, &self.points());

    }

    fn point(&self, _index: usize) -> (f64, f64)
    {

        return self.position;

    }

}

#[cfg(test)]
mod aabb_tests
{

    use super::*;
    use crate::{float_equal, Shape};

    #[test]
    fn test_num_axes()
    {

        let aabb = AABB::new((1.0, 2.0), 3.0, 2.0);

        assert_eq!(aabb.num_axes(), 2);

    }

    #[test]
    fn test_get_axis()
    {

        let aabb = AABB::new((4.0, 10.0), 2.0, 5.0);

        let axis1 = aabb.get_axis(0, (1.0, 0.0));
        let axis2 = aabb.get_axis(1, (13.0, 20.0));

        assert!(float_equal(axis1.0, 1.0));
        assert!(float_equal(axis1.1, 0.0));
        assert!(float_equal(axis2.0, 0.0));
        assert!(float_equal(axis2.1, 1.0));

    }

    #[test]
    fn test_project()
    {

        let aabb = AABB::new((1.0, 2.0), 3.0, 4.0);

        let axis1 = (1.0, 0.0);
        let axis2 = (1.0, -1.0);

        let projection1 = aabb.project(axis1, true);
        let projection2 = aabb.project(axis2, false);

        assert!(float_equal(projection1.0, 1.0));
        assert!(float_equal(projection1.1, 4.0));
        assert!(float_equal(projection2.0, -5.0));
        assert!(float_equal(projection2.1, 2.0));

    }

    #[test]
    fn test_needs_closest()
    {

        let aabb = AABB::new((1.0, 2.0), 3.0, 2.0);

        assert!(!aabb.needs_closest(1));

    }

}