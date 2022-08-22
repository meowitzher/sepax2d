#[cfg(feature = "serde")]
use serde::{Serialize, ser::SerializeStruct, Serializer, Deserialize, Deserializer};

use crate::circle::Circle;

/// A struct representing a capsule, i.e. a rotated rectangle capped by half circles.
/// The position is located in the center of the rectangle, with the arm vector denoting
/// half of the capsule's length. Radius represents half of the capsule's width, which
/// corresponds to the radius of the half circles.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::prelude::*;
/// 
/// let capsule = Capsule::new((0.0, 0.0), (0.0, 1.0), 2.0);
/// //A capsule formed from a rectangle with vertices (-1, 1), (1,1), (-1,-1), and (1,-1)
/// 
/// let hexagon = Polygon::from_vertices((1.5, 0.0), vec![(0.0, 0.0), (0.0, 1.0), (1.0, 2.0), (2.0, 1.0), (2.0, 0.0)]);
/// let square = AABB::new((2.0, 3.0), 4.0, 0.5);
/// 
/// assert!(sat_overlap(&capsule, &hexagon));
/// assert!(!sat_overlap(&square, &capsule));
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Capsule
{

    pub position: (f32, f32),
    arm: (f32, f32),
    perp: (f32, f32),
    pub radius: f32

}

impl Capsule
{

    /// Create a new capsule with the given center position, arm, and radius. 
    pub fn new(position: (f32, f32), arm: (f32, f32), radius: f32) -> Capsule 
    {

        return Capsule { position, arm, radius, perp: Capsule::set_perp(arm, radius) };

    }

    fn set_perp(arm: (f32, f32), radius: f32) -> (f32, f32)
    {

        let length = f32::sqrt((arm.0 * arm.0) + (arm.1 * arm.1));
        let mut perp = (-arm.1, arm.0);

        if length > f32::EPSILON
        {

            perp = ((perp.0 * radius) / length, (perp.1 * radius) / length);

        }
        
        return perp;

    }

    /// Used to change the radius of the capsule.
    pub fn set_radius(&mut self, radius: f32)
    {

        self.radius = radius;

        self.perp = Capsule::set_perp(self.arm, self.radius);

    }

    /// Used to access the radius of the capsule.
    pub fn radius(&self) -> f32
    {

        return self.radius;

    }

    /// Used to change the arm of the capsule. Remember that the arm
    /// is half of the capsule's length, not the entire length.
    pub fn set_arm(&mut self, arm: (f32, f32))
    {

        self.arm = arm;
        
        self.perp = Capsule::set_perp(self.arm, self.radius);

    }

    ///Used to access the arm vector of the capsule.
    pub fn arm(&self) -> (f32, f32)
    {

        return self.arm;

    }

    ///Used to access the vector perpendicular to the capsule. Set
    ///by the values of radius and arm.
    pub fn perp(&self) -> (f32, f32)
    {

        return self.perp;

    }

    fn points(&self) -> [(f32, f32); 6]
    {

        //TODO: Determine if this needs to be optimized or if the compiler does it for us
        return 
        [
            
            self.arm,
            (-self.arm.0, -self.arm.1),
            (self.arm.0 + self.perp.0, self.arm.1 + self.perp.1),
            (self.arm.0 - self.perp.0, self.arm.1 - self.perp.1),
            (-self.arm.0 + self.perp.0, -self.arm.1 + self.perp.1),
            (-self.arm.0 - self.perp.0, -self.arm.1 - self.perp.1)

        ];

    }
    
}

impl crate::Shape for Capsule
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

        return 4;

    }


    fn get_axis(&self, index: usize, target: (f32, f32)) -> (f32, f32)
    {

        return match index
        {

            0 => self.arm,
            1 => self.perp,
            2 => (target.0 - self.position.0 - self.arm.0, target.1 - self.position.1 - self.arm.1),
            _ => (target.0 - self.position.0 + self.arm.0, target.1 - self.position.1 + self.arm.1)
            
        };


    }

    fn project(&self, axis: (f32, f32), normalize: bool) -> (f32, f32)
    {

        let circle1 = Circle::new((self.position.0 + self.arm.0, self.position.1 + self.arm.1), self.radius);
        let circle2 = Circle::new((self.position.0 - self.arm.0, self.position.1 - self.arm.1), self.radius);

        let projection1 = circle1.project(axis, normalize);
        let projection2 = circle2.project(axis, normalize);

        return (f32::min(projection1.0, projection2.0), f32::max(projection1.1, projection2.1));

    }

    fn needs_closest(&self, index: usize) -> bool
    {

        return match index
        {

            0 => false,
            1 => false,
            2 => true,
            _ => true
            
        };

    }

    fn get_closest(&self, target: (f32, f32)) -> (f32, f32)
    {

        return crate::closest(self.position, target, &self.points());

    }

    fn point(&self, index: usize) -> (f32, f32)
    {

        return match index
        {

            0 => (self.position.0 + self.arm.0 + self.perp.0, self.position.1 + self.arm.1 + self.perp.1),
            1 => (self.position.0 + self.arm.0 + self.perp.0, self.position.1 + self.arm.1 + self.perp.1),
            2 => (self.position.0 + self.arm.0, self.position.1 + self.arm.1),
            _ => (self.position.0 - self.arm.0, self.position.1 - self.arm.1)
            
        };

    }

}

impl crate::Rotate for Capsule
{

    fn rotate(&mut self, angle: f32)
    {

        let sin = f32::sin(angle);
        let cos = f32::cos(angle);

        self.rotate_sincos(sin, cos);

    }

    fn rotate_sincos(&mut self, sin: f32, cos: f32)
    {

        let arm = crate::rotate!(sin, cos, self.arm);
        self.set_arm(arm);

    }

}

#[cfg(feature = "serde")]
#[derive(Serialize, Deserialize)]
#[serde(rename = "Capsule")]
struct Cap
{

    position: (f32, f32),
    arm: (f32, f32),
    radius: f32

}

#[cfg(feature = "serde")]
impl Serialize for Capsule
{

    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error> where S: Serializer
    {

        let mut state = serializer.serialize_struct("Capsule", 3)?;
        state.serialize_field("position", &self.position)?;
        state.serialize_field("arm", &self.arm)?;
        state.serialize_field("radius", &self.radius)?;
        state.end()

    }

}

#[cfg(feature = "serde")]
impl <'de> Deserialize<'de> for Capsule
{

    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error> where D: Deserializer<'de>
    {

        let raw = <Cap>::deserialize(deserializer)?;
        return Ok(Capsule::new(raw.position, raw.arm, raw.radius));

    }

}

#[cfg(test)]
mod capsule_tests
{

    use super::*;
    use crate::{float_equal, Shape};

    #[test]
    fn test_num_axes()
    {

        let capsule = Capsule::new((1.0, 2.0), (2.0, 3.0), 3.0);

        assert_eq!(capsule.num_axes(), 4);

    }

    #[test]
    fn test_get_axis()
    {

        let capsule = Capsule::new((0.0, 0.0), (0.0, 2.0), 1.0);

        let axis1 = capsule.get_axis(0, (4.0, 2.0));
        let axis2 = capsule.get_axis(1, (4.0, 2.0));
        let axis3 = capsule.get_axis(2, (4.0, 2.0));
        let axis4 = capsule.get_axis(3, (4.0, 2.0));

        assert!(float_equal(axis1.0, 0.0));
        assert!(float_equal(axis1.1, 2.0));
        assert!(float_equal(axis2.0, -1.0));
        assert!(float_equal(axis2.1, 0.0));
        assert!(float_equal(axis3.0, 4.0));
        assert!(float_equal(axis3.1, 0.0));
        assert!(float_equal(axis4.0, 4.0));
        assert!(float_equal(axis4.1, 4.0));

    }

    #[test]
    fn test_project()
    {

        let capsule = Capsule::new((1.0, 2.0), (1.0, 1.0), 2.0);

        let axis1 = (0.0, 1.0);
        let axis2 = (-1.0, 2.0);

        let projection1 = capsule.project(axis1, true);
        let projection2 = capsule.project(axis2, false);

        assert!(float_equal(projection1.0, -1.0));
        assert!(float_equal(projection1.1, 5.0));
        assert!(float_equal(projection2.0, 2.0 - 2.0 * f32::sqrt(5.0)));
        assert!(float_equal(projection2.1, 4.0 + 2.0 * f32::sqrt(5.0)));

    }

    #[test]
    fn test_needs_closest()
    {

        let capsule = Capsule::new((8.0, 10.0), (-2.5, 3.0), 3.4);

        assert!(!capsule.needs_closest(0));
        assert!(!capsule.needs_closest(1));
        assert!(capsule.needs_closest(2));
        assert!(capsule.needs_closest(3));

    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_serde()
    {

        let capsule: Capsule = ron::from_str("Capsule(position: (0.0, 0.0), arm: (10.0, 0.0), radius: 5.0)").unwrap();

        assert!(float_equal(capsule.perp.0, 0.0));
        assert!(float_equal(capsule.perp.1, 5.0));

        let de = ron::to_string(&capsule).unwrap();

        assert_eq!(de, "(position:(0.0,0.0),arm:(10.0,0.0),radius:5.0)");

    }

}
