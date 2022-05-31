/// A struct representing a circle via a position and radius.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::circle::Circle;
/// use sepax2d::polygon::Polygon;
/// use sepax2d::{sat_overlap, sat_collision};
/// 
/// let circle1 = Circle::new((0.0, 0.0), 2.0);
/// let circle2 = Circle::new((2.0, 2.0), 2.0);
/// 
/// let resolution = sat_collision(&circle1, &circle2);
/// let difference = 2.0 * (f64::sqrt(2.0) - 1.0); // 2 root 2 - 2
/// assert!(resolution.0 - difference < f64::EPSILON && resolution.0 - difference > -f64::EPSILON);
/// assert!(resolution.1 - difference < f64::EPSILON && resolution.1 - difference > -f64::EPSILON);
/// 
/// let polygon = Polygon::from_vertices((0.0, 3.0), vec![(0.0, 0.0), (1.0, 0.0), (1.0, -1.5), (0.0, -1.5)]);
/// assert!(sat_overlap(&polygon, &circle1));
/// assert!(sat_overlap(&circle1, &polygon));
/// ```
pub struct Circle
{

    pub position: (f64, f64),
    pub radius: f64

}

impl Circle
{

    /// Create a new circle with given position and radius.
    pub fn new(position: (f64, f64), radius: f64) -> Circle
    {

        return Circle { position, radius };

    }

}

impl super::Shape for Circle
{

    fn position(&self) -> (f64, f64)
    {

        return self.position;

    }

    fn num_axes(&self) -> usize
    {

        return 1;

    }


    fn get_axis(&self, _index: usize, target: (f64, f64)) -> (f64, f64)
    {

        return (target.0 - self.position.0, target.1 - self.position.1);

    }

    fn project(&self, axis: (f64, f64), normalize: bool) -> (f64, f64)
    {

        let projection = (self.position.0 * axis.0) + (self.position.1 * axis.1);
        
        //The projection of the circle along the axis is found by simply adding distance r
        //along it. However, if the vector is not a unit vector, then adding r equates to
        //adding r times the length.
        let mut magnitude = 1.0;
        if !normalize
        {

            magnitude = f64::sqrt((axis.0 * axis.0) + (axis.1 * axis.1));

        }

        return (projection - (self.radius * magnitude), projection + (self.radius * magnitude));

    }

    fn needs_closest(&self) -> bool
    {

        return true;

    }

    fn get_closest(&self, _target: (f64, f64)) -> (f64, f64)
    {

        return self.position;

    }

}

#[cfg(tests)]
mod circle_tests
{
}