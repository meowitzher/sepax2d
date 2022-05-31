/// A struct representing a circle via a position and radius.
/// 
/// # Examples
/// 
/// ```
/// 
/// ```
pub struct Circle
{

    pub position: (f64, f64),
    pub radius: f64

}

impl Circle
{

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