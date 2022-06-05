#![allow(clippy::needless_return)]

pub mod polygon;
pub mod circle;
pub mod aabb;
pub mod capsule;

/// A trait describing the behavior needed to implement SAT overlap and collision
/// for a given shape.
pub trait Shape
{

    /// The location of the shape in 2D space. 
    fn position(&self) -> (f32, f32);

    /// Set the location of the shape.
    fn set_position(&mut self, position: (f32, f32));

    /// The number of axes the shape provides for testing. For polygons, it is
    /// the same as the number of vertices, but for circles it is simply one. 
    fn num_axes(&self) -> usize;

    /// The method used to access the axes during the SAT calculations. This is
    /// used to avoid the memory allocation of a new vector or array each time 
    /// we calculate collisions.
    fn get_axis(&self, index: usize, target: (f32, f32)) -> (f32, f32);

    /// Getting the minimum and maximum projection of the shape onto the given axis
    /// to look for overlap. Normalize denotes whether or not the axis passed in is
    /// a unit vector to avoid repeating calculations.
    fn project(&self, axis: (f32, f32), normalize: bool) -> (f32, f32);

    /// Determine whether or not the shape needs access to the closest vertex of 
    /// another shape to check collisions.
    fn needs_closest(&self, index: usize) -> bool;

    /// Gets the closest vertex/primary point/position to the given target, NOT the closest point
    /// on the shape.
    fn get_closest(&self, target: (f32, f32)) -> (f32, f32);

    /// The point corresponding to the given axis, if applicable. Otherwise, position.
    fn point(&self, index: usize) -> (f32, f32);

}

/// Returns true if the given shapes overlap, and false if they do not. Does not work for
/// degenerate polygons.
/// 
/// This method performs a floating point comparison with Rust's built in epsilon constant, so it may
/// return the incorrect answer for shapes which are very small or very close together.
/// 
/// Requires both shapes to be convex.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::sat_overlap;
/// use sepax2d::polygon::Polygon;
/// 
/// let square = Polygon::from_vertices((1.0, 1.0), vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]);
/// let triangle = Polygon::from_vertices((0.0, 0.0), vec![(2.0, 2.0), (0.0, -2.0), (-1.0, 0.0)]);
/// 
/// assert!(sat_overlap(&square, &triangle));
/// ```
pub fn sat_overlap(left: &impl Shape, right: &impl Shape) -> bool
{

    if !shape_overlap(left, right, false).0
    {

        return false;

    }

    if !shape_overlap(right, left, false).0
    {

        return false;

    }

    return true;

}

/// Returns the vector that needs to be added to the second shape's position to resolve a collision with 
/// the first shape. Does not work for degenerate polygons.
/// 
/// If the shapes are not colliding, it returns the zero vector.
/// 
/// This method performs a floating point comparison with Rust's built in epsilon constant, so it may
/// return the incorrect answer for shapes which are very small or very close together.
/// 
/// Requires both shapes to be convex.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::sat_collision;
/// use sepax2d::polygon::Polygon;
/// use sepax2d::aabb::AABB;
/// 
/// let square = Polygon::from_vertices((1.0, 1.0), vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]);
/// let triangle = Polygon::from_vertices((-3.5, 1.0), vec![(4.0, 0.0), (0.0, 6.0), (-4.0, 0.0)]);
/// 
/// let aabb = AABB::new((0.0, 2.0), 2.0, 2.0);
/// 
/// let resolution = sat_collision(&square, &triangle);
/// //resolution = (-0.5, 0.0) up to floating point error
/// 
/// assert!(resolution.0 + 0.5 < f32::EPSILON && resolution.0 + 0.5 > -f32::EPSILON);
/// assert!(resolution.1 < f32::EPSILON && resolution.1 > -f32::EPSILON);
/// 
/// let aabb_resolution = sat_collision(&aabb, &triangle);
/// //resolution = (-0.5, 0.0) up to floating point error. aabb is a different way to
/// //represent the same shape as square
/// 
/// assert!(resolution.0 + 0.5 < f32::EPSILON && resolution.0 + 0.5 > -f32::EPSILON);
/// assert!(resolution.1 < f32::EPSILON && resolution.1 > -f32::EPSILON);
/// 
/// ```
pub fn sat_collision(left: &impl Shape, right: &impl Shape) -> (f32, f32)
{

    let l_overlap = shape_overlap(left, right, true);
    let r_overlap = shape_overlap(right, left, true);

    //Ensure that the vector points from left to right
    let r_flipped = (true, r_overlap.1, (-r_overlap.2.0, -r_overlap.2.1));

    let overlap = if l_overlap.1 < r_flipped.1 { l_overlap } else { r_flipped };

    return (overlap.1 * overlap.2.0, overlap.1 * overlap.2.1);

}

/// Returns true if the given shape contains the specified point, and false if 
/// it does not. Does not work for degenerate polygons.
/// 
/// This method performs a floating point comparison with Rust's built in epsilon constant, so it may
/// return the incorrect answer for shapes which are very small or very close to the point.
/// 
/// Requires the shape to be convex.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::contains_point;
/// use sepax2d::polygon::Polygon;
/// 
/// let square = Polygon::from_vertices((1.0, 1.0), vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]);
/// let triangle = Polygon::from_vertices((0.0, 0.0), vec![(2.0, 2.0), (0.0, -2.0), (-1.0, 0.0)]);
/// 
/// assert!(contains_point(&triangle, (0.5, 0.5)));
/// assert!(!contains_point(&square, (-2.0, 2.0)));
/// ```
pub fn contains_point(shape: &impl Shape, point: (f32, f32)) -> bool
{

    let polygon = polygon::Polygon::from_vertices(point, vec![(0.0, 0.0)]);

    return shape_overlap(shape, &polygon, false).0;

}

fn shape_overlap(axes: &impl Shape, projected: &impl Shape, normalize: bool) -> (bool, f32, (f32, f32))
{

    let mut min_overlap = f32::MAX;
    let mut min_axis = (0.0, 0.0);

    let num_axes = axes.num_axes();
    for i in 0..num_axes
    {

        let closest = if axes.needs_closest(i) { projected.get_closest(axes.point(i)) } else { (0.0, 0.0) };
        let mut axis = axes.get_axis(i, closest);

        //If we are just checking for overlap, we can skip normalizing the axis. However,
        //we need to normalize to find the minimum penetration vector.
        if normalize
        {

            let length = f32::sqrt((axis.0 * axis.0) + (axis.1 * axis.1));
                    
            if length > f32::EPSILON
            {

                axis = (axis.0 / length, axis.1 / length);

            }

        }

        let (min_l, max_l) = axes.project(axis, normalize);
        let (min_r, max_r) = projected.project(axis, normalize);

        //If there is no overlap, we can return early
        if min_l > max_r - f32::EPSILON || min_r > max_l - f32::EPSILON
        {

            return (false, 0.0, (0.0, 0.0)); 

        }

        let overlap = f32::min(max_l, max_r) - f32::max(min_l, min_r);
        if overlap < min_overlap
        {

            min_overlap = overlap;
            min_axis = axis;

        }

    }

    //Ensure that the chosen axis of penetration points from axes to projected
    let axes_position = axes.position();
    let projected_position = projected.position();
    let difference = (projected_position.0 - axes_position.0, projected_position.1 - axes_position.1);
    if (difference.0 * min_axis.0 + difference.1 * min_axis.1) < -f32::EPSILON
    {

        min_axis.0 *= -1.0;
        min_axis.1 *= -1.0;

    }

    return (true, min_overlap, min_axis);

}

fn project(position: (f32, f32), axis: (f32, f32), points: &[(f32, f32)]) -> (f32, f32)
{

    let mut min = f32::MAX;
    let mut max = f32::MIN;

    for (x, y) in points
    {

        let position = (position.0 + *x, position.1 + *y);

        let projection = (position.0 * axis.0) + (position.1 * axis.1);

        min = f32::min(min, projection);
        max = f32::max(max, projection);

    }

    return (min, max);

}

fn closest(position: (f32, f32), target: (f32, f32), points: &[(f32, f32)]) -> (f32, f32)
{

    let mut point = (0.0, 0.0);
    let mut min = f32::MAX;

    for (x, y) in points
    {

        let position = (position.0 + *x, position.1 + *y);
        let dist_square = (position.0 - target.0) * (position.0 - target.0) + (position.1 - target.1) * (position.1 - target.1);

        if dist_square < min
        {

            point = position;
            min = dist_square;

        }

    }

    return point;

}

#[allow(dead_code)]
fn float_equal(left: f32, right: f32) -> bool
{

    return (left - right).abs() < f32::EPSILON;

}

#[cfg(test)]
mod sat_tests
{

    use super::*;
    use super::polygon::Polygon;
    use super::circle::Circle;
    use super::aabb::AABB;
    use super::capsule::Capsule;

    #[test]
    fn test_sat_overlap()
    {

        //Polygons
        let square = Polygon::from_vertices((0.0, 0.0), vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]);
        let triangle = Polygon::from_vertices((2.0, 2.0), vec![(-1.0, 1.0), (0.0, -1.0), (1.0, 1.0)]);
        let pentagon = Polygon::from_vertices((-3.0, 0.0), vec![(2.0, 0.0), (4.0, 1.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]);
        let triangle2 = Polygon::from_vertices((4.0, 3.0), vec![(-2.0, 1.0), (-1.0, -2.0), (2.0, 0.0)]);

        assert!(sat_overlap(&pentagon, &pentagon));
        assert!(sat_overlap(&square, &triangle));
        assert!(sat_overlap(&triangle, &square));
        assert!(sat_overlap(&square, &pentagon));
        assert!(sat_overlap(&pentagon, &square));
        assert!(sat_overlap(&square, &triangle2));
        assert!(sat_overlap(&triangle2, &square));
        assert!(sat_overlap(&triangle, &triangle2));
        assert!(sat_overlap(&triangle, &triangle2));

        //Circles
        let circle1 = Circle::new((-1.0, -1.0), 2.0);
        let circle2 = Circle::new((-2.0, -1.0), 1.1);

        assert!(sat_overlap(&circle1, &square));
        assert!(sat_overlap(&circle2, &pentagon));
        assert!(sat_overlap(&pentagon, &circle1));
        assert!(sat_overlap(&circle1, &circle2));

    }

    #[test]
    fn test_sat_no_overlap()
    {

        //Polygons
        let triangle = Polygon::from_vertices((2.0, 2.0), vec![(-1.0, 1.0), (0.0, -1.0), (1.0, 1.0)]);
        let pentagon = Polygon::from_vertices((-3.0, 0.0), vec![(2.0, 0.0), (4.0, 1.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]);
        let triangle2 = Polygon::from_vertices((4.0, 3.0), vec![(-2.0, 1.0), (-1.0, -2.0), (2.0, 0.0)]);

        assert!(!sat_overlap(&pentagon, &triangle));
        assert!(!sat_overlap(&triangle, &pentagon));
        assert!(!sat_overlap(&pentagon, &triangle2));
        assert!(!sat_overlap(&triangle2, &pentagon));

        //Circles
        let circle1 = Circle::new((-1.0, -1.0), 2.0);
        let circle2 = Circle::new((-2.0, -1.0), 1.1);
        let circle3 = Circle::new((2.0, 5.0), 1.0);

        assert!(!sat_overlap(&triangle, &circle1));
        assert!(!sat_overlap(&circle2, &triangle2));
        assert!(!sat_overlap(&circle1, &circle3));
        assert!(!sat_overlap(&circle3, &circle2));

    }

    #[test]
    fn test_sat_collision()
    {

        //Polygons
        let square = Polygon::from_vertices((0.0, 0.0), vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]);
        let rectangle = Polygon::from_vertices((1.0, -2.0), vec![(0.0, 0.0), (1.0, 0.0), (1.0, 2.1), (0.0, 2.1)]);
        let triangle = Polygon::from_vertices((0.0, 0.0), vec![(0.5, 0.6), (0.0, -1.0), (-0.5, 0.6)]);
        let triangle2 = Polygon::from_vertices((-0.4, -0.4), vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);

        let resolution = sat_collision(&square, &rectangle);
        assert!(float_equal(resolution.0, 0.0));
        assert!(float_equal(resolution.1, -0.1));

        let resolution1 = sat_collision(&rectangle, &square);
        assert!(float_equal(resolution1.0, 0.0));
        assert!(float_equal(resolution1.1, 0.1));

        let resolution2 = sat_collision(&square, &triangle);
        assert!(float_equal(resolution2.0, -0.5));
        assert!(float_equal(resolution2.1, 0.0));

        let resolution3 = sat_collision(&square, &triangle2);
        assert!(float_equal((resolution3.0 * -1.0) + (resolution3.1), 0.0)); //Assert that the minimum axis of penetration is (1,1)

        //Circles
        let circle = Circle::new((2.0, -0.5), 1.0);
        let circle2 = Circle::new((0.0, -0.5), 1.1);

        let resolution4 = sat_collision(&circle, &square);
        assert!(float_equal(resolution4.0, 0.0));
        assert!(float_equal(resolution4.1, 0.5));

        let resolution5 = sat_collision(&circle, &triangle2);
        assert!(float_equal(resolution5.0, 0.0));
        assert!(float_equal(resolution5.1, 0.0));

        let resolution6 = sat_collision(&circle2, &circle);
        assert!(float_equal(resolution6.0, 0.1));
        assert!(float_equal(resolution6.1, 0.0));

    }

    #[test]
    fn test_capsule_collision()
    {

        let capsule = Capsule::new((0.0, 0.0), (0.0, 2.0), 2.0);

        let triangle = Polygon::from_vertices((0.0, 5.0), vec![(0.0, -2.0), (-1.0, 2.0), (1.0, 2.0)]);
        let rectangle = AABB::new((-4.0, 4.0), 2.5, 0.5);
        let circle1 = Circle::new((2.0, -2.5), 1.0);
        let circle2 = Circle::new((3.0, 0.0), 1.5);

        assert!(sat_overlap(&capsule, &circle1));
        assert!(!sat_overlap(&rectangle, &capsule));

        let resolution1 = sat_collision(&circle2, &capsule);
        let resolution2 = sat_collision(&capsule, &triangle);

        assert!(float_equal(resolution1.0, -0.5));
        assert!(float_equal(resolution1.1, 0.0));
        assert!(float_equal(resolution2.0, 0.0));
        assert!(float_equal(resolution2.1, 1.0));

    }

    #[test]
    fn test_contains_point()
    {

        let capsule = Capsule::new((0.0, 0.0), (0.0, 2.0), 2.0);
        let triangle = Polygon::from_vertices((0.0, 5.0), vec![(0.0, -2.0), (-1.0, 2.0), (1.0, 2.0)]);
        let rectangle = AABB::new((-4.0, 4.0), 2.5, 0.5);
        let circle = Circle::new((2.0, -2.5), 1.0);

        assert!(contains_point(&capsule, (1.0, 1.0)));
        assert!(contains_point(&triangle, (0.0, 5.0)));
        assert!(contains_point(&rectangle, (-2.0, 4.1)));
        assert!(contains_point(&circle, (2.5, -3.0)));


        assert!(!contains_point(&capsule, (2.0, 4.0)));
        assert!(!contains_point(&triangle, (0.0, 0.0)));
        assert!(!contains_point(&rectangle, (-4.0, 3.9)));
        assert!(!contains_point(&circle, (1.5, -3.5)));

    }

}