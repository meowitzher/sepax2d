pub mod polygon;

use polygon::Polygon;

/// Returns true if the given polygons overlap, and false if they do not. Does not work for
/// points or line segments.
/// 
/// This method performs a floating point comparison with Rust's built in epsilon constant, so it may
/// return the incorrect answer for polygons which are very small or very close together.
/// 
/// Requires both polygons to be convex.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::{sat_overlap, polygon::Polygon};
/// 
/// let square = Polygon::from_vertices((1.0, 1.0), vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]);
/// let triangle = Polygon::from_vertices((0.0, 0.0), vec![(2.0, 2.0), (0.0, -2.0), (-1.0, 0.0)]);
/// 
/// assert!(sat_overlap(&square, &triangle));
/// ```
pub fn sat_overlap(left: &Polygon, right: &Polygon) -> bool
{

    if !polygon_overlap(left, right, false).0
    {

        return false;

    }

    if !polygon_overlap(right, left, false).0
    {

        return false;

    }

    return true;

}

/// Returns the vector that needs to be added to right's position to resolve a collision with left.
/// Does not work for points or line segments.
/// 
/// If the polygons are not colliding, it returns the zero vector.
/// 
/// This method performs a floating point comparison with Rust's built in epsilon constant, so it may
/// return the incorrect answer for polygons which are very small or very close together.
/// 
/// Requires both polygons to be convex.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::{sat_collision, polygon::Polygon};
/// 
/// let square = Polygon::from_vertices((1.0, 1.0), vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]);
/// let triangle = Polygon::from_vertices((-3.5, 1.0), vec![(4.0, 0.0), (0.0, 6.0), (-4.0, 0.0)]);
/// 
/// let resolution = sat_collision(&square, &triangle);
/// //resolution = (-0.5, 0.0) up to floating point error
/// 
/// assert!(resolution.0 + 0.5 < f64::EPSILON && resolution.0 + 0.5 > -f64::EPSILON);
/// assert!(resolution.1 < f64::EPSILON && resolution.1 > -f64::EPSILON);
/// ```
pub fn sat_collision(left: &Polygon, right: &Polygon) -> (f64, f64)
{

    let l_overlap = polygon_overlap(left, right, true);
    let r_overlap = polygon_overlap(right, left, true);

    //Ensure that the vector points from left to right
    let r_flipped = (true, r_overlap.1, (-r_overlap.2.0, -r_overlap.2.1));

    let overlap = if l_overlap.1 < r_flipped.1 { l_overlap } else { r_flipped };

    return (overlap.1 * overlap.2.0, overlap.1 * overlap.2.1);

}

fn polygon_overlap(axes: &Polygon, projected: &Polygon, normalize: bool) -> (bool, f64, (f64, f64))
{

    let mut min_overlap = f64::MAX;
    let mut min_axis = (0.0, 0.0);

    //We do not check line or point collisions
    if axes.vertices.len() > 2 && projected.vertices.len() > 2
    {

        if let Some((first_x, first_y)) = axes.vertices.first()
        {

            let mut previous = (*first_x, *first_y);

            for (x, y) in axes.vertices.iter().cycle().skip(1).take(axes.vertices.len())
            {

                let side = (*x - previous.0, *y - previous.1);
                let mut axis = (-side.1, side.0);

                //If we are just checking for overlap, we can skip normalizing the axis
                if normalize
                {

                    let length = f64::sqrt((axis.0 * axis.0) + (axis.1 * axis.1));
                    
                    if length > f64::EPSILON
                    {

                        axis = (axis.0 / length, axis.1 / length);

                    }

                }

                let (min_l, max_l) = axes.project(axis);
                let (min_r, max_r) = projected.project(axis);

                //If there is no overlap, we can return early
                if min_l > max_r - f64::EPSILON || min_r > max_l - f64::EPSILON
                {

                    return (false, 0.0, (0.0, 0.0)); 

                }

                let overlap = f64::min(max_l, max_r) - f64::max(min_l, min_r);
                if overlap < min_overlap
                {

                    min_overlap = overlap;
                    min_axis = axis;

                }

                previous = (*x, *y);

            }

        }

    }

    //Ensure that the chosen axis of penetration points from axes to projected
    let difference = (projected.position.0 - axes.position.0, projected.position.1 - axes.position.1);
    if (difference.0 * min_axis.0 + difference.1 * min_axis.1) < -f64::EPSILON
    {

        min_axis.0 *= -1.0;
        min_axis.1 *= -1.0;

    }

    return (true, min_overlap, min_axis);

}

#[cfg(test)]
mod sat_tests
{

    use super::*;

    fn float_equal(left: f64, right: f64) -> bool
    {

        return (left - right).abs() < f64::EPSILON;

    }

    #[test]
    fn test_sat_overlap()
    {

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

    }

    #[test]
    fn test_sat_no_overlap()
    {

        let triangle = Polygon::from_vertices((2.0, 2.0), vec![(-1.0, 1.0), (0.0, -1.0), (1.0, 1.0)]);
        let pentagon = Polygon::from_vertices((-3.0, 0.0), vec![(2.0, 0.0), (4.0, 1.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]);
        let triangle2 = Polygon::from_vertices((4.0, 3.0), vec![(-2.0, 1.0), (-1.0, -2.0), (2.0, 0.0)]);

        assert!(!sat_overlap(&pentagon, &triangle));
        assert!(!sat_overlap(&triangle, &pentagon));
        assert!(!sat_overlap(&pentagon, &triangle2));
        assert!(!sat_overlap(&triangle2, &pentagon));

    }

    #[test]
    fn test_sat_collision()
    {

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

    }

}