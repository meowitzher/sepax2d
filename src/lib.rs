//! [![Crates.io](https://img.shields.io/crates/v/sepax2d.svg)](https://crates.io/crates/sepax2d)
//! [![MIT/Apache 2.0](https://img.shields.io/badge/license-MIT%2FApache-blue.svg)](./LICENSE)
//!
//! # sepax2d
//! A safe Rust crate for finding and resolving collisions of convex shapes using the Separating Axis Theorem (SAT) in two dimensions.
//!
//! ### Usage
//!
//! Add the following to the `[dependencies]` section of your `Cargo.toml`:
//!
//! ```toml
//! sepax2d = "0.3"
//! ```
//!
//! Import the types of shapes you'd like to use and create some new shapes:
//!
//! ```rust
//! use sepax2d::prelude::*;
//! #
//! # let top_left = (1.0, 1.0);
//! # let width = 5.0;
//! # let height = 5.0;
//! # let center = (2.0, 0.0);
//! # let radius = 2.0;
//! # let arm = (1.0, 0.0);
//! # let position = (-1.0, -1.0);
//! # let u_vector = (2.0, 1.0);
//! # let v_vector = (0.0, -2.0);
//!
//! let rectangle = AABB::new(top_left, width, height);
//! let circle = Circle::new(center, radius);
//! let capsule = Capsule::new(center, arm, radius);
//! let parallelogram = Parallelogram::new(position, u_vector, v_vector);
//!
//! //The vertices of a polygon are position vectors
//! //relative to the shape's position, i.e. if position
//! //is (1, 2), then the absolute location of the
//! //first vertex is (1, 0).
//! let triangle = Polygon::from_vertices
//! (
//!
//!     position,
//!     vec![(0.0, -2.0), (-1.0, 2.0), (1.0, 2.0)]
//!
//! );
//! ```
//!
//! Use the `sat_overlap(&left, &right)` method to get a `bool` indicating whether or not the two given shapes overlap.
//! Any struct implementing the `Shape` trait can be used.
//!
//! ```rust
//! # use sepax2d::prelude::*;
//! #
//! # let top_left = (1.0, 1.0);
//! # let width = 5.0;
//! # let height = 5.0;
//! # let center = (2.0, 0.0);
//! # let radius = 2.0;
//! # let arm = (1.0, 0.0);
//! # let position = (-1.0, -1.0);
//! # let rectangle = AABB::new(top_left, width, height);
//! # let circle = Circle::new(center, radius);
//! # let capsule = Capsule::new(center, arm, radius);
//! #
//! # let triangle = Polygon::from_vertices
//! # (
//! #
//! #    position,
//! #    vec![(0.0, -2.0), (-1.0, 2.0), (1.0, 2.0)]
//! #
//! # );
//! #
//! assert!(sat_overlap(&circle, &capsule));
//! assert!(!sat_overlap(&triangle, &rectangle));
//! ```
//!
//! Use the `sat_collision(&left, &right)` method to get a `(f32, f32)` which represents the shift needed to add to `right`'s
//! position in order to resolve the overlap of the two shapes.
//!
//! ```rust
//! # use sepax2d::prelude::*;
//! #
//! # let top_left = (1.0, 1.0);
//! # let width = 5.0;
//! # let height = 5.0;
//! # let center = (2.0, 0.0);
//! # let radius = 2.0;
//! # let arm = (1.0, 0.0);
//! # let position = (-1.0, -1.0);
//! #
//! # let rectangle = AABB::new(top_left, width, height);
//! # let circle = Circle::new(center, radius);
//! # let capsule = Capsule::new(center, arm, radius);
//! #
//! # let mut triangle = Polygon::from_vertices
//! # (
//! #
//! #    position,
//! #    vec![(0.0, -2.0), (-1.0, 2.0), (1.0, 2.0)]
//! #
//! # );
//! #
//! let resolution = sat_collision(&circle, &triangle);
//!
//! let position = triangle.position();
//! triangle.set_position((position.0 + resolution.0, position.1 + resolution.1));
//!
//! assert!(!sat_overlap(&circle, &triangle));
//! ```
//!
//! Use the `contains_point(&shape, point)` method to get a `bool` indicating whether or not the specified point
//! is inside the given shape.
//!
//! ```rust
//! # use sepax2d::prelude::*;
//!
//! let rect = AABB::new((0.0, 0.0), 5.0, 2.0);
//!
//! assert!(contains_point(&rect, (1.0, 1.0)));
//! assert!(!contains_point(&rect, (10.0, -1.0)));
//! ```
//!
//! `Polygon`, `Circle`, `Capsule`, and `Parallelogram` shapes implement the `Rotate` trait, which allows you to rotate them
//! around their `position`.
//!
//! ```rust
//! # use sepax2d::prelude::*;
//! # let position = (-1.0, -1.0);
//!
//! let mut triangle = Polygon::from_vertices
//! (
//!
//!    position,
//!    vec![(-1.0, 0.0), (0.0, 2.0), (1.0, 0.0)]
//!
//! );
//!
//! triangle.rotate(std::f32::consts::FRAC_PI_2)
//! //New vertices: [(0.0, -1.0), (-2.0, 0.0), (0.0, 1.0)]
//! ```
//!
//! You can use the `intersects_line`, `intersects_ray`, and `intersects_segment` methods to
//! check whether a shape intersects with the corresponding type of line.
//!
//! ```rust
//! # use sepax2d::prelude::*;
//!
//! let triangle = Polygon::from_vertices((0.0, vec![(0.0, 0.0), (1.0, 1.0), (-1.0, 1.0)]));
//!
//! assert!(intersects_segment(&triangle, (2.0, 0.5), (-2.0, 0.5)));
//! ```
//!
//! ### Features
//!
//! Enable the `serde` feature for (De)Serialization of supported shapes!

#![allow(clippy::needless_return)]

pub mod polygon;
pub mod circle;
pub mod aabb;
pub mod capsule;
pub mod parallelogram;

pub mod line;

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

/// A trait indicating that a shape can be rotated around its position. Applicable
/// to all shapes other than AABB.
pub trait Rotate
{

    /// Rotate the shape by the given angle, with the rotation counterclockwise when
    /// the Y-axis points up.
    fn rotate(&mut self, angle: f32);

    /// Rotate the shape using the given sine and cosine of an angle. Use this when
    /// you are rotating multiple shapes by the same angle and don't want to re-calculate
    /// the trig functions.
    fn rotate_sincos(&mut self, sin: f32, cos: f32);

}

//Helper macro to rotate the given 2D vector v by the rotation matrix with sine s and cosine c
#[macro_export]
macro_rules! rotate
{

    ($s: expr, $c: expr, $v: expr) =>
    {

        ($c * $v.0 - $s * $v.1, $s * $v.0 + $c * $v.1   )

    };

}

/// Returns true if the given shapes overlap, and false if they do not. Does not work for
/// degenerate shapes.
///
/// This method performs a floating point comparison with Rust's built in epsilon constant, so it may
/// return the incorrect answer for shapes which are very small or very close together.
///
/// Requires both shapes to be convex.
///
/// # Examples
///
/// ```
/// use sepax2d::prelude::*;
///
/// let square = Polygon::from_vertices((1.0, 1.0), vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]);
/// let triangle = Polygon::from_vertices((0.0, 0.0), vec![(2.0, 2.0), (0.0, -2.0), (-1.0, 0.0)]);
///
/// assert!(sat_overlap(&square, &triangle));
/// ```
pub fn sat_overlap(left: &(impl Shape + ?Sized), right: &(impl Shape + ?Sized)) -> bool
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
/// the first shape. Does not work for degenerate shapes.
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
/// use sepax2d::prelude::*;
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
pub fn sat_collision(left: &(impl Shape + ?Sized), right: &(impl Shape + ?Sized)) -> (f32, f32)
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
/// use sepax2d::prelude::*;
///
/// let square = Polygon::from_vertices((1.0, 1.0), vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]);
/// let triangle = Polygon::from_vertices((0.0, 0.0), vec![(2.0, 2.0), (0.0, -2.0), (-1.0, 0.0)]);
///
/// assert!(contains_point(&triangle, (0.5, 0.5)));
/// assert!(!contains_point(&square, (-2.0, 2.0)));
/// ```
pub fn contains_point(shape: &(impl Shape + ?Sized), point: (f32, f32)) -> bool
{

    let polygon = polygon::Polygon::from_vertices(point, vec![(0.0, 0.0)]);

    return shape_overlap(shape, &polygon, false).0;

}

fn shape_overlap(axes: &(impl Shape + ?Sized), projected: &(impl Shape + ?Sized), normalize: bool) -> (bool, f32, (f32, f32))
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

        let overlap = f32::min(max_l - min_r, max_r - min_l);
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

    return (left - right).abs() < 0.00001;

}

#[cfg(test)]
mod sat_tests
{

    use super::*;
    use super::prelude::*;

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
    fn test_parallelogram_collision()
    {

        let gram = Parallelogram::new((0.0, -0.5), (2.0, 1.0), (-1.0, -1.0));

        let triangle = Polygon::from_vertices((0.0, 5.0), vec![(0.0, -2.0), (-1.0, 2.0), (1.0, 2.0)]);
        let rectangle = AABB::new((-1.0, 0.0), 4.0, 2.0);
        let circle = Circle::new((2.0, -2.5), 1.0);

        assert!(!sat_overlap(&gram, &triangle));
        assert!(!sat_overlap(&gram, &circle));

        assert!(sat_overlap(&gram, &rectangle));

        let resolution = sat_collision(&rectangle, &gram);
        println!("{:?}", resolution);
        
        assert!(float_equal(resolution.0, 0.0));
        assert!(float_equal(resolution.1, -0.5));

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

    #[test]
    fn test_rotate()
    {

        let mut triangle = Polygon::from_vertices((0.0, 0.0), vec![(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)]);
        let mut capsule = Capsule::new((2.0, 1.0), (2.0, 0.0), 4.0);
        let mut gram = Parallelogram::new((3.0, 4.0), (2.0, 1.0), (-1.0, 1.0));

        capsule.rotate(std::f32::consts::FRAC_PI_4);

        assert!(float_equal(capsule.arm().0, 2.0 / f32::sqrt(2.0)));
        assert!(float_equal(capsule.arm().1, 2.0 / f32::sqrt(2.0)));
        assert!(float_equal(capsule.perp().0, -4.0 / f32::sqrt(2.0)));
        assert!(float_equal(capsule.perp().1, 4.0 / f32::sqrt(2.0)));

        let sin = f32::sin(std::f32::consts::PI);
        let cos = f32::cos(std::f32::consts::PI);

        triangle.rotate_sincos(sin, cos);

        assert!(float_equal(triangle.vertices[1].0, -2.0));
        assert!(float_equal(triangle.vertices[1].1, 0.0));
        assert!(float_equal(triangle.vertices[2].0, 0.0));
        assert!(float_equal(triangle.vertices[2].1, -1.0));

        gram.rotate_sincos(sin, cos);

        assert!(float_equal(gram.u.0, -2.0));
        assert!(float_equal(gram.u.1, -1.0));
        assert!(float_equal(gram.v.0, 1.0));
        assert!(float_equal(gram.v.1, -1.0));

    }

}

pub mod prelude
{

    pub use crate::{sat_overlap, sat_collision, contains_point};

    pub use crate::Shape;
    pub use crate::Rotate;

    pub use crate::polygon::Polygon;
    pub use crate::circle::Circle;
    pub use crate::aabb::AABB;
    pub use crate::capsule::Capsule;
    pub use crate::parallelogram::Parallelogram;

    pub use crate::line::intersects_line;
    pub use crate::line::intersects_ray;
    pub use crate::line::intersects_segment;

}