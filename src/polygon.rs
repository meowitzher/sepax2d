#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize};

/// A polygon with a position and finitely many vertices given in either clockwise or
/// counterclockwise orientation.
/// 
/// The polygon is represented as a star shape: to find the absolute position of the 
/// n-th vertex, add the `position` to the n-th element of `vertices`.
/// 
/// The position is NOT used as a vertex by default: if you want your position to
/// be a vertex, then define a vertex with offset `(0.0, 0.0)`.
/// 
/// The order of the vertices container determines the order in which the vertices are connected.
/// Both clockwise and counterclockwise orientations are valid.
/// 
/// # Examples
/// 
/// ```
/// use sepax2d::prelude::*;
/// 
/// let mut triangle = Polygon::new((3.0, 2.0));
/// triangle.add((-1.0, -1.0));
/// triangle.add((0.0, 1.0));
/// triangle.add((1.0, -1.0)); 
/// //Triangle with vertices (2, 1), (3, 3), and (4, 1)
/// 
/// assert!(triangle.is_convex());
/// 
/// let vertices = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)];
/// let rectangle = Polygon::from_vertices((0.0, 0.0), vertices); 
/// //Rectangle with vertices (0,0), (2,0), (2,1), and (0, 1)
/// 
/// assert!(rectangle.is_convex());
/// 
/// let concave_vertices = vec![(0.0, 0.0), (2.0, 0.0), (0.0, 1.0), (2.0, 1.0)];
/// let concave_shape = Polygon::from_vertices((0.0, 0.0), concave_vertices); 
/// //Non-convex hour-glass shape with the same vertices
/// 
/// assert!(!concave_shape.is_convex());
/// ```
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polygon
{

    pub position: (f32, f32),
    pub vertices: Vec::<(f32,f32)>

}

impl Polygon
{

    /// Define a new polygon with no vertices at the given position. 
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::prelude::*;
    /// 
    /// let point = Polygon::new((17.0, 5.0));
    /// ```
    pub fn new(position: (f32, f32)) -> Polygon
    {

        return Polygon { position, vertices: Vec::new() };

    }

    /// Returns `true` when the polygon is convex, i.e. when none of its points lie on opposite sides of one
    /// of its perimeter segements, and `false` when it is concave.
    /// 
    /// Returns true automatically when the shape has two or fewer vertices, as points and line segments are
    /// trivially convex.
    /// 
    /// This method performs a floating point comparison with Rust's built in epsilon constant, so it may
    /// return the incorrect answer for polygons which are almost complex or almost concave.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::prelude::*;
    /// 
    /// let vertices = vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)];
    /// let square = Polygon::from_vertices((0.0, 0.0), vertices);
    /// 
    /// assert!(square.is_convex());
    /// 
    /// let concave_vertices = vec![(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 1.0), (0.0, 2.0)];
    /// let concave_shape = Polygon::from_vertices((0.0, 0.0), concave_vertices);
    /// 
    /// assert!(!concave_shape.is_convex());
    /// ```
    pub fn is_convex(&self) -> bool
    {

        if self.vertices.len() > 2
        {

            if let Some((first_x, first_y)) = self.vertices.first()
            {

                //We do not need to use position, as convexity is independent of position
                let mut previous = (*first_x, *first_y);

                //For each side, project the polygon onto its perpendicular axis
                for (i, (x, y)) in self.vertices.iter().cycle().skip(1).take(self.vertices.len()).enumerate()
                {

                    let side = (*x - previous.0, *y - previous.1);
                    let axis = (-side.1, side.0);

                    let (min, max) = self.side_projection(axis, i + 1);

                    if min < -f32::EPSILON && max > f32::EPSILON
                    {

                        //There are points on both sides of a polygon edge, which must mean it is not convex
                        return false;

                    }

                    previous = (*x, *y);

                }

            }

        }

        //If the polygon doesn't have enough vertices, then it is trivially convex
        return true;

    }

    fn side_projection(&self, axis: (f32, f32), start: usize) -> (f32, f32)
    {

        let mut min = f32::MAX;
        let mut max = f32::MIN;

        //Project all vertices not touching the given side onto the given axis
        for (x, y) in self.vertices.iter().cycle().skip(start + 1).take(self.vertices.len() - 2)
        {

            //Dot product the current position vector to the axis of projection
            let projection = (*x * axis.0) + (*y * axis.1);

            min = f32::min(min, projection);
            max = f32::max(max, projection);

        }

        return (min, max);

    }

    /// Adds a vertex to the given shape, setting it as the last vertex.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::prelude::*;
    /// 
    /// let vertices = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)];
    /// let mut polygon = Polygon::from_vertices((0.0, 0.0), vertices); //Triangle with vertices (0, 0), (2, 0), and (2, 2)
    /// 
    /// polygon.add((0.0, 2.0)) //Square with vertices (0, 0), (2, 0), (2, 2), and (0, 2)
    /// ```
    pub fn add(&mut self, vertex: (f32, f32))
    {

        self.vertices.push(vertex);

    }

    /// Creates a polygon from the given vertices.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::prelude::*;
    /// 
    /// let vertices = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)];
    /// let rectangle = Polygon::from_vertices((0.0, 0.0), vertices); 
    /// //Rectangle with vertices (0,0), (2,0), (2,1), and (0, 1)
    /// ```
    pub fn from_vertices(position: (f32, f32), vertices: Vec<(f32, f32)>) -> Polygon
    {

        return Polygon { position, vertices };

    }

    /// Attempts to create a convex polygon from the given vertices. It returns the poloygon if
    /// it is convex, and `None` if it is concave. 
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::prelude::*;
    /// 
    /// let vertices = vec![(-1.0, 1.0), (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)];
    /// let square = Polygon::convex_from_vertices((0.0, 0.0), vertices);
    /// 
    /// assert!(square.is_some());
    /// 
    /// let concave_vertices = vec![(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 1.0), (0.0, 2.0)];
    /// let concave_shape = Polygon::convex_from_vertices((0.0, 0.0), concave_vertices);
    /// 
    /// assert!(concave_shape.is_none());
    /// ```
    pub fn convex_from_vertices(position: (f32, f32), vertices: Vec<(f32, f32)>) -> Option<Polygon>
    {

        let polygon = Polygon { position, vertices };

        if polygon.is_convex()
        {

            return Some(polygon);

        }

        return None;

    }

}

impl crate::Shape for Polygon
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

        //We only perform collision check for polygons, not lines or points
        if self.vertices.len() <= 1
        {

            return 0;

        }

        if self.vertices.len() == 2
        {

            return 1;

        }

        return self.vertices.len();

    }


    fn get_axis(&self, index: usize, _target: (f32, f32)) -> (f32, f32)
    {

        if self.vertices.len() <= 1 || index >= self.vertices.len()
        {

            return (0.0, 0.0);

        }

        let vertex = self.vertices[index];
        let next = self.vertices[(index + 1) % self.vertices.len()];

        let side = (next.0 - vertex.0, next.1 - vertex.1);
        let axis = (-side.1, side.0);

        return axis;

    }

    fn project(&self, axis: (f32, f32), _normalize: bool) -> (f32, f32)
    {

        return crate::project(self.position, axis, &self.vertices);

    }

    fn needs_closest(&self, _index: usize) -> bool
    {

        return false;

    }

    fn get_closest(&self, target: (f32, f32)) -> (f32, f32)
    {

        return crate::closest(self.position, target, &self.vertices);

    }

    fn point(&self, index: usize) -> (f32, f32)
    {

        if index >= self.vertices.len()
        {

            return self.position;

        }

        let vertex = self.vertices[index];
        return (self.position.0 + vertex.0, self.position.1 + vertex.1);

    }

}

impl crate::Rotate for Polygon
{

    fn rotate(&mut self, angle: f32)
    {

        let sin = f32::sin(angle);
        let cos = f32::cos(angle);

        self.rotate_sincos(sin, cos);

    }

    fn rotate_sincos(&mut self, sin: f32, cos: f32)
    {

        for v in self.vertices.iter_mut()
        {

            *v = crate::rotate!(sin, cos, v);

        }

    }

}

#[cfg(test)]
mod polygon_tests
{

    use super::*;
    use crate::{float_equal, Shape};

    #[test]
    fn test_add_is_convex()
    {

        let mut square = Polygon::new((0.0, 0.0));
        square.add((0.0, 0.0));
        square.add((1.0, 0.0));
        square.add((1.0, 1.0));
        square.add((0.0, 1.0));

        assert!(square.is_convex());

        //Ensure that it works for polygons with more vertices than necessary
        let mut triangle = Polygon::new((-32.0, 160.0));
        triangle.add((-5.0, 0.0));
        triangle.add((-2.0, -2.0));
        triangle.add((5.0, 0.0));
        triangle.add((0.0, 0.0));

        assert!(triangle.is_convex());

    }

    #[test]
    fn test_add_not_convex()
    {

        let mut four = Polygon::new((0.0, 0.0));
        four.add((0.0, 0.0));
        four.add((1.0, 0.0));
        four.add((0.0, 1.0));
        four.add((1.0, 1.0));

        assert!(!four.is_convex());

        let mut five = Polygon::new((-32.0, 160.0));
        five.add((-5.0, 0.0));
        five.add((-2.0, -2.0));
        five.add((5.0, 0.0));
        five.add((5.0, 10.0));
        five.add((0.0, 0.0));

        assert!(!five.is_convex());

    }

    #[test]
    fn test_from_vertices_is_convex()
    {

        let vertices = vec![(8.0, 0.0), (12.0, 4.0), (6.0, 8.0), (0.0, 4.0), (0.0, 0.0)];
        let pentagon = Polygon::from_vertices((3.0, 2.0), vertices);

        assert!(pentagon.is_convex());

    }

    #[test]
    fn test_from_vertices_not_convex()
    {

        let vertices = vec![(-4.0, 0.0), (-3.0, 2.0), (0.0, 1.0), (3.0, 2.0), (4.0, 0.0), (0.0, 0.0)];
        let concave = Polygon::from_vertices((32.0, 4.56), vertices);

        assert!(!concave.is_convex());

    }

    #[test]
    fn test_convex_from_vertices_is_convex()
    {

        let vertices = vec![(1.0, 1.0), (1.0, 2.0), (0.0, 3.0), (-1.0, 2.0), (-1.0, 1.0), (0.0, 0.0)];
        let hexagon = Polygon::convex_from_vertices((0.0, 0.0), vertices);

        assert!(hexagon.is_some());

    }

    #[test]
    fn test_convex_from_vertices_not_convex()
    {

        let vertices = vec![(1.0, 1.0), (-1.0, 1.0), (0.0, 3.0), (-1.0, 2.0), (2.0, 1.0), (0.0, 0.0)];
        let concave = Polygon::convex_from_vertices((0.0, 0.0), vertices);

        assert!(concave.is_none());

    }

    #[test]
    fn test_num_axes()
    {

        let empty = Polygon::new((1.0, 2.0));
        let point = Polygon::from_vertices((0.0, 3.0), vec![(2.0, 1.0)]);
        let line = Polygon::from_vertices((0.0, 0.0), vec![(0.0, 1.0), (0.0, 0.0)]);

        let vertices = vec![(1.0, 1.0), (1.0, 2.0), (0.0, 3.0), (-1.0, 2.0), (-1.0, 1.0), (0.0, 0.0)];
        let hexagon = Polygon::from_vertices((0.0, 0.0), vertices);

        assert_eq!(empty.num_axes(), 0);
        assert_eq!(point.num_axes(), 0);
        assert_eq!(line.num_axes(), 1);
        assert_eq!(hexagon.num_axes(), 6);


    }

    #[test]
    fn test_get_axis()
    {

        let vertices = vec![(1.0, 1.0), (1.0, 2.0), (0.0, 3.0), (-1.0, 2.0), (-1.0, 1.0), (0.0, 0.0)];
        let hexagon = Polygon::from_vertices((0.0, 0.0), vertices);

        let axis1 = hexagon.get_axis(1, (0.0, 0.0));
        let axis3 = hexagon.get_axis(3, (110.0, 11.0));
        let axis4 = hexagon.get_axis(5, (-1.0, 1.0 / 16.0));

        assert!(float_equal(axis1.0, -1.0));
        assert!(float_equal(axis1.1, -1.0));
        assert!(float_equal(axis3.0, 1.0));
        assert!(float_equal(axis3.1, 0.0));
        assert!(float_equal(axis4.0, -1.0));
        assert!(float_equal(axis4.1, 1.0));

    }

    #[test]
    fn test_project()
    {

        let vertices = vec![(1.0, 1.0), (1.0, 2.0), (0.0, 3.0), (-1.0, 2.0), (-1.0, 1.0), (0.0, 0.0)];
        let hexagon = Polygon::from_vertices((0.0, 0.0), vertices);

        let axis = (1.0, 2.0);
        let projection = hexagon.project(axis, false);

        assert!(float_equal(projection.0, 0.0));
        assert!(float_equal(projection.1, 6.0));

    }

    #[test]
    fn test_needs_closest()
    {

        let square = Polygon::from_vertices((1.0, 2.0), vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);

        assert!(!square.needs_closest(2));

    }

}