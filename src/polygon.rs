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
/// use sepax2d::polygon::Polygon;
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
pub struct Polygon
{

    pub position: (f64, f64),
    pub vertices: Vec::<(f64,f64)>

}

impl Polygon
{

    /// Define a new polygon with no vertices at the given position. 
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::polygon::Polygon;
    /// 
    /// let point = Polygon::new((17.0, 5.0));
    /// ```
    pub fn new(position: (f64, f64)) -> Polygon
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
    /// use sepax2d::polygon::Polygon;
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

                    if min < -f64::EPSILON && max > f64::EPSILON
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

    fn side_projection(&self, axis: (f64, f64), start: usize) -> (f64, f64)
    {

        let mut min = f64::MAX;
        let mut max = f64::MIN;

        //Project all vertices not touching the given side onto the given axis
        for (x, y) in self.vertices.iter().cycle().skip(start + 1).take(self.vertices.len() - 2)
        {

            //Dot product the current position vector to the axis of projection
            let projection = (*x * axis.0) + (*y * axis.1);

            min = f64::min(min, projection);
            max = f64::max(max, projection);

        }

        return (min, max);

    }

    /// Adds a vertex to the given shape, setting it as the last vertex.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::polygon::Polygon;
    /// 
    /// let vertices = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)];
    /// let mut polygon = Polygon::from_vertices((0.0, 0.0), vertices); //Triangle with vertices (0, 0), (2, 0), and (2, 2)
    /// 
    /// polygon.add((0.0, 2.0)) //Square with vertices (0, 0), (2, 0), (2, 2), and (0, 2)
    /// ```
    pub fn add(&mut self, vertex: (f64, f64))
    {

        self.vertices.push(vertex);

    }

    /// Creates a polygon from the given vertices.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::polygon::Polygon;
    /// 
    /// let vertices = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)];
    /// let rectangle = Polygon::from_vertices((0.0, 0.0), vertices); 
    /// //Rectangle with vertices (0,0), (2,0), (2,1), and (0, 1)
    /// ```
    pub fn from_vertices(position: (f64, f64), vertices: Vec<(f64, f64)>) -> Polygon
    {

        return Polygon { position, vertices };

    }

    /// Attempts to create a convex polygon from the given vertices. It returns the poloygon if
    /// it is convex, and `None` if it is concave. 
    /// 
    /// # Examples
    /// 
    /// ```
    /// use sepax2d::polygon::Polygon;
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
    pub fn convex_from_vertices(position: (f64, f64), vertices: Vec<(f64, f64)>) -> Option<Polygon>
    {

        let polygon = Polygon { position, vertices };

        if polygon.is_convex()
        {

            return Some(polygon);

        }

        return None;

    }

}

impl super::Shape for Polygon
{

    fn position(&self) -> (f64, f64)
    {

        return self.position;

    }

    fn num_axes(&self) -> usize
    {

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


    fn get_axis(&self, index: usize) -> (f64, f64)
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

    fn project(&self, axis: (f64, f64)) -> (f64, f64)
    {

        let mut min = f64::MAX;
        let mut max = f64::MIN;

        for (x, y) in self.vertices.iter()
        {

            let position = (self.position.0 + *x, self.position.1 + *y);

            let projection = (position.0 * axis.0) + (position.1 * axis.1);

            min = f64::min(min, projection);
            max = f64::max(max, projection);

        }

        return (min, max);

    }

}

#[cfg(test)]
mod polygon_tests
{

    use super::*;

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

}