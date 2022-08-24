use crate::Shape;

/// Checks if the given shape intersects the infinite line located at 
/// line_position pointing in the direction of the vector line_direction. 
pub fn intersects_line(shape: &(impl Shape + ?Sized), line_position: (f32, f32), line_direction: (f32, f32)) -> bool
{

    let axis = (-line_direction.1, line_direction.0);
    let location = (axis.0 * line_position.0) + (axis.1 * line_position.1);

    let projection = shape.project(axis, false);

    return projection.0 < location - f32::EPSILON && projection.1 + f32::EPSILON > location;

}

#[cfg(test)]
mod line_tests
{

    use crate::prelude::*;

    #[test]
    fn test_line_intersection()
    {

        let triangle = Polygon::from_vertices((0.0, 0.0), vec![(0.0, 0.0), (-1.0, 1.0), (1.0, 1.0)]);
        let circle = Circle::new((2.0, 2.0), 1.0);
        let capsule = Capsule::new((10.0, 5.0), (-3.0, 2.0), 2.0);

        assert!(intersects_line(&triangle, (0.0, 0.0), (0.0, 1.0)));
        assert!(intersects_line(&triangle, (10.0, -9.0), (-1.0, 1.0)));
        assert!(!intersects_line(&triangle, (0.0, 2.0), (1.0, 0.0)));

        assert!(intersects_line(&circle, (2.0, 2.0), (-0.35, 1.15)));
        assert!(intersects_line(&circle, (2.9, 0.0), (0.0, 1.0)));
        assert!(!intersects_line(&circle, (0.0, 5.6), (1.0, -1.0)));

        assert!(intersects_line(&capsule, (5.0, 0.0), (1.0, 1.0)));
        assert!(intersects_line(&capsule, (0.0, 8.9), (1.0, 0.0)));
        assert!(!intersects_line(&capsule, (0.0, 0.9), (-1.0, 0.0)));

    }

}