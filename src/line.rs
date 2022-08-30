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

/// Checks if the given shape intersects the infinite ray locates at
/// ray_position point in the direction of the vector ray_direction.
pub fn intersects_ray(shape: &(impl Shape + ?Sized), ray_position: (f32, f32), ray_direction: (f32, f32)) -> bool
{
    
    if !intersects_line(shape, ray_position, ray_direction)
    {
        
        return false;
        
    }
    
    let mut normalized = (0.0, 0.0);
    
    let num_axes = shape.num_axes();
    for i in 0..num_axes
    {
        
        let closest = if shape.needs_closest(i)
        {
          
            if normalized == (0.0, 0.0)
            {
             
                let length = f32::sqrt(ray_direction.0 * ray_direction.0 + ray_direction.1 * ray_direction.1);
                normalized = (ray_direction.0 / length, ray_direction.1 / length);
                
            }
            
            let point = shape.point(i);
            let distance = (point.0 * normalized.0) + (point.1 * normalized.1);
            
            if distance < f32::EPSILON
            {
                
                ray_position
                
            }
            else
            {
                
                (ray_position.0 + distance * normalized.0, ray_position.1 + distance * normalized.1)
            
            }
            
        }
        else
        {
            
            (0.0, 0.0)
            
        };
        
        let axis = shape.get_axis(i, closest);
        
        let position = (axis.0 * ray_position.0) + (axis.1 * ray_position.1);
        let projection = shape.project(axis, false);

        let direction = (axis.0 * ray_direction.0) + (axis.1 * ray_direction.1);
        
        if (direction < -f32::EPSILON && position < projection.0 - f32::EPSILON) ||
           (direction > f32::EPSILON && position > projection.1 + f32::EPSILON)
        {
            
            return false;
            
        }
        
    }
    
    return true;
    
}

/// Checks if the given shape intersects the line segment between
/// the two given points.
pub fn intersects_segment(shape: &(impl Shape + ?Sized), line_start: (f32, f32), line_end: (f32, f32)) -> bool
{
    
    let direction = (line_end.0 - line_start.0, line_end.1 - line_start.1);
    
    if !intersects_line(shape, line_start, direction)
    {

        return false;
        
    }
    
    let mut normalized = (0.0, 0.0);
    
    let num_axes = shape.num_axes();
    for i in 0..num_axes
    {
        
        let closest = if shape.needs_closest(i)
        {
          
            if normalized == (0.0, 0.0)
            {
             
                let length = f32::sqrt(direction.0 * direction.0 + direction.1 * direction.1);
                normalized = (direction.0 / length, direction.1 / length);
                
            }
            
            let point = shape.point(i);
            let distance = (point.0 * normalized.0) + (point.1 * normalized.1);
            
            if distance < f32::EPSILON
            {
                
                line_start
                
            }
            else if distance > 1.0 + f32::EPSILON
            {
                
                line_end
                
            }
            else
            {
                
                (line_start.0 + distance * normalized.0, line_start.1 + distance * normalized.1)
            
            }
            
        }
        else
        {
            
            (0.0, 0.0)
            
        };
        
        let axis = shape.get_axis(i, closest);
        
        let start = (axis.0 * line_start.0) + (axis.1 * line_start.1);
        let end = (axis.0 * line_end.0) + (axis.1 * line_end.1);
        let projection = shape.project(axis, false);
        
        if f32::min(start, end) > projection.1 + f32::EPSILON || f32::max(start, end) < projection.0 - f32::EPSILON
        {

            return false;
            
        }
        
    }
    
    return true;
    
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

    #[test]
    fn test_ray_intersection()
    {
        
        let triangle = Polygon::from_vertices((0.0, 0.0), vec![(0.0, 0.0), (-1.0, 1.0), (1.0, 1.0)]);
        let circle = Circle::new((2.0, 2.0), 1.0);
        let capsule = Capsule::new((10.0, 5.0), (-3.0, 2.0), 2.0);

        assert!(intersects_ray(&triangle, (0.0, 0.0), (0.0, 1.0)));
        assert!(intersects_ray(&triangle, (10.0, -9.0), (-1.0, 1.0)));
        assert!(!intersects_ray(&triangle, (10.0, -9.0), (1.0, -1.0)));
        assert!(!intersects_ray(&triangle, (0.0, 2.0), (1.0, 0.0)));

        assert!(intersects_ray(&circle, (2.0, 2.0), (-0.35, 1.15)));
        assert!(intersects_ray(&circle, (2.0, 2.0), (0.35, -1.15)));
        assert!(intersects_ray(&circle, (2.9, 0.0), (0.0, 1.0)));
        assert!(!intersects_ray(&circle, (2.0, 0.0), (0.0, -1.0)));
        assert!(!intersects_ray(&circle, (0.0, 5.6), (1.0, -1.0)));

        assert!(intersects_ray(&capsule, (5.0, 0.0), (1.0, 1.0)));
        assert!(intersects_ray(&capsule, (0.0, 8.9), (1.0, 0.0)));
        assert!(!intersects_ray(&capsule, (0.0, 8.9), (-1.0, 0.0)));
        assert!(!intersects_ray(&capsule, (0.0, 0.9), (-1.0, 0.0)));
        
    }

    #[test]
    fn test_segment_intersection()
    {
        
        let pentagon = Polygon::from_vertices((0.0, 0.0), vec![(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (0.5, 2.0), (-1.0, 1.0)]);
        let circle = Circle::new((0.0, 0.0), 1.0);
        let gram = Parallelogram::new((2.0, 3.0), (1.0, 2.0), (2.0, 1.0));
        
        assert!(intersects_segment(&pentagon, (4.0, 1.0), (-2.0, 1.0)));
        assert!(intersects_segment(&pentagon, (1.1, -0.1), (1.1, 0.4)));
        assert!(!intersects_segment(&pentagon, (1.0, 1.75), (2.0, 1.75)));
        
        assert!(intersects_segment(&circle, (0.5, 0.5), (1.0, 0.5)));
        assert!(intersects_segment(&circle, (-0.9, 0.01), (-1.0, 3.0)));
        assert!(!intersects_segment(&circle, (0.0, 1.6), (1.6, 0.0)));
        
        assert!(intersects_segment(&gram, (2.0, 4.0), (3.0, 3.0)));
        assert!(intersects_segment(&gram, (3.8, 4.0), (3.9, 1.0)));
        assert!(!intersects_segment(&gram, (2.0, 3.1), (3.0, 5.1)));
        
    }

}